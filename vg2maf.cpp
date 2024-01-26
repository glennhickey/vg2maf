#include "vg2maf.hpp"

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <fstream>
#include <memory>
#include <unistd.h>
#include <getopt.h>
#include <omp.h>

#include "handlegraph/path_handle_graph.hpp"
#include "bdsg/packed_graph.hpp"
#include "bdsg/hash_graph.hpp"
#include "bdsg/snarl_distance_index.hpp"
#include "bdsg/overlays/overlay_helper.hpp"
#include "vg/io/vpkg.hpp"
#include "stream_index.hpp"
extern "C" {
#include "taf.h"
#include "sonLib.h"
}

// number of nodes to scan before sending to parallel batch
static const int64_t node_buffer_size = 250000;
// number of bases for each gam index query
static const int64_t gam_idx_query_bp = 10000;
// biggest hole allowable in id range passed to gam index
static const int64_t gam_idx_max_gap = 2;
// maximum read depth (to prevent memory issues)
static const int64_t gam_max_depth = 1000;

//#define debug

using namespace std;
using namespace handlegraph;
using namespace bdsg;

// copied from vg/src/path.cpp
static void reverse_complement_mapping_in_place(vg::Mapping* m, const function<int64_t(nid_t)>& node_length) {
            
    vg::Position* pos = m->mutable_position();
    pos->set_is_reverse(!pos->is_reverse());
    int length = 0;
    for (const auto& edit : m->edit()) {
        length += edit.from_length();
    }
    pos->set_offset(node_length(pos->node_id()) - pos->offset() - length);
    
    size_t swap_size = m->edit_size() / 2;
    for (size_t i = 0, j = m->edit_size() - 1; i < swap_size; i++, j--) {
        vg::Edit* e1 = m->mutable_edit(i);
        vg::Edit* e2 = m->mutable_edit(j);
        
        int64_t from_length_tmp = e1->from_length();
        int64_t to_length_tmp = e1->to_length();
        string sequence_tmp = e1->sequence();
        
        e1->set_from_length(e2->from_length());
        e1->set_to_length(e2->to_length());
        e1->set_sequence(reverse_complement(e2->sequence()));
        
        e2->set_from_length(from_length_tmp);
        e2->set_to_length(to_length_tmp);
        e2->set_sequence(reverse_complement(sequence_tmp));
    }
    
    
    if (m->edit_size() % 2) {
        vg::Edit* e = m->mutable_edit(swap_size);
        reverse_complement_in_place(*e->mutable_sequence());
    }
}

unordered_map<int64_t, unordered_map<int64_t, string>> get_insertion_index(const vector<vg::Mapping>& mappings) {

    unordered_map<int64_t, unordered_map<int64_t, string>>idx;
        
    for (int64_t row = 0; row < mappings.size(); ++row) {
        const vg::Mapping& mapping = mappings.at(row);
        int64_t offset = mapping.position().offset();
        for (int64_t i = 0; i < mapping.edit_size(); ++i) {
            const vg::Edit& edit = mapping.edit(i);
            if (edit.from_length() < edit.to_length()) {
                if (mapping.position().is_reverse()) {
                    cerr << "skipping backward mapping for insertion index!" << endl;
                    continue;
                }
                // note: offsetting on from_length, since we only want to include the "Dangling" inserted bit
                idx[offset + edit.from_length()][row] =  edit.sequence().substr(edit.from_length());   
            }                        
            offset += edit.from_length();
        }        
    }

    return idx;
}

unordered_map<int64_t, unordered_map<int64_t, string>> align_insertion_index(const unordered_map<int64_t, unordered_map<int64_t, string>>& in_idx) {

    // todo: replace this with abPOA.
    // for now, we just leave all insertions unaligned.
    // note: that taffy norm will *not* realign insertions in the middle of blocks.
    //       could explore chopping (instead of abPOA) but it kind of sounds like more hassle right now

    unordered_map<int64_t, unordered_map<int64_t, string>> out_idx;
        
    for (const pair<int64_t, unordered_map<int64_t, string>>& in_elem : in_idx) {
        const unordered_map<int64_t, string>& in_map = in_elem.second;
        unordered_map<int64_t, string>& out_map = out_idx[in_elem.first];
        int64_t width = 0;
        vector<int64_t> rows;
        // compute the number of columns (just sum of all lengths since we're not aligning)
        for (const auto& ie : in_map) {
            width += ie.second.length();
            rows.push_back(ie.first);
        }
        std::sort(rows.begin(), rows.end());
        // allocate the output alignment rows
        for (const auto& ie : in_map) {
            out_map[ie.first].resize(width);
        }
        // write the unaligned output alignment
        int64_t cur_row = 0;
        int64_t cur_offset = 0;
        // do it column by column
        for (int64_t col = 0; col < width; ++col) {
            bool shift_row = false;
            for (int64_t r_idx = 0; r_idx < rows.size(); ++r_idx) {
                int64_t i_row = rows[r_idx];
                const string& i_string = in_map.at(i_row);
                // write the sequence for each row
                if (r_idx == cur_row) {
                    out_map[i_row][col] = i_string[cur_offset];
                    ++cur_offset;
                    if (cur_offset == i_string.length()) {
                        shift_row = true;
                    }
                } else {
                    // and the rest is unaligned
                    out_map[i_row][col] = '-';
                }
            }
            if (shift_row) {
                cur_offset = 0;
                ++cur_row;
                shift_row = false;
            }
        }
        assert(cur_offset == 0 && cur_row == rows.size());
    }

#ifdef debug
    for (auto xx : out_idx) {
        for (auto yy : xx.second) {
            cerr  << "INS[" << xx.first << "][" << yy.first <<"]=" << yy.second << endl;
        }
    }
#endif

    return out_idx;
}

vector<handle_t> get_ref_traversal(PathPositionHandleGraph& graph, path_handle_t ref_path_handle, handle_t start_handle, handle_t end_handle) {

    if (start_handle == end_handle) {
        return {start_handle};
    }
    
    vector<handle_t> ref_trav;
    vector<step_handle_t> start_step_handles = graph.steps_of_handle(start_handle);
    for (step_handle_t start_step : start_step_handles) {
        if (graph.get_path_handle_of_step(start_step) == ref_path_handle) {
            bool forward = false;
            // look for next id on path going into snarl
            // todo: is there a crazy case where we need strand (thank vg call code deals with all this)
            if (graph.has_next_step(start_step)) {
                nid_t next_id = graph.get_id(graph.get_handle_of_step(graph.get_next_step(start_step)));
                graph.follow_edges(start_handle, false, [&](handle_t follow_handle) {
                    if (graph.get_id(follow_handle) == next_id) {
                        forward = true;
                    }
                });
            }
            if (forward == false) {
                bool found = false;
                if (graph.has_previous_step(start_step)) {
                    nid_t prev_id = graph.get_id(graph.get_handle_of_step(graph.get_previous_step(start_step)));
                    graph.follow_edges(start_handle, false, [&](handle_t follow_handle) {
                        if (graph.get_id(follow_handle) == prev_id) {
                            found = true;
                    }
                    });
                }
                assert(found);
            }
            ref_trav.clear();
            if (forward) {
                step_handle_t end_of_path = graph.path_end(ref_path_handle);
                for (step_handle_t cur_step = start_step; cur_step != end_of_path; cur_step = graph.get_next_step(cur_step)) {
                    handle_t cur_handle = graph.get_handle_of_step(cur_step);
                    ref_trav.push_back(cur_handle);
                    if (graph.get_id(cur_handle) == graph.get_id(end_handle)) {
                        return ref_trav;
                    }
                }
            } else {
                step_handle_t end_of_path = graph.path_front_end(ref_path_handle);
                for (step_handle_t cur_step = start_step; cur_step != end_of_path; cur_step = graph.get_previous_step(cur_step)) {
                    handle_t cur_handle = graph.get_handle_of_step(cur_step);
                    ref_trav.push_back(graph.flip(cur_handle));
                    if (graph.get_id(cur_handle) == graph.get_id(end_handle)) {
                        return ref_trav;
                    }                
                }
            }
        }
    }
    cerr << "[vg2maf] warning: Unable to find reference path (" << graph.get_path_name(ref_path_handle) << ") through snarl "
         << graph.get_id(start_handle) << ":" << graph.get_is_reverse(start_handle) << " -> "
         << graph.get_id(end_handle) << ":" << graph.get_is_reverse(end_handle)
         << ". Will export anyway but blocks may be badly unsorted" << endl;
    return vector<handle_t>();
}

Alignment* convert_node(PathPositionHandleGraph& graph, const vector<vg::Alignment>& gam_alignments,
                        handle_t handle, path_handle_t ref_path_handle) {

    // we don't care about the chain orientation, everything below is based on the underlying node being forward
    if (graph.get_is_reverse(handle)) {
        handle = graph.flip(handle);
    }
    
    vector<step_handle_t> steps = graph.steps_of_handle(handle);
    if (steps.empty()) {
        //cerr << "[vg2maf] warning: Skipping node " << graph.get_id(handle) << " because there are no paths on it" << endl;
        return NULL;
    }
    
    Alignment* alignment = (Alignment*)st_calloc(1, sizeof(Alignment));
    alignment->row_number = steps.size();
    alignment->column_number = graph.get_length(handle);

    // convert each step to row
    vector<Alignment_Row*> rows;

    // index for stitching in insertion alignments, which require extra columns to be inserted
    // this is a basically 
    unordered_map<int64_t, unordered_map<int64_t, string>> ins_alignments;

    nid_t node_id = graph.get_id(handle);
    string node_sequence = graph.get_sequence(handle);

    unordered_map<Alignment_Row*, string> base_qualities;

    // convert the gam indxes to rows
    if (!gam_alignments.empty()) {
        // query the index for our node [todo: is this too simplisitic due to repeated queries?]
        vector<vg::Mapping> mappings;
        vector<string> names;
        vector<int64_t> start_positions;
        vector<int64_t> sequence_lengths;
        vector<bool> mapping_reversed;
        vector<int64_t> alignment_index;
#ifdef debug
        cerr << "doing gam index query on node " << node_id << endl;
#endif
        for (int64_t ai = 0; ai < gam_alignments.size(); ++ai) {
            const vg::Alignment& aln = gam_alignments[ai];
            // this is the position w.r.t the read alignment
            int64_t pos = 0;
            for (int64_t i = 0; i < aln.path().mapping_size(); ++i) {
                const vg::Mapping& mapping = aln.path().mapping(i);
                if (mapping.position().node_id() == node_id) {
                    mappings.push_back(mapping);
                    alignment_index.push_back(ai);
                    mapping_reversed.push_back(mapping.position().is_reverse());
                    // easier just to deal with forward mapping from here down
                    // we remember it was reversed in the mapping_reversed vector for the strand field
                    if (mapping.position().is_reverse()) {
                        reverse_complement_mapping_in_place(&mappings.back(), [&graph](nid_t mnid) {
                            return graph.get_length(graph.get_handle(mnid));
                        });
                    }
                    string name = aln.name();
                    if (name.empty()) {
                        name = "aln"; 
                    }
                    names.push_back(name);
                    start_positions.push_back(pos);
                }
                for (int64_t j = 0; j < mapping.edit_size(); ++j) {
                    pos += mapping.edit(j).to_length();
                }
            }
            int64_t n_mappings = mappings.size() - sequence_lengths.size();
            // manually count it up, to support case where sequence string not in alignment (gaf?)
            for (int64_t i = 0; i < n_mappings; ++i) {
                if (!aln.sequence().empty()) {
                    //assert(aln.sequence().length() == pos);
                }
                sequence_lengths.push_back(pos);
            }
        }

        // collect all the inserted sequences, indexed by column
        auto ins_idx = get_insertion_index(mappings);

        // make a dummy alignment [todo: switch to actual alignment]
        ins_alignments = align_insertion_index(ins_idx);

        // add enough alignment columns for all insertions
        for (const auto& elem : ins_alignments) {
            alignment->column_number += elem.second.begin()->second.length();
        }

#ifdef debug
        if (mappings.size() > 0) {
            cerr << "adding " << mappings.size() << " mappings for node " << node_id << endl;
        }
#endif
        // add in insertion gaps
        auto insert_insertion_gaps = [&](int64_t node_offset, int64_t row,
                                         int64_t& column, char* bases) {
            if (!ins_alignments.at(node_offset).count(row)) {
                int64_t gaps = ins_alignments.at(node_offset).begin()->second.length();
                for (int64_t i = 0; i < gaps; ++i) {
                    bases[column++] = '-';
                }
            }
        };

        // copy our mappings into the alignment rows
        for (int64_t i = 0; i < mappings.size(); ++i) {
            const vg::Alignment& aln = gam_alignments[alignment_index[i]];
            vg::Mapping& mapping = mappings[i];
            Alignment_Row* row = (Alignment_Row*)st_calloc(1, sizeof(Alignment_Row));
            ++alignment->row_number;
            row->sequence_name = stString_copy(names[i].c_str());
            row->length = 0;
            row->sequence_length = sequence_lengths[i];
            row->start = start_positions[i];
            row->strand = mapping_reversed[i] ? 0 : 1;
            row->bases = (char*)st_calloc(alignment->column_number + 1, sizeof(char));
            string* row_qualities = nullptr;
            if (!aln.quality().empty()) {
                // default to maximum phred character
                base_qualities[row].resize(alignment->column_number, (char)0x7e);
                row_qualities = &base_qualities[row];
            }
            // add the opening gaps
            int64_t col = 0;
            int64_t node_offset = 0;
            for (; node_offset < mappings[i].position().offset(); ++col, ++node_offset) {
                if (ins_alignments.count(node_offset)) {
                    insert_insertion_gaps(node_offset, i, col, row->bases);
                }
                row->bases[col] = '-';
            }
            
            // add in the sequence
            for (int64_t j = 0; j < mappings[i].edit_size(); ++j) {
                const vg::Edit& edit = mappings[i].edit(j);
                if (edit.from_length() == edit.to_length()) {
                    //match
                    for (int64_t k = 0; k < edit.from_length(); ++k) {
                        if (ins_alignments.count(node_offset)) {
                            insert_insertion_gaps(node_offset, i, col, row->bases);
                        }
                        if (!edit.sequence().empty()) {
                            row->bases[col] = edit.sequence()[k];                            
                        } else {
                            row->bases[col] = node_sequence[node_offset];
                        }
                        if (!aln.quality().empty()) {
                            (*row_qualities)[col] = (char)aln.quality()[row->length];
                        }
                        ++col;
                        ++node_offset;
                        ++row->length;
                    }
                } else if (edit.to_length() == 0 && edit.from_length() > 0) {
                    // delete
                    for (int64_t k = 0; k < edit.from_length(); ++k) {
                        if (ins_alignments.count(node_offset)) {
                            insert_insertion_gaps(node_offset, i, col, row->bases);
                        }                        
                        row->bases[col++] = '-';
                        ++node_offset;
                    }
                } else {
                    // insert
                    assert(edit.from_length() < edit.to_length());
                    // add the common part [todo: this should probably not be aligned automaticall going forward]
                    for (int64_t k = 0; k < edit.from_length(); ++k) {
                        if (!aln.quality().empty()) {
                            (*row_qualities)[col] = (char)aln.quality()[row->length];
                        }
                        row->bases[col++] = edit.sequence()[k];
                        ++node_offset;
                        ++row->length;
                    }
                    // look up the insertion from the index
                    const unordered_map<int64_t, string>& row_alignments = ins_alignments.at(node_offset);
                    // find the row
                    const string& row_string = row_alignments.at(i);
                    for (int64_t k = 0; k < row_string.length(); ++k) {
                        if (row_string[k] != '-') {
                            if (!aln.quality().empty()) {
                                (*row_qualities)[col] = (char)aln.quality()[row->length];
                            }
                            ++row->length;
                        }
                        row->bases[col++] = row_string[k];
                    }
                }
            }

            // add the closing gaps
            for (; node_offset <= node_sequence.length(); ++node_offset) {
                if (ins_alignments.count(node_offset)) {
                    insert_insertion_gaps(node_offset, i, col, row->bases);
                }
                if (node_offset < node_sequence.length()) {
                    row->bases[col++] = '-';
                }
            }
            assert(col == alignment->column_number);
            row->bases[col] = '\0';
            
            if (row->strand == 0) {
                row->start = row->sequence_length - row->start - row->length;
            }

            rows.push_back(row);
        }
    }

#ifdef debug
    cerr << "alignment rows after adding mappings: " << rows.size() << endl;
#endif

    for (step_handle_t& step_handle : steps) {
        Alignment_Row* row = (Alignment_Row*)st_calloc(1, sizeof(Alignment_Row));
        path_handle_t step_path_handle = graph.get_path_handle_of_step(step_handle);
        row->sequence_name = stString_copy(graph.get_path_name(step_path_handle).c_str());
        row->start = graph.get_position_of_step(step_handle);
        row->length = node_sequence.length();
        row->sequence_length = graph.get_path_length(step_path_handle);
        // todo: check this strand logic
        handle_t handle_of_step = graph.get_handle_of_step(step_handle);
        row->strand = graph.get_is_reverse(handle_of_step) ? 0 : 1;
        if (row->strand == 0) {
            // if reverse strand, MAF wants coordinates from end of path. 
            row->start = row->sequence_length - row->start - row->length;
        }
        int64_t gaps = 0;
        for (const pair<int64_t, unordered_map<int64_t, string>>& ie : ins_alignments) {
            gaps += ie.second.begin()->second.length();
        }
        row->bases = (char*)st_calloc(node_sequence.length() + gaps + 1, sizeof(char));
        // calloc should do this but just in case
        row->bases[node_sequence.length() + gaps] = '\0';
        int64_t maf_col = 0;
        // copy the node sequence in base by base
        for (int64_t col = 0; col <= node_sequence.length(); ++col) {
            // add insertion gaps if the column is in the insertion index
            if (ins_alignments.count(col)) {
                int64_t gaps = ins_alignments[col].begin()->second.length();
                for (int64_t k = 0; k < gaps; ++k) {
                    row->bases[maf_col++] = '-';
                }
                
            }
            if (col < node_sequence.length()) {
                row->bases[maf_col++] = node_sequence[col];
            }
        }
        assert(maf_col == alignment->column_number);
        rows.push_back(row);
    }

#ifdef debug
    cerr << "alignment rows after adding paths: " << rows.size() << endl;
#endif

    // sort the rows
    std::sort(rows.begin(), rows.end(), [&](const Alignment_Row* row1, const Alignment_Row* row2) {
        int cmp = strcmp(row1->sequence_name, row2->sequence_name);
        return cmp < 0 || (cmp == 0 && row1->start < row2->start);
    });

    // put them in the alignment, starting with ref path
    alignment->row = NULL;    
    string ref_path_name = graph.get_path_name(ref_path_handle);
    for (Alignment_Row* row : rows) {
        if (strcmp(row->sequence_name, ref_path_name.c_str()) == 0) {
            alignment->row = row;
            break;
        }
    }

    assert(alignment->row_number == rows.size());
    Alignment_Row* cur_row = alignment->row;
    for (Alignment_Row* row : rows) {
        if (cur_row == NULL) {
            alignment->row = row;
            cur_row = row;
        } else if (row != alignment->row) {
            cur_row->n_row = row;
            cur_row = row;            
        }
    }

    // allocate the tags array which as it's required my taf
    alignment->column_tags = (Tag**)st_calloc(alignment->column_number, sizeof(Tag*));

    // add in the base qualities
    if (!base_qualities.empty()) {
        // transpose our row qualities into column qualities
        // todo: we can do columns from the get-go but would need to refactor
        // so that rows don't get sorted after being built
        vector<string> column_quality_strings(alignment->column_number);
        for (int64_t i = 0; i < alignment->column_number; ++i) {
            column_quality_strings[i].resize(alignment->row_number);
        }
        string empty_quality(alignment->column_number, char(0x7e)); // default to *max* ascii phred quality
        cur_row = alignment->row;
        for (int64_t i = 0; cur_row; cur_row = cur_row->n_row, ++i) {
            string& row_quals = base_qualities.count(cur_row) ? base_qualities[cur_row] : empty_quality;
            assert(row_quals.length() == alignment->column_number);
            for (int64_t j = 0; j < row_quals.length(); ++j) {
                column_quality_strings[j][i] = (char)phred_byte_to_ascii(row_quals[j]);
            }                
        }
        // copy the column qualities into tags
        for (int64_t i = 0; i < alignment->column_number; ++i) {
            alignment->column_tags[i] = tag_construct((char*)TAF_BASE_QUALITY_TAG_KEY,
                                                      (char*)column_quality_strings[i].c_str(), NULL);
        }        
    }

    return alignment; // caller must run alignment_destruct()    
}

void convert_node_range(PathPositionHandleGraph& graph, GAMInfo* gam_info, const vector<handle_t>& node_buffer,
                        const vector<int64_t>& sorted_nodes, int64_t range_start, int64_t range_end,
                        path_handle_t ref_path_handle, vector<Alignment*>& out_alignment_buffer) {

    // perform range query on gam index
    nid_t first_id = graph.get_id(node_buffer[sorted_nodes[range_start]]);
    nid_t last_id = graph.get_id(node_buffer[sorted_nodes[range_end - 1]]);
    vector<vg::Alignment> alignments;
    int64_t total_aln_length = 0;
    bool warned = false;
    if (gam_info != nullptr) {
        gam_info->index.find(gam_info->cursor, first_id, last_id, [&](const vg::Alignment& aln) {
            if (total_aln_length / gam_idx_query_bp <= gam_max_depth) {
                alignments.push_back(aln);                
            } else if (warned) {
                cerr << "[vg2maf] warning: dropping alignments found in node range " << range_start <<"-" << range_end
                     << " to avoid making a giant buffer (coverage limit " << gam_max_depth << " exceeded)." << endl;
                warned = true;
            }
            total_aln_length += aln.sequence().length();
        });
    }

    // convert nodes individually, letting them sort out which alignments are relevant
    // (this should still be more efficient that doing individual gam index queries)
    for (int64_t i = range_start; i < range_end; ++i) {
        int64_t buffer_index = sorted_nodes[i];
        const handle_t& handle = node_buffer[buffer_index];
        out_alignment_buffer[buffer_index] = convert_node(graph, alignments, handle, ref_path_handle);
    }
}

void traverse_snarl(PathPositionHandleGraph& graph, SnarlDistanceIndex& distance_index,
                    net_handle_t snarl, path_handle_t ref_path_handle,
                    bool ref_path_reversed, vector<handle_t>& out_handles) {

    net_handle_t start_bound = distance_index.get_bound(snarl, false, true);
    net_handle_t end_bound = distance_index.get_bound(snarl, true, false);

    handle_t start_handle = distance_index.get_handle(start_bound, &graph);
    handle_t end_handle = distance_index.get_handle(end_bound, &graph);
    if (ref_path_reversed) {
        std::swap(start_handle, end_handle);
        start_handle = graph.flip(start_handle);
        end_handle = graph.flip(end_handle);
    }

#ifdef debug
    cerr << "snarl goes from " << graph.get_id(start_handle) << ":" << graph.get_is_reverse(start_handle) << " to "
         << graph.get_id(end_handle) << ":" << graph.get_is_reverse(end_handle) << endl;
#endif    

    // quick hack to try to enforce ordering on reference path
    // (does not apply of no reference traversal)
    vector<handle_t> ref_path = get_ref_traversal(graph, ref_path_handle, start_handle, end_handle);
    unordered_map<handle_t, int64_t> ref_order;
    for (int64_t i = 0; i < ref_path.size(); ++i) {
        ref_order[ref_path[i]] = i;
    }
    int64_t cur_ref_pos = 1; // we dont search start, so set to 1 and not 0
    
    // use start-to-end bfs search to order the blocks
    deque<handle_t> bfs_queue;
    unordered_set<handle_t> visited = {start_handle, end_handle, graph.flip(start_handle)};    
    graph.follow_edges(start_handle, false, [&](handle_t other_handle) {
        bfs_queue.push_back(other_handle);                
    });

    while (!bfs_queue.empty()) {
        handle_t handle = bfs_queue.front();
        bfs_queue.pop_front();
        if (!visited.count(handle)) {
            if (ref_order.count(handle) && ref_order.at(handle) != cur_ref_pos) {
                // hack to enforce ordering on ref path
                bfs_queue.push_back(handle);
            } else {
                out_handles.push_back(handle);
                visited.insert(handle);
                graph.follow_edges(handle, false, [&](handle_t other_handle) {
                    bfs_queue.push_back(other_handle);
                });
                graph.follow_edges(handle, true, [&](handle_t other_handle) {
                    bfs_queue.push_back(other_handle);
                });
                if (ref_order.count(handle)) {
                    assert(cur_ref_pos == ref_order.at(handle));
                    ++cur_ref_pos;
                }
            }
        }
    }
}

void convert_chain(PathPositionHandleGraph& graph, SnarlDistanceIndex& distance_index, vector<GAMInfo*>& gam_info,
                   net_handle_t chain, const string& ref_path, bool progress, const pair<int64_t, int64_t>& chain_idx,
                   bool taf_output, LW* output) {

    net_handle_t start_bound = distance_index.get_bound(chain, false, true);
    net_handle_t end_bound = distance_index.get_bound(chain, true, false);

    handle_t start_handle = distance_index.get_handle(start_bound, &graph);
    handle_t end_handle = distance_index.get_handle(end_bound, &graph);

#ifdef debug
    cerr << "chain goes from " << graph.get_id(start_handle) << ":" << graph.get_is_reverse(start_handle) << " to "
         << graph.get_id(end_handle) << ":" << graph.get_is_reverse(end_handle) << endl;
#endif

    unordered_map<path_handle_t, int64_t> start_ref_paths;
    unordered_map<path_handle_t, int64_t> end_ref_paths;    
    graph.for_each_step_on_handle(start_handle, [&](step_handle_t step_handle) {
        path_handle_t step_path_handle = graph.get_path_handle_of_step(step_handle);
        string path_name = graph.get_path_name(step_path_handle);
        if (path_name.compare(0, ref_path.length(), ref_path) == 0) {
            if (start_ref_paths.count(step_path_handle)) {
                cerr << "[vg2maf] warning: multiple steps of reference path " << graph.get_path_name(step_path_handle)
                     << " found on chain start node " << graph.get_id(start_handle)
                     << ". MAF blocks and rows may be out of order in output" << endl;
            }
            start_ref_paths[step_path_handle] = graph.get_position_of_step(step_handle);
        }
    });
    if (!start_ref_paths.empty()) {
        graph.for_each_step_on_handle(end_handle, [&](step_handle_t step_handle) {
            path_handle_t step_path_handle = graph.get_path_handle_of_step(step_handle);
            if (start_ref_paths.count(step_path_handle)) {
                if (end_ref_paths.count(step_path_handle)) {
                    cerr << "[vg2maf] warning: multiple steps of reference path " << graph.get_path_name(step_path_handle)
                         << " found on chain end node " << graph.get_id(start_handle)
                         << ". MAF blocks and rows may be out of order in output" << endl;
                }
                end_ref_paths[step_path_handle] = graph.get_position_of_step(step_handle);
            }
        });                
    }
        
    if (end_ref_paths.empty()) {
        cerr <<"[vg2maf] warning: Skipping chain because no reference path found through chain from " 
             << graph.get_id(start_handle) << ":" << graph.get_is_reverse(start_handle) << " to "
             << graph.get_id(end_handle) << ":" << graph.get_is_reverse(end_handle) << endl;
        return;
    }

    if (end_ref_paths.size() > 1) {
        cerr << "[vg2maf] warning: " << end_ref_paths.size() << " reference paths found through chain from "
             << graph.get_id(start_handle) << ":" << graph.get_is_reverse(start_handle) << " to "
             << graph.get_id(end_handle) << ":" << graph.get_is_reverse(end_handle) << ". Just choosing first" << endl;
        cerr << "the paths are "; for (const auto xx : end_ref_paths) cerr << graph.get_path_name(xx.first) << ", "; cerr << endl;
    }

    path_handle_t ref_path_handle = end_ref_paths.begin()->first;

    if (progress) {
        cerr << "[vg2maf]: Converting chain " << chain_idx.first << " / " << chain_idx.second
             << " on " << graph.get_path_name(ref_path_handle) 
             << " from " << graph.get_id(start_handle) << ":" << graph.get_is_reverse(start_handle)
             << " to " << graph.get_id(end_handle) << ":" << graph.get_is_reverse(end_handle)
             << endl;
    }

    //todo: this may fall down with cyclic ref
    bool ref_path_reversed = end_ref_paths.begin()->second < start_ref_paths.at(ref_path_handle);   
    if (ref_path_reversed) {
        // todo: why does this seem to have no effect on iteration below.
        chain = distance_index.flip(chain);
    }

    vector<net_handle_t> chain_childs;
    distance_index.for_each_child(chain, [&](net_handle_t net_handle) {
        chain_childs.push_back(net_handle);
    });
    if (ref_path_reversed) {
        std::reverse(chain_childs.begin(), chain_childs.end());
    }

    vector<handle_t> node_buffer;
    vector<Alignment*> alignment_buffer;
    Alignment* prev_alignment = nullptr;
    
    // convert the chain, one node/snarl at a time
    for (int64_t i = 0; i < chain_childs.size(); ++i) {
        net_handle_t net_handle = chain_childs[i];
        // queue up all child nodes in order we want to convert
        if (distance_index.is_node(net_handle)) {
            node_buffer.push_back(distance_index.get_handle(net_handle, &graph));
        } else if (distance_index.is_snarl(net_handle)) {
            traverse_snarl(graph, distance_index, net_handle, ref_path_handle, ref_path_reversed, node_buffer);
        } else if (distance_index.is_chain(net_handle)) {
            cerr << "TODO: CHILD CHAIN" << endl;
        } else {
            assert(false);
        }

        if (node_buffer.size() >= node_buffer_size || i == chain_childs.size() - 1) {
            // sort the nodes into ranges
            // todo: this logic only required for gam index acces, but think it's fast
            // enough that it doesn't need to be conditional on -g
            // sorted_nodes stores offset in node_buffer
            vector<int64_t> sorted_nodes(node_buffer.size());
            std::iota(sorted_nodes.begin(), sorted_nodes.end(), 0);
            std::sort(sorted_nodes.begin(), sorted_nodes.end(), [&graph,&node_buffer](int64_t x1, int64_t x2) {
                return graph.get_id(node_buffer[x1]) < graph.get_id(node_buffer[x2]);
            });
            // ranges stores offsets in sorted_nodes
            vector<int64_t> ranges;
            int64_t range_size = numeric_limits<int64_t>::max();
            nid_t prev_id = 0;
            for (int64_t j = 0; j < sorted_nodes.size(); ++j) {
                nid_t cur_id = graph.get_id(node_buffer[sorted_nodes[j]]);
                if (range_size > gam_idx_query_bp || (j > 0 && cur_id - prev_id > gam_idx_max_gap)) {
                    ranges.push_back(j);
                    range_size = 0;
                }
                range_size += graph.get_length(node_buffer[sorted_nodes[j]]);
                prev_id = cur_id;
            }

            alignment_buffer.resize(node_buffer.size());
#pragma omp parallel for schedule(dynamic, 1)
            for (int64_t j = 0; j < ranges.size(); ++j) {
                int tid = omp_get_thread_num();
                // these are offsets (open end) in sorted_nodes that define our batch
                int64_t range_start = ranges[j];
                int64_t range_end = j < ranges.size() - 1 ? ranges[j+1] : sorted_nodes.size();
                convert_node_range(graph, gam_info[tid], node_buffer, sorted_nodes, range_start, range_end,
                                   ref_path_handle, alignment_buffer);
            }

            // write them in series
            for (int64_t j = 0; j < alignment_buffer.size(); ++j) {
                if (alignment_buffer[j]) {
                    if (taf_output) {
                        taf_write_block(prev_alignment, alignment_buffer[j], false, 10000, output);
                    } else {
                        maf_write_block(alignment_buffer[j], output);
                    }
                    if (prev_alignment) {
                        alignment_destruct(prev_alignment, true);
                    }
                    prev_alignment = alignment_buffer[j];
                }
            }
            
            node_buffer.clear();
        }
    }
    assert(node_buffer.empty());
    if (prev_alignment) {
        alignment_destruct(prev_alignment, true);
    }
}

unique_ptr<PathHandleGraph> load_graph(istream& graph_stream) {

    char magic_bytes[4];
    graph_stream.read(magic_bytes, 4);
    uint32_t magic_number = ntohl(*((uint32_t*) magic_bytes));
    graph_stream.clear();
    graph_stream.seekg(0, ios::beg);

    PathHandleGraph* graph;
    if (magic_number == PackedGraph().get_magic_number()) {
        graph = new PackedGraph();
    } else if (magic_number == HashGraph().get_magic_number()) {
        graph = new HashGraph();
    }  else {
        cerr << "Unable to parse input graph with magic number " << magic_number << endl;
        exit(1);
    }
    dynamic_cast<SerializableHandleGraph*>(graph)->deserialize(graph_stream);

    return unique_ptr<PathHandleGraph>(graph);
}
