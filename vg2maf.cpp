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

//#define debug

using namespace std;
using namespace handlegraph;
using namespace bdsg;


// from hal2vg/clip-vg.cpp
static unique_ptr<PathHandleGraph> load_graph(istream& graph_stream);

// start with really barebones placeholde implementation
struct GAMInfo {
    vg::GAMIndex index;
    vg::GAMIndex::cursor_t cursor;
};

// from vg/src/path.cpp
static void reverse_complement_mapping_in_place(vg::Mapping* m, const function<int64_t(nid_t)>& node_length);
static unordered_map<int64_t, unordered_map<int64_t, string>> get_insertion_index(const vector<vg::Mapping>& mappings);
static unordered_map<int64_t, unordered_map<int64_t, string>> align_insertion_index(const unordered_map<int64_t, unordered_map<int64_t, string>>& in_idx);
static void convert_node(PathPositionHandleGraph& graph, GAMInfo* gam_info, handle_t handle,
                         path_handle_t ref_path_handle, LW* output);
static void convert_chain(PathPositionHandleGraph& graph, SnarlDistanceIndex& distance_index, GAMInfo* gam_info,
                          net_handle_t chain, const string& ref_path);
static void convert_snarl(PathPositionHandleGraph& graph, SnarlDistanceIndex& distance_index, GAMInfo* gam_info,
                          net_handle_t snarl, path_handle_t ref_path_handle, LW* output);

void help(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <graph> " << endl
       << "Chop out path intervals from a vg graph" << endl
       << endl
       << "options: " << endl
       << "    -p, --progress            Print progress" << endl
       << "    -d, --dist FILE           Distance index from vg index -j [REQUIRED]" << endl
       << "    -r, --ref-prefix NAME     Prefix of reference path(s) [REQUIRED]" << endl
       << "    -g, --gam FILE            Sorted GAM file. Must have .gai index. Make both with \"vg gamsort x.gam -i x.sort.gam.gai > x.sort.gam\"" << endl
       << "    -t, --threads N           Number of threads to use [default: all available]" << endl      
       << endl;
}    

int main(int argc, char** argv) {

    string ref_path_prefix;
    string distance_index_filename;
    string gam_filename;
    bool progress = false;
    int c;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"progress", no_argument, 0, 'p'},
            {"ref-path", required_argument, 0, 'r'},
            {"dist", required_argument, 0, 'd'},
            {"gam", required_argument, 0, 'g'},
            {"threads", required_argument, 0, 't'},            
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hpr:d:g:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'p':
            progress = true;
            break;
        case 'r':
            ref_path_prefix = optarg;
            break;
        case 'd':
            distance_index_filename = optarg;
            break;
        case 'g':
            gam_filename = optarg;
            break;
        case 't':
        {
            int num_threads = stoi(optarg);
            if (num_threads <= 0) {
                cerr << "[vg2maf] error: Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                exit(1);
            }
            omp_set_num_threads(num_threads);
            break;
        }                                    
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 1) {
        help(argv);
        return 1;
    }
    if (distance_index_filename.empty()) {
        cerr << "[vg2maf] error: Distance index must be specified with -d" << endl;
        return 1;
    }
    if (ref_path_prefix.empty()) {
        cerr << "[vg2maf] error: Reference path must be speficied with -r" << endl;
        return 1;
    }

    string graph_filename = argv[optind++];
    ifstream graph_stream(graph_filename);
    if (!graph_stream) {
        cerr << "[vg2maf] error: Unable to open input graph " << graph_filename << endl;
        return 1;
    }    
    unique_ptr<PathHandleGraph> base_graph = load_graph(graph_stream);
    graph_stream.close();
    if (progress) {
        cerr << "[vg2maf]: Loaded graph" << endl;
    }
    bdsg::ReferencePathOverlayHelper overlay_helper;
    PathPositionHandleGraph* graph = overlay_helper.apply(base_graph.get());
    if (progress && dynamic_cast<PathPositionHandleGraph*>(base_graph.get()) == nullptr) {
        cerr << "[vg2maf]: Applied position overlay" << endl;
    }
    
    SnarlDistanceIndex distance_index;
    distance_index.deserialize(distance_index_filename);
    if (progress) {
        cerr << "[vg2maf]: Loaded distance index" << endl;
    }

    unique_ptr<GAMInfo> gam_info;
    ifstream gam_file;
    if (!gam_filename.empty()) {
        gam_info.reset(new GAMInfo());
        gam_file.open(gam_filename);
        if (!gam_file) {
            cerr << "[vg2maf] error: Unable to open gam file " << gam_filename << endl;
            return 1;
        }
        gam_info->cursor = vg::GAMIndex::cursor_t(gam_file);
        string gam_index_filename = gam_filename + ".gai";
        ifstream gam_index_file(gam_index_filename);
        if (!gam_index_file) {
            cerr << "[vg2maf] error: Unable to open gam index " << gam_index_filename << endl;
            return 1;
        }
        gam_info->index.load(gam_index_file);
        if (progress) {
            cerr << "[vg2maf]: Loaded GAM index" << endl;
        }
    }

    // iterate the top-level chains
    distance_index.for_each_child(distance_index.get_root(), [&](net_handle_t net_handle) {
        if (distance_index.is_chain(net_handle)) {
            convert_chain(*graph, distance_index, gam_info.get(), net_handle, ref_path_prefix);
        }        
    });

    return 0;
}

// copied from vg/src/path.cpp
void reverse_complement_mapping_in_place(vg::Mapping* m, const function<int64_t(nid_t)>& node_length) {
            
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
        // compute the number of columns (just sum of all lengths since we're not aligning)
        for (const auto& ie : in_map) {
            width += ie.second.length();
        }
        // allocate the output alignment rows
        for (const auto& ie : in_map) {
            out_map[ie.first].resize(width);
        }
        // write the unaligned output alignment
        int64_t cur_row = 0;
        int64_t cur_offset = 0;
        // do it column by column
        for (int64_t col = 0; col < width; ++col) {
            for (const auto& ie : in_map) {
                int64_t i_row = ie.first;
                // write the sequence for each row
                if (cur_row == i_row) {
                    out_map[i_row][col] = ie.second[cur_offset];
                    ++cur_offset;
                    if (cur_offset == ie.second.length()) {
                        cur_offset = 0;
                        ++cur_row;
                    }
                } else {
                    // and the rest is unaligned
                    out_map[i_row][col] = '-';
                }
            }            
        }
        assert(cur_offset == 0 && cur_row == in_elem.second.size());
    }

    return out_idx;
}

void convert_node(PathPositionHandleGraph& graph, GAMInfo* gam_info, handle_t handle, path_handle_t ref_path_handle, LW* output) {

    // we don't care about the chain orientation, everything below is based on the underlying node being forward
    if (graph.get_is_reverse(handle)) {
        handle = graph.flip(handle);
    }
    
    vector<step_handle_t> steps = graph.steps_of_handle(handle);
    if (steps.empty()) {
        //cerr << "[vg2maf] warning: Skipping node " << graph.get_id(handle) << " because there are no paths on it" << endl;
        return;
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

    // convert the gam indxes to rows
    if (gam_info) {
        // query the index for our node [todo: is this too simplisitic due to repeated queries?]
        vector<vg::Mapping> mappings;
        vector<string> names;
        vector<int64_t> start_positions;
        vector<int64_t> sequence_lengths;
        vector<bool> mapping_reversed;
#ifdef debug
        cerr << "doing gam index query on node " << node_id << endl;
#endif
        gam_info->index.find(gam_info->cursor, node_id, [&](const vg::Alignment& aln) {            
            // this is the position w.r.t the read alignment
            int64_t pos = 0;
            for (int64_t i = 0; i < aln.path().mapping_size(); ++i) {
                const vg::Mapping& mapping = aln.path().mapping(i);
                if (mapping.position().node_id() == node_id) {
                    mappings.push_back(mapping);
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
        });

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

        // copy our mappings into the alignment rows
        for (int64_t i = 0; i < mappings.size(); ++i) {
            vg::Mapping& mapping = mappings[i];
            Alignment_Row* row = (Alignment_Row*)st_calloc(1, sizeof(Alignment_Row));
            row->sequence_name = stString_copy(names[i].c_str());
            row->length = 0;
            row->sequence_length = sequence_lengths[i];
            row->start = start_positions[i];
            row->strand = mapping_reversed[i] ? 0 : 1;
            row->bases = (char*)st_calloc(alignment->column_number, sizeof(char));
            // add the opening gaps
            int64_t col = 0;
            int64_t node_offset = 0;
            for (; col < mappings[i].position().offset(); ++col) {
                row->bases[col] = '-';
            }
            node_offset = col;
            // add in the sequence
            for (int64_t j = 0; j < mappings[i].edit_size(); ++j) {
                const vg::Edit& edit = mappings[i].edit(j);
                if (edit.from_length() == edit.to_length()) {
                    //match
                    for (int64_t k = 0; k < edit.from_length(); ++k) {
                        if (!edit.sequence().empty()) {
                            row->bases[col] = edit.sequence()[k];                            
                        } else {
                            row->bases[col] = node_sequence[node_offset];
                        }
                        ++col;
                        ++node_offset;
                        ++row->length;
                    }
                } else if (edit.to_length() == 0 && edit.from_length() > 0) {
                    // delete
                    for (int64_t k = 0; k < edit.from_length(); ++k) {
                        row->bases[col++] = '-';
                        ++node_offset;
                    }
                    
                } else {
                    // insert
                    assert(edit.from_length() < edit.to_length());
                    // add the common part [todo: this should probably not be aligned automaticall going forward]
                    for (int64_t k = 0; k < edit.from_length(); ++k) {
                        row->bases[col++] = edit.sequence()[k];
                        ++node_offset;
                        ++row->length;
                    }
                    // look up the insertion from the index
                    const unordered_map<int64_t, string>& row_alignments = ins_alignments.at(col);
                    // find the row
                    const string& row_string = row_alignments.at(i);
                    for (int64_t k = 0; k < row_string.length(); ++k) {
                        row->bases[col++] = row_string[k];
                        if (row_string[k] != '-') {
                            ++row->length;
                        }
                    }
                }
            }

            // add the closing gaps
            for (int64_t j = node_offset; j < node_sequence.length(); ++j) {
                row->bases[j] = '-';
            }
            
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
        for (int64_t col = 0; col < node_sequence.length(); ++col) {
            // add insertion gaps if the column is in the insertion index
            if (ins_alignments.count(col)) {
                int64_t gaps = ins_alignments[col].begin()->second.length();
                for (int64_t k = 0; k < gaps; ++k) {
                    row->bases[maf_col++] = '-';
                }
                
            }
            row->bases[maf_col++] = node_sequence[col];
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

    // write the alignment    
    maf_write_block(alignment, output);

    alignment_destruct(alignment, true);
}

void convert_snarl(PathPositionHandleGraph& graph, SnarlDistanceIndex& distance_index, GAMInfo* gam_info,
                   net_handle_t snarl, path_handle_t ref_path_handle, LW* output) {

    net_handle_t start_bound = distance_index.get_bound(snarl, false, true);
    net_handle_t end_bound = distance_index.get_bound(snarl, true, false);

    handle_t start_handle = distance_index.get_handle(start_bound, &graph);
    handle_t end_handle = distance_index.get_handle(end_bound, &graph);

#ifdef debug
    cerr << "snarl goes from " << graph.get_id(start_handle) << ":" << graph.get_is_reverse(start_handle) << " to "
         << graph.get_id(end_handle) << ":" << graph.get_is_reverse(end_handle) << endl;
#endif

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
            convert_node(graph, gam_info, handle, ref_path_handle, output);
            visited.insert(handle);
            graph.follow_edges(handle, false, [&](handle_t other_handle) {
                bfs_queue.push_back(other_handle);
            });
            graph.follow_edges(handle, true, [&](handle_t other_handle) {
                bfs_queue.push_back(other_handle);
            });
        }
    }
}

void convert_chain(PathPositionHandleGraph& graph, SnarlDistanceIndex& distance_index, GAMInfo* gam_info,
                   net_handle_t chain, const string& ref_path) {

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
        cerr <<"[vg2maf] warning: No reference path found through chain from " 
             << graph.get_id(start_handle) << ":" << graph.get_is_reverse(start_handle) << " to "
             << graph.get_id(end_handle) << ":" << graph.get_is_reverse(end_handle) << endl;
    }

    if (end_ref_paths.size() > 1) {
        cerr << "[vg2maf] warning: " << end_ref_paths.size() << " reference paths found through chain from "
             << graph.get_id(start_handle) << ":" << graph.get_is_reverse(start_handle) << " to "
             << graph.get_id(end_handle) << ":" << graph.get_is_reverse(end_handle) << ". Just choosing first" << endl;
        cerr << "the paths are "; for (const auto xx : end_ref_paths) cerr << graph.get_path_name(xx.first) << ", "; cerr << endl;
    }

    path_handle_t ref_path_handle = end_ref_paths.begin()->first;
    //todo: this may fall down with cyclic ref
    bool ref_path_reversed = end_ref_paths.begin()->second < start_ref_paths.at(ref_path_handle);   
    if (ref_path_reversed) {
        // todo: why does this seem to have no effect on iteration below.
        chain = distance_index.flip(chain);
    }

    LW *output = LW_construct(stdout, false);

    Tag* tag = tag_construct((char*)"version", (char*)"1", NULL);
    maf_write_header(tag, output);

    vector<net_handle_t> chain_childs;
    distance_index.for_each_child(chain, [&](net_handle_t net_handle) {
        chain_childs.push_back(net_handle);
    });
    if (ref_path_reversed) {
        std::reverse(chain_childs.begin(), chain_childs.end());
    }

    // convert the chain, one node/snarl at a time
    for (int64_t i = 0; i < chain_childs.size(); ++i) {
        net_handle_t net_handle = chain_childs[i];
        if (ref_path_reversed) {
            net_handle = distance_index.flip(net_handle);
        }
        if (distance_index.is_node(net_handle)) {            
            convert_node(graph, gam_info, distance_index.get_handle(net_handle, &graph), ref_path_handle, output);
        } else if (distance_index.is_snarl(net_handle)) {
            convert_snarl(graph, distance_index, gam_info, net_handle, ref_path_handle, output);
        } else if (distance_index.is_chain(net_handle)) {
            cerr << "TODO: CHILD CHAIN" << endl;
        } else {
            assert(false);
        }
    }
    
    LW_destruct(output, false);
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
