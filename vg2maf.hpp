#pragma once

#include "handlegraph/path_handle_graph.hpp"
#include "handlegraph/path_position_handle_graph.hpp"
#include "stream_index.hpp"
#include "bdsg/snarl_distance_index.hpp"

extern "C" {
#include "taf.h"
#include "sonLib.h"
}
#include "abpoa.h"

using namespace std;
using namespace handlegraph;
using namespace bdsg;

// start with really barebones placeholder implementation
struct GAMInfo {
    vg::GAMIndex index;
    vg::GAMIndex::cursor_t cursor;
};


// extract all insertions from the edits in a list of mappings
// the index returned is of form [col][row] -> inserted string
// col is the position in the node
// row is the position in the mappings array
unordered_map<int64_t, unordered_map<int64_t, string>> get_insertion_index(const vector<vg::Mapping>& mappings);

// turn the insertion index (above) into an alignment.
// for now, it just adds gaps and leaves all insertions unaligned
// todo: run through multiple aligner like abpoa
unordered_map<int64_t, unordered_map<int64_t, string>> align_insertion_index(const unordered_map<int64_t, unordered_map<int64_t, string>>& in_idx);

// actually align the insertion index (new default). above function kept around for tests for now. 
unordered_map<int64_t, unordered_map<int64_t, string>> abpoa_align_insertion_index(
    abpoa_para_t* abpoa_params, const unordered_map<int64_t, unordered_map<int64_t, string>>& in_idx);

// scan through a given path, returning the list of nodes that traverses the two handles
// todo: there's a lot of code in vg for this, perhaps reuse?
vector<handle_t> get_ref_traversal(PathPositionHandleGraph& graph, path_handle_t ref_path_handle,
                                   handle_t start_handle, handle_t end_handle);

// convert a node to maf
// alignment object must be freed with alignmenet_destruct(alignment, true)
Alignment* convert_node(PathPositionHandleGraph& graph, const vector<vg::Alignment>& gam_alignments, handle_t handle,
                        path_handle_t ref_path_handle, abpoa_para_t* abpoa_params);

// converts a batch of nodes to maf
// this is used in attempt to aggregate gam index queries into ranges.
// ranges_start and range_end give an open-ended interval in sorted_nodes
// sorted_nodes stores offsets in node_buffer and out_alignment_buffer
void convert_node_range(PathPositionHandleGraph& graph, GAMInfo* gam_info, const vector<handle_t>& node_buffer,
                        const vector<int64_t>& sorted_nodes, int64_t range_start, int64_t range_end,
                        path_handle_t ref_path_handle, abpoa_para_t* abpoa_params, vector<Alignment*>& out_alignment_buffer);


// convert a chain to maf, by scanning its children in order
void convert_chain(PathPositionHandleGraph& graph, SnarlDistanceIndex& distance_index, vector<GAMInfo*>& gam_info,
                   net_handle_t chain, const string& ref_path, bool progress, const pair<int64_t, int64_t>& chain_idx,
                   abpoa_para_t* abpoa_params, bool taf_output, LW* output);

// return the handles inside a snarl in the order that we want them in the maf
// todo: this function probably needs some work to effectively put complex regions through taffy norm
void traverse_snarl(PathPositionHandleGraph& graph, SnarlDistanceIndex& distance_index,
                    net_handle_t snarl, path_handle_t ref_path_handle,
                    bool ref_path_reversed, vector<handle_t>& out_handles);

// make an abpoa parameters object, using defautls from cactus
abpoa_para_t* construct_abpoa_params();

// use ascii phred inside TAF (as it's a text-based format)
// from https://en.wikipedia.org/wiki/FASTQ_format
// The byte representing quality runs from 0x21 (lowest quality; '!' in ASCII) to 0x7e (highest quality; '~' in ASCII). 
inline unsigned char phred_byte_to_ascii(unsigned char b) {
    return b >= 0x7e - 0x21 ? 0x7e : b + 0x21;
}

