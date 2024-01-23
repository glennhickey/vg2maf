#pragma once

#include "handlegraph/path_handle_graph.hpp"
#include "handlegraph/path_position_handle_graph.hpp"
#include "stream_index.hpp"
#include "bdsg/snarl_distance_index.hpp"

extern "C" {
#include "taf.h"
#include "sonLib.h"
}

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

// scan through a given path, returning the list of nodes that traverses the two handles
// todo: there's a lot of code in vg for this, perhaps reuse?
vector<handle_t> get_ref_traversal(PathPositionHandleGraph& graph, path_handle_t ref_path_handle,
                                   handle_t start_handle, handle_t end_handle);

// convert a node to maf
// if gam_info is given, extract all relevant mappings and add them too
// alignment object must be freed with alignmenet_destruct(alignment, true)
Alignment* convert_node(PathPositionHandleGraph& graph, GAMInfo* gam_info, handle_t handle,
                        path_handle_t ref_path_handle);

// convert a chain to maf, by scanning its children in order
void convert_chain(PathPositionHandleGraph& graph, SnarlDistanceIndex& distance_index, vector<GAMInfo*>& gam_info,
                   net_handle_t chain, const string& ref_path);

// return the handles inside a snarl in the order that we want them in the maf
// todo: this function probably needs some work to effectively put complex regions through taffy norm
void traverse_snarl(PathPositionHandleGraph& graph, SnarlDistanceIndex& distance_index,
                    net_handle_t snarl, path_handle_t ref_path_handle, vector<handle_t>& out_handles);
