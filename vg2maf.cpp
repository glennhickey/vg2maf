#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <fstream>
#include <memory>
#include <unistd.h>
#include <getopt.h>


using namespace std;

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

using namespace std;
using namespace handlegraph;
using namespace bdsg;

// from hal2vg/clip-vg.cpp
static unique_ptr<PathHandleGraph> load_graph(istream& graph_stream);

// start with really barebones placeholde implementation
static void convert_node(PathPositionHandleGraph& graph, vg::GAMIndex gam_index, handle_t handle,
                         path_handle_t ref_path_handle, LW* output);
static void convert_chain(PathPositionHandleGraph& graph, SnarlDistanceIndex& distance_index, vg::GAMIndex* gam_index,
                          net_handle_t chain, const string& ref_path);
static void convert_snarl(PathPositionHandleGraph& graph, SnarlDistanceIndex& distance_index, vg::GAMIndex* gam_index,
                          net_handle_t snarl, path_handle_t ref_path_handle, LW* output);

void help(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <graph> <distance-index>" << endl
       << "Chop out path intervals from a vg graph" << endl
       << endl
       << "options: " << endl
       << "    -p, --progress            Print progress" << endl
       << "    -d, --dist FILE           Distance index from vg index -j [REQUIRED]" << endl
       << "    -r, --ref-prefix NAME     Name or prefix of reference path [REQUIRED]" << endl
       << "    -g, --gam FILE            GAM file. Must have .gai index from \"vg gamsort x.gam -i x.sort.gai > x.sort.gam\"" << endl
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

    unique_ptr<vg::GAMIndex> gam_index;
    if (!gam_filename.empty()) {
        string gam_index_filename = gam_filename + ".gai";
        ifstream gam_index_file(gam_index_filename);
        if (!gam_index_file) {
            cerr << "vg2maf]: Unable to open gam index " << gam_index_filename << endl;
            return 1;
        }
        gam_index.reset(new vg::GAMIndex());
        gam_index->load(gam_index_file);
        if (progress) {
            cerr << "[vg2maf]: Loaded GAM index" << endl;
        }
    }

    // iterate the top-level chains
    distance_index.for_each_child(distance_index.get_root(), [&](net_handle_t net_handle) {
        if (distance_index.is_chain(net_handle)) {
            convert_chain(*graph, distance_index, gam_index.get(), net_handle, ref_path_prefix);
        }        
    });

    return 0;
}

void convert_node(PathPositionHandleGraph& graph, vg::GAMIndex* gam_index, handle_t handle, path_handle_t ref_path_handle, LW* output) {

    vector<step_handle_t> steps = graph.steps_of_handle(handle);
    if (steps.empty()) {
        cerr << "[vg2maf] warning: Skipping node " << graph.get_id(handle) << " because there are no paths on it" << endl;
        return;
    }
    
    Alignment* alignment = (Alignment*)st_calloc(1, sizeof(Alignment));
    alignment->row_number = steps.size();
    alignment->column_number = graph.get_length(handle);
    alignment->column_tags = (Tag**)st_calloc(alignment->column_number, sizeof(Tag*));

    // convert each step to row
    vector<Alignment_Row*> rows;

    for (step_handle_t& step_handle : steps) {
        Alignment_Row* row = (Alignment_Row*)st_calloc(1, sizeof(Alignment_Row));
        path_handle_t step_path_handle = graph.get_path_handle_of_step(step_handle);
        row->sequence_name = stString_copy(graph.get_path_name(step_path_handle).c_str());
        row->start = graph.get_position_of_step(step_handle);
        row->length = alignment->column_number;
        row->sequence_length = graph.get_path_length(step_path_handle);
        // todo: check this strand logic
        row->strand = graph.get_is_reverse(handle) == graph.get_is_reverse(graph.get_handle_of_step(step_handle));
        if (row->strand == 0) {
            row->start += row->length - 1;
        }
        row->bases = stString_copy(graph.get_sequence(handle).c_str());
        rows.push_back(row);
    }

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
        } else if (row != alignment->row) {
            cur_row->n_row = row;
        }
        cur_row = row;
    }

    // write the alignment    
    maf_write_block(alignment, output);

    alignment_destruct(alignment, true);
}

void convert_snarl(PathPositionHandleGraph& graph, SnarlDistanceIndex& distance_index, vg::GAMIndex* gam_index,
                   net_handle_t snarl, path_handle_t ref_path_handle, LW* output) {

    net_handle_t start_bound = distance_index.get_bound(snarl, false, true);
    net_handle_t end_bound = distance_index.get_bound(snarl, true, false);

    handle_t start_handle = distance_index.get_handle(start_bound, &graph);
    handle_t end_handle = distance_index.get_handle(end_bound, &graph);

    cerr << "snarl goes from " << graph.get_id(start_handle) << ":" << graph.get_is_reverse(start_handle) << " to "
         << graph.get_id(end_handle) << ":" << graph.get_is_reverse(end_handle) << endl;

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
            convert_node(graph, gam_index, handle, ref_path_handle, output);
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

void convert_chain(PathPositionHandleGraph& graph, SnarlDistanceIndex& distance_index, vg::GAMIndex* gam_index,
                   net_handle_t chain, const string& ref_path) {

    net_handle_t start_bound = distance_index.get_bound(chain, false, true);
    net_handle_t end_bound = distance_index.get_bound(chain, true, false);

    handle_t start_handle = distance_index.get_handle(start_bound, &graph);
    handle_t end_handle = distance_index.get_handle(end_bound, &graph);

    cerr << "chain goes from " << graph.get_id(start_handle) << ":" << graph.get_is_reverse(start_handle) << " to "
         << graph.get_id(end_handle) << ":" << graph.get_is_reverse(end_handle) << endl;

    set<pair<path_handle_t, bool>> start_ref_paths;
    set<pair<path_handle_t, bool>> end_ref_paths;
    graph.for_each_step_on_handle(start_handle, [&](step_handle_t step_handle) {
        path_handle_t step_path_handle = graph.get_path_handle_of_step(step_handle);
        string path_name = graph.get_path_name(step_path_handle);
        if (path_name.compare(0, ref_path.length(), ref_path) == 0) {
            bool reversed = graph.get_is_reverse(graph.get_handle_of_step(step_handle)) == graph.get_is_reverse(start_handle);
            start_ref_paths.insert(make_pair(step_path_handle, reversed));
        }
    });
    if (!start_ref_paths.empty()) {
        graph.for_each_step_on_handle(end_handle, [&](step_handle_t step_handle) {
            path_handle_t step_path_handle = graph.get_path_handle_of_step(step_handle);
            string path_name = graph.get_path_name(step_path_handle);
            bool reversed = graph.get_is_reverse(graph.get_handle_of_step(step_handle)) == graph.get_is_reverse(end_handle);
            auto val = make_pair(step_path_handle, reversed);
            if (path_name.compare(0, ref_path.length(), ref_path) == 0 && start_ref_paths.count(val)) {
                end_ref_paths.insert(val);
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
    }

    path_handle_t ref_path_handle = end_ref_paths.begin()->first;
    bool ref_path_reversed = end_ref_paths.begin()->second;
    if (ref_path_reversed) {
        // line the chain up to the path -- todo: may not be practical once we smarten up logic
        // but mc graphs always have forward reference paths, which hopefully persists to chain level
        handle_t temp_handle = start_handle;
        start_handle = graph.flip(end_handle);
        end_handle = graph.flip(temp_handle);
    }

    LW *output = LW_construct(stdout, false);

    Tag* tag = tag_construct((char*)"version", (char*)"1", NULL);
    maf_write_header(tag, output);

    // convert the chain, one node/snarl at a time
    distance_index.for_each_child(chain, [&](net_handle_t net_handle) {
        if (distance_index.is_node(net_handle)) {            
            convert_node(graph, gam_index, distance_index.get_handle(net_handle, &graph), ref_path_handle, output);
        } else if (distance_index.is_snarl(net_handle)) {
            convert_snarl(graph, distance_index, gam_index, net_handle, ref_path_handle, output);
        } else if (distance_index.is_chain(net_handle)) {
            cerr << "TODO: CHILD CHAIN" << endl;
        } else {
            assert(false);
        }
    });
    
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
