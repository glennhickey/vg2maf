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

//#define debug

using namespace std;
using namespace handlegraph;
using namespace bdsg;


// from hal2vg/clip-vg.cpp
static unique_ptr<PathHandleGraph> load_graph(istream& graph_stream);


void help(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <graph> " << endl
       << "Chop out path intervals from a vg graph" << endl
       << endl
       << "options: " << endl
       << "    -p, --progress            Print progress" << endl
       << "    -d, --dist FILE           Distance index from vg index -j" << endl
       << "    -r, --ref-prefix NAME     Prefix of reference path(s) [REQUIRED]" << endl
       << "    -g, --gam FILE            Sorted GAM file. Must have .gai index. Make both with \"vg gamsort x.gam -i x.sort.gam.gai > x.sort.gam\"" << endl
       << "    -u, --unaligned-ins       Leave GAM insertions unaligned\n"
       << "    -o, --only-ref            Only convert nodes on reference path(s)" << endl
       << "    -t, --threads N           Number of threads to use [default: all available]" << endl      
       << endl;
}    

int main(int argc, char** argv) {

    string ref_path_prefix;
    string distance_index_filename;
    string gam_filename;
    bool align_insertions = true;
    bool taf_output = false;
    bool ref_only = false;
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
            {"unaligned-ins", no_argument, 0, 'u'},
            {"only-ref", no_argument, 0, 'o'},
            {"threads", required_argument, 0, 't'},            
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hpr:d:g:uot:",
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
        case 'u':
            align_insertions = false;
            break;
        case 'o':
            ref_only = true;
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
    if (progress) {
        cerr << "[vg2maf]: Using " << get_thread_count() << (get_thread_count() > 1 ? " threads" : " thread") << endl;
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
    
    if (distance_index_filename.empty()) {
        distance_index_filename = graph_filename.substr(0, graph_filename.find_last_of(".")) + ".dist";
        cerr << "[vg2maf]: Assuming distance index is " << distance_index_filename << endl;
    }
    SnarlDistanceIndex distance_index;
    distance_index.deserialize(distance_index_filename);
    if (progress) {
        cerr << "[vg2maf]: Loaded distance index" << endl;
    }

    vector<unique_ptr<GAMInfo>> gam_info_vec;
    gam_info_vec.resize(get_thread_count());
    vector<GAMInfo*> gam_info_ptrs(get_thread_count(), nullptr);
    vector<ifstream> gam_file_vec(get_thread_count());
    if (!gam_filename.empty()) {
        for (int i = 0; i < gam_info_vec.size(); ++i) {
            unique_ptr<GAMInfo>& gam_info = gam_info_vec[i];
            ifstream& gam_file = gam_file_vec[i];
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
            if (progress && i == 0) {
                cerr << "[vg2maf]: Loaded GAM index" << endl;
            }
            gam_info_ptrs[i] = gam_info.get();
        }
    }

    LW *output = LW_construct(stdout, false);

    Tag* tag = tag_construct((char*)"version", (char*)"1", NULL);
    if (taf_output) {
        taf_write_header(tag, output);
    } else {
        maf_write_header(tag, output);
    }
    tag_destruct(tag);

    abpoa_para_t* abpoa_params = nullptr;
    if (align_insertions) {
        abpoa_params = construct_abpoa_params();
    }

    // iterate the top-level chains
    int64_t i = 0;
    int64_t num_chains = 0;
    if (progress) {
        distance_index.for_each_child(distance_index.get_root(), [&](net_handle_t net_handle) {
            if (distance_index.is_chain(net_handle)) {
                ++num_chains;
            }
        });
    }
    distance_index.for_each_child(distance_index.get_root(), [&](net_handle_t net_handle) {
        if (distance_index.is_chain(net_handle)) {
            convert_chain(*graph, distance_index, gam_info_ptrs, net_handle, ref_path_prefix, progress,
                          make_pair(++i, num_chains), abpoa_params, ref_only, output);
        }        
    });

    LW_destruct(output, false);

    if (progress) {
        cerr << "[vg2maf]: Finished" << endl;
    }

    if (abpoa_params) {
        abpoa_free_para(abpoa_params);
    }
    return 0;
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
