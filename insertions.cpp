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

#include "abpoa.h"

using namespace std;

//#define debug

// char <--> uint8_t conversion copied over from abPOA example
// AaCcGgTtNn ==> 0,1,2,3,4
static unsigned char nst_nt4_table[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

// 65,97=>A, 67,99=>C, 71,103=>G, 84,85,116,117=>T, else=>N
static const char nst_nt256_table[256] = {
       'A', 'C', 'G', 'T',  'N', '-', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', '-',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N'
};

static inline char msa_to_base(uint8_t n) {
    return (char)nst_nt256_table[n];
}

static inline uint8_t msa_to_byte(char c) {
    return nst_nt4_table[(int)c];
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

unordered_map<int64_t, unordered_map<int64_t, string>> abpoa_align_insertion_index(
    abpoa_para_t* abpoa_params, const unordered_map<int64_t, unordered_map<int64_t, string>>& in_idx) {

    unordered_map<int64_t, unordered_map<int64_t, string>> out_idx;

    // init abpoa (todo: move this up (maybe even to main?)
    abpoa_t *ab = abpoa_init();
    
    for (const auto& ie : in_idx) {
        if (ie.second.size() == 1) {
            // nothing to do with single row, just copy it over
            out_idx[ie.first] = ie.second;
        } else {
            unordered_map<int64_t, string>& oe = out_idx[ie.first];
            int seq_no = ie.second.size();
            int* seq_lens = (int*)st_calloc(ie.second.size(), sizeof(int));
            // allocate the poa input buffer
            uint8_t **bseqs = (uint8_t**)st_malloc(sizeof(uint8_t*) * seq_no);
            int64_t row_idx = 0;
            vector<int64_t> row_nums;
            for (const auto& ri : ie.second) {
                seq_lens[row_idx] = (int)ri.second.length();
                bseqs[row_idx] = (uint8_t*)st_malloc(sizeof(uint8_t) * ri.second.length());
                // copy in the row bases, doing the abpoa conversion
                for (int64_t col_idx = 0; col_idx < ri.second.length(); ++col_idx) {
                    bseqs[row_idx][col_idx] = msa_to_byte(ri.second[col_idx]);
                }
                ++row_idx;
                row_nums.push_back(ri.first);
            }
                        
            abpoa_msa(ab, abpoa_params, seq_no, NULL, seq_lens, bseqs, NULL, NULL);

            // copy the alignment back into the map
            row_idx = 0;
            for (const auto& ri : ie.second) {
                string& out_row = oe[row_nums[row_idx]];
                out_row.resize(ab->abc->msa_len);
                for (int64_t col = 0; col < ab->abc->msa_len; ++col) {
                    out_row[col] = msa_to_base(ab->abc->msa_base[row_idx][col]);
                }
                ++row_idx;
            }

            free(seq_lens);
            free(bseqs);

        }
    }

#ifdef debug
    cerr << " abpoa input" << endl;
    for (auto xx : in_idx) {
        for (auto yy : xx.second) {
            cerr  << "INPUT[" << xx.first << "][" << yy.first <<"]=" << yy.second << endl;
        }
    }

    
    cerr << " abpoa output" << endl;
    for (auto xx : out_idx) {
        for (auto yy : xx.second) {
            cerr  << "INS[" << xx.first << "][" << yy.first <<"]=" << yy.second << endl;
        }
    }
#endif

    // todo: push up with init
    abpoa_free(ab);

    return out_idx;
}

// from cactus/bar/impl/poaBarAligner.c
abpoa_para_t* construct_abpoa_params() {
    abpoa_para_t *abpt = abpoa_init_para();

    // output options
    abpt->out_msa = 1; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
    abpt->out_cons = 0; // generate consensus sequence, set 0 to disable

    // alignment mode. 0:global alignment, 1:local, 2:extension
    // only global works
    abpt->align_mode = ABPOA_GLOBAL_MODE;

    // banding parameters
    abpt->wb = 300;
    abpt->wf = 0.05;

    // gap scoring model
    abpt->gap_open1 = 400;
    abpt->gap_ext1 = 30;
    abpt->gap_open2 = 1200;
    abpt->gap_ext2 = 1;
    
    // seeding paramters
    abpt->disable_seeding = 1;
    abpt->k = 19;
    abpt->w = 10;
    abpt->min_w = 500;

    // progressive toggle
    abpt->progressive_poa = 1;

    // generate the substitution matrix
    abpt->use_score_matrix = 0;
    abpoa_post_set_para(abpt);

    // optionally override the substitution matrix
    char *submat_string = stString_copy("91 -114 -61 -123 -100 -114 100 -125 -61 -100 -61 -125 100 -114 -100 -123 -61 -114 91 -100 -100 -100 -100 -100 100");
    if (submat_string && strlen(submat_string) > 0) {
        // Note, this will be used to explicitly override abpoa's subsitution matrix just before aligning
        abpt->use_score_matrix = 1;
        assert(abpt->m == 5);
        int count = 0;
        for (char* val = strtok(submat_string, " "); val != NULL; val = strtok(NULL, " ")) {
            abpt->mat[count++] = atoi(val);
        }
        assert(count == 25);
        int i; abpt->min_mis = 0, abpt->max_mat = 0;
        for (i = 0; i < abpt->m * abpt->m; ++i) {
            if (abpt->mat[i] > abpt->max_mat)
                abpt->max_mat = abpt->mat[i];
            if (-abpt->mat[i] > abpt->min_mis) 
                abpt->min_mis = -abpt->mat[i];
        }
    }
    free(submat_string);
    return abpt;
}
