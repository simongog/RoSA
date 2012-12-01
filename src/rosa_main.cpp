#include <sdsl/int_vector.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/select_support.hpp>
#include <sdsl/rrr_vector.hpp>
#include <sdsl/gap_vector.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/wt_huff.hpp>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/wt.hpp>
#include <sdsl/util.hpp>
#include <sdsl/testutils.hpp>

/* typedefs for different index configuration
 *
 */
#ifdef FM_HUFF
#define INDEX_SUF "fm_huff"
#include "in_memory_index.hpp"
typedef in_memory_index<csa_wt<wt_huff<>,64000,64000> > tIDX;
#endif
#ifdef FM_HUFF_RRR
#define INDEX_SUF "fm_huff_rrr"
#include "in_memory_index.hpp"
#include <sdsl/rrr_vector.hpp>
typedef rrr_vector<63> tRRR;
typedef in_memory_index<csa_wt<wt_huff<tRRR, tRRR::rank_1_type, tRRR::select_1_type, tRRR::select_0_type>,64000,64000> > tIDX;
#endif
#ifdef FM_HUFF_RRR127
#define INDEX_SUF "fm_huff_rrr127"
#include "in_memory_index.hpp"
#include <sdsl/rrr_vector.hpp>
typedef rrr_vector<127> tRRR;
typedef in_memory_index<csa_wt<wt_huff<tRRR, tRRR::rank_1_type, tRRR::select_1_type, tRRR::select_0_type>,64000,64000> > tIDX;
#endif
#ifdef FM_RL
#define INDEX_SUF "fm_rlmn"
#include "in_memory_index.hpp"
typedef in_memory_index<csa_wt<wt_rlmn<sd_vector<>,
        sd_vector<>::rank_1_type, sd_vector<>::select_1_type>,64000,64000> >  tIDX;
#endif
#ifdef ROSA_BV
#define INDEX_SUF "rosa_bv"
#include "rosa.hpp"
typedef rosa<> tIDX;
#endif
#ifdef ROSA_RRR
#define INDEX_SUF "rosa_rrr"
#include "rosa.hpp"
typedef rosa<rrr_vector<> > tIDX;
#endif
#ifdef ROSA_RRR63
#define INDEX_SUF "rosa_rrr63"
#include "rosa.hpp"
typedef rosa<rrr_vector<63> > tIDX;
#endif
#ifdef ROSA_GAP
#define INDEX_SUF "rosa_gap"
#include "rosa.hpp"
typedef rosa<gap_vector<> > tIDX;
#endif
#ifndef INDEX_SUF
#define INDEX_SUF "rosa_sd"
#include "rosa.hpp"
typedef rosa<sd_vector<> > tIDX;
#endif

#include "rosa.hpp"
#include "pattern_file.hpp"



#include <iostream>// for input output
#include <cstdlib> // for atoll
#include <string>  // for filenames

#include <getopt.h>

using namespace sdsl;
using namespace std;



typedef bit_vector::size_type size_type;


template<class tIndex>
void generate_index(tIndex& idx, const char* file_name, size_type b, bool output_tikz=false, bool delete_tmp=false,
                    const char* tmp_file_dir="./", const char* output_dir=".")
{
    typedef bit_vector::size_type size_type;
    cerr << "call generate_index " << endl;

    string int_idx_file_name = tIndex::get_int_idx_filename(file_name,b, output_dir);
    string ext_idx_file_name = tIndex::get_ext_idx_filename(file_name,b, output_dir);

    ifstream int_idx_stream(int_idx_file_name.c_str());
    ifstream ext_idx_stream(ext_idx_file_name.c_str());
    if (!int_idx_stream or !ext_idx_stream) {
        if (!int_idx_stream) {
            cerr << "internal index is not stored in " << int_idx_file_name << endl;
        }
        if (!ext_idx_stream) {
            cerr << "external index is not stored in " << ext_idx_file_name << endl;
        }
        tIndex index(file_name, b, output_tikz, delete_tmp, tmp_file_dir, output_dir);
        util::store_to_file(index, int_idx_file_name.c_str());
    } else {
        cerr << "internal and external indexes exist." << endl;
        cerr << "Locations:" << endl;
        cerr << int_idx_file_name << endl;
        cerr << ext_idx_file_name << endl;
        close_stream_if_open(int_idx_stream);
        close_stream_if_open(ext_idx_stream);
    }
}


template<class tIndex>
void delete_index_files(const char* file_name, const char* output_dir, size_type b)
{
    string int_idx_file_name = tIndex::get_int_idx_filename(file_name,b, output_dir);
    string ext_idx_file_name = tIndex::get_ext_idx_filename(file_name,b, output_dir);

    std::remove(int_idx_file_name.c_str());
    std::remove(ext_idx_file_name.c_str());
    tIndex::remove_tmp_files(file_name, output_dir);
}


template<class tIndex>
void benchmark(const char* file_name, size_type b, const char* pattern_file_name, const char* tmp_file_dir, const char* output_dir,
               bool repeated_in_memory_search=false,
               bool only_in_memory=false, bool only_external_memory=false
              )
{
    std::cout << "# index_suf = " << INDEX_SUF << std::endl;
    std::cout<< "# repeated_in_memory_search = " << repeated_in_memory_search << std::endl;
    std::cout<< "# pattern_file_name = " << pattern_file_name << std::endl;
    stop_watch sw;
    {
        tIndex index;
        generate_index(index, file_name, b, false, false, tmp_file_dir, output_dir);
    }

    tIndex index;
    string int_idx_file_name = tIndex::get_int_idx_filename(file_name,b, output_dir);
    std::cerr << "load index from file " << int_idx_file_name << std::endl;
    sw.start();
    util::load_from_file(index, int_idx_file_name.c_str());
    index.set_file_name(string(file_name));
    index.set_output_dir(string(output_dir));
    sw.stop();
    std::cout << "# file_name = " << index.file_name << std::endl;
    std::cout << "# b = " << b << std::endl;
    cout << "# index_load_rtime = " << sw.get_real_time() << endl;
    cout << "# index_load_utime = " << sw.get_user_time() << endl;
    if (!only_external_memory) {
        pattern_file pf(pattern_file_name);
        pf.reset();
        cout << "# pattern_len = " << pf.pattern_len << endl;
        cout << "# pattern_swaps = " << pf.swaps << endl;
        size_type cnt=0, cnt1=0, cnt2=0, cnt3=0, rb, lb, d;
        cout << "# in_memory_queries = " << pf.pattern_cnt << endl;
        sw.start();
        for (size_type i=0; i < pf.pattern_cnt; ++i) {
            const unsigned char* c = (const unsigned char*)pf.get_next_pattern();
            cnt += index.get_interval(c, pf.pattern_len, lb, rb, d);
#ifdef ROSA_MAIN_DEBUG
            cout << i <<" " << d << "-[" << lb << "," << rb << "]" << endl;
            if (i == 0) {
                cout << "pattern=" << c <<endl;
            }
#endif
            cnt1 += lb; cnt2 += rb; cnt3 += d;
        }
        sw.stop();
        cout << "# cnt = " << cnt << endl;
        cout << "# cnt1 = " << cnt1 << endl;
        cout << "# cnt2 = " << cnt2 << endl;
        cout << "# cnt3 = " << cnt3 << endl;
        cout << "# avg_depth = " << ((double)cnt3)/pf.pattern_cnt << endl;
        cout << "# rtime = " << sw.get_real_time() << endl;
        cout << "# utime = " << sw.get_user_time() << endl;
        cout << "# stime = " << sw.get_sys_time() << endl;
    }

    if (!only_in_memory) {
        pattern_file pf(pattern_file_name);
        pf.reset();
        cout << "# full_queries = " << pf.pattern_cnt << endl;
        index.reset_counters();
        size_type cnt4=0;
        size_type cnt5=0;
        sw.start();
        for (size_type i=0, count=0; i < pf.pattern_cnt; ++i) {
//std::cerr << "i= " << i << std::endl;
//std::cerr << "pattern_len = " << pf.pattern_len << std::endl;
            count = index.count((const unsigned char*)pf.get_next_pattern(), pf.pattern_len, repeated_in_memory_search);
            cnt4 += count;
            cnt5 += (count > 0);
        }
        sw.stop();
        cout << "# cnt4 = " << cnt4 << endl;
        cout << "# cnt5 = " << cnt5 << endl;
        cout << "# rtime_full = " << sw.get_real_time() << endl;
        cout << "# utime_full = " << sw.get_user_time() << endl;
        cout << "# stime_full = " << sw.get_sys_time() << endl;
#ifdef OUTPUT_STATS
        double disk_access_per_query = 0.0;
        double gap_disk_access_per_query = 0.0;
        double avg_int_match_depth = 0.0;
        double ratio_of_int_matched_queries = 0.0;
        if (pf.pattern_cnt > 0) {
            disk_access_per_query = index.count_disk_access/(double)pf.pattern_cnt;
            gap_disk_access_per_query = index.count_gap_disk_access/(double)pf.pattern_cnt;
            avg_int_match_depth = index.count_int_steps/(double)pf.pattern_cnt;
        }
        cout << "# disk_access_per_query = " << disk_access_per_query << endl;
        cout << "# gap_disk_access_per_query = " << gap_disk_access_per_query << endl;
        cout << "# avg_int_match_depth = " << avg_int_match_depth << endl;
        if (index.count_queries) {
            ratio_of_int_matched_queries = index.count_int_match/(double)index.count_queries;
        }
        cout << "# ratio_of_int_matched_queries = " << ratio_of_int_matched_queries << endl;
        size_type ext_queries = index.count_queries-index.count_int_match;
        double block_length_per_ext_query = 0;
        if (ext_queries) {
            block_length_per_ext_query = index.count_block_length/(double)ext_queries;
        }
        cout << "# block_length_per_ext_query = " << block_length_per_ext_query << endl;
#endif
    }
}

void display_usage(char* command)
{
    cout << "usage of " << command << endl;
    cout << "  " << command << " --input_file=[input_file] "<< endl;
    cout << " input_file       : file name of the input text" << endl;
    cout << " threshold        : block size threshold b; default=4096" << endl;
    cout << " generate_index   : generates the index" << endl;
    cout << " generate_patterns: generate a pattern file; default=OFF" << endl;
    cout << " pattern_len      : length of each generated pattern; default=20" << endl;
    cout << " pattern_number   : number of generated patterns; default=100" << endl;
    cout << " pattern_swaps    : number of swaps applied to the patterns; default=0" << endl;
    cout << " pattern_min_occ  : minimum number of occurrences of a pattern in the text*; default=0" << endl;
    cout << " pattern_max_occ  : maximum number of occurrences of a pattern in the text*; default=0" << endl;
    cout << " benchmark        : run benchmark; default=0" << endl;
    cout << " benchmark_int    : run benchmark for in-memory part only; default=0" << endl;
    cout << " benchmark_ext    : run benchmark for external memory part only; default=0" << endl;
    cout << " pattern_file     : the pattern file for the benchmark" << endl;
    cout << " output_Hk        : output H_k for k=0..t of the text; default t=0" << endl;
    cout << " output_mem_info  : output information about the memory usage of the index" << endl;
    cout << " output_statistics: output statistics about the data structure" << endl;
    cout << " output_tikz      : output tikz code in the files: output_dir/basename(input_file).[f|b]wd_idx.tex" << endl;
    cout << " output_trie_nodes: output the number of nodes of a LOF-SA trie" << endl;
    cout << " delete_tmp_file  : delete the temporary files after the construction; default=0" << endl;
    cout << " tmp_file_dir     : directory for the temporary files; default=./" << endl;
    cout << " output_dir       : directory for the constructed indexes; default=dirname(input_file)" << endl;
    cout << " delete_index_file: delete the generated index files; default=0" << endl;
    cout << " verbose          : activate verbose mode." << endl;
    cout << " interactive      : enter interactive mode after creating/loading the index." << endl;
    cout << " greedy           : do a greedy parse." << endl;
	cout << " reconstruct_text : reconstruct text from factorization and condensed BWT." << endl;
	cout << " factor_occ_freq  : outputs for all x the number of factors which occur x times." << endl;
    cout << endl;
    cout << "* if pattern_min_occ=pattern_max_occ=0 this restriction is not used in the " << endl;
    cout << "  pattern generation process." << endl;
}

int main(int argc, char* argv[])
{
    std::string input_file_name = "";
    std::string pattern_file_name = "";
    std::string tmp_file_dir = "./";
    std::string output_dir = "./";

    size_type b = 4096;
    size_type plen = 20;
    size_type pno  = 100;
    size_type swapno  = 0;
    size_type pmin_occ = 0;
    size_type pmax_occ = 0;
    int repeated_in_memory_search = 0;
    int output_mem_info = 0;
    int output_tikz = 0;
    int output_statistics = 0;
    int output_trie_nodes = 0;
    int run_benchmark = 0;
    int run_benchmark_int = 0;
    int run_benchmark_ext = 0;
    int generate_patterns = 0;
    int do_generate_index = 0;
    int delete_tmp_file = 0;
    int delete_index_file = 0;
    size_type output_Hk = 0;
    size_type max_k = 0;
    int verbose = false;
    int interactive = 0;
    int greedy = 0;
	int reconstruct_text = 0;
	int factor_occ_freq = 0;

    int c;
    while (1) {
        static struct option long_options[] = {
            {"input_file", required_argument, 0, 'i'},
            {"threshold",  required_argument, 0, 'b'},
            {"tmp_file_dir", required_argument, 0, 't'},
            {"output_dir", required_argument, 0, 'o'},
            {"pattern_file", required_argument, 0, 'q'},
            {"pattern_len", required_argument, 0, 'm'},
            {"pattern_number", required_argument, 0, 'n'},
            {"pattern_swaps", required_argument, 0, 's'},
            {"pattern_min_occ", required_argument, 0, 'x'},
            {"pattern_max_occ", required_argument, 0, 'y'},
            {"generate_patterns", no_argument, &generate_patterns, 1},
            {"generate_index", no_argument, &do_generate_index, 1},
            {"repeated_in_memory_search", no_argument, &repeated_in_memory_search, 1},
            {"output_mem_info", no_argument, &output_mem_info, 1},
            {"output_tikz", no_argument, &output_tikz, 1},
            {"output_statistics", no_argument, &output_statistics, 1},
            {"output_trie_nodes", no_argument, &output_trie_nodes, 1},
            {"output_Hk", required_argument, 0, 'h'},
            {"benchmark", no_argument, &run_benchmark, 1},
            {"benchmark_int", no_argument, &run_benchmark_int, 1},
            {"benchmark_ext", no_argument, &run_benchmark_ext, 1},
            {"delete_tmp_file", no_argument, &delete_tmp_file, 1},
            {"delete_index_file", no_argument, &delete_index_file, 1},
            {"verbose", no_argument, &verbose, 1},
            {"interactive",no_argument, &interactive, 1},
            {"greedy",no_argument, &greedy, 1},
			{"reconstruct_text", no_argument, &reconstruct_text, 1},
			{"factor_occ_freq", no_argument, &factor_occ_freq, 1},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "i:b:", long_options, &option_index);
        if (c == -1) {  // detect the end of the options
            break;
        }
        string arg;
        switch (c) {
                // If this option set a flag, do nothing else now
                if (long_options[option_index].flag != 0)
                    break;
            case 'i':
                input_file_name = string(optarg);
                break;
            case 'q':
                pattern_file_name = string(optarg);
                break;
            case 't':
                tmp_file_dir = string(optarg);
                break;
            case 'o':
                output_dir = string(optarg);
                break;
            case 'b':
                b = atoll(optarg);
                break;
            case 'm':
                plen = atoll(optarg);
                break;
            case 'n':
                pno = atoll(optarg);
                break;
            case 's':
                swapno = atoll(optarg);
                break;
            case 'h':
                output_Hk = 1;
                max_k = atoll(optarg);
                break;
            case 'x':
                pmin_occ = atoll(optarg);
                break;
            case 'y':
                pmax_occ = atoll(optarg);
                break;
        }
    }

    if (input_file_name.size() == 0) {
        display_usage(argv[0]);
        return 1;
    }

    util::verbose = verbose;

    size_type input_size = (size_type)util::get_file_size(input_file_name.c_str());
    if (b > input_size) {
        std::cerr << "ERROR: input size ("<< input_size <<") is smaller than the block threshold b=" << b << endl;
        return 1;
    }


    if (do_generate_index) {
        tIDX index;
        generate_index(index, input_file_name.c_str(), b, false, delete_tmp_file, tmp_file_dir.c_str(), output_dir.c_str());
        string int_idx_file_name = tIDX::get_int_idx_filename(input_file_name.c_str(),b, output_dir.c_str());
        std::cerr << "load index from file " << int_idx_file_name << std::endl;
        util::load_from_file(index, int_idx_file_name.c_str());
        return 0;
    }

    if (output_tikz) {
        tIDX index;
        generate_index(index, input_file_name.c_str(), b, true, true, tmp_file_dir.c_str(), output_dir.c_str());
        return 0;
    }

    if (interactive) {
        tIDX index;
        string int_idx_file_name = tIDX::get_int_idx_filename(input_file_name.c_str(), b, output_dir.c_str());
        if (util::load_from_file(index, int_idx_file_name.c_str())) {
            std::cout<<"Entering interactive mode. Please enter patterns or Ctrl-D to exit."<<std::endl;
            std::string s;
            int mode=0;
            std::cout<<"Enter \"l\" for location mode and press any other key for counting mode."<<std::endl;
            std::cout<<"Followed by enter"<<std::endl;
            std::getline(cin, s);
            if ("l" == s) {
                mode = 1;
            }
            while (std::getline(cin, s)) {
                if (0 == mode) {
                    size_type cnt = index.count((const unsigned char*)s.c_str(), s.size(), false);
                    std::cout<<">>>>> pattern occurs "<<cnt<<" times"<<std::endl;
                } else {
                    std::vector<size_type> res;
                    size_type cnt = index.locate((const unsigned char*)s.c_str(), s.size(), res);
                    std::cout<<">>>>> pattern occurs "<<cnt<<" times"<<std::endl;
                    std::cout<<"locations:";
                    for (size_t i=0; i<res.size(); ++i) cout << " " << res[i];
                    cout << endl;
                }
            }
            return 0;
        } else
            return 1;
    }

    if (delete_index_file) {
        delete_index_files<tIDX>(input_file_name.c_str(), output_dir.c_str(), b);
        return 0;
    }

    if (output_mem_info or output_statistics or output_Hk or output_trie_nodes or greedy or reconstruct_text or factor_occ_freq) {
        tIDX index;
        string int_idx_file_name = tIDX::get_int_idx_filename(input_file_name.c_str(), b, output_dir.c_str());
        if (util::load_from_file(index, int_idx_file_name.c_str())) {
            index.set_file_name(input_file_name);
            index.set_output_dir(output_dir);
            if (output_mem_info) {
                util::write_structure<JSON_FORMAT>(index, std::cout);
            } else if (output_statistics) {
                index.statistics();
            } else if (output_Hk) {
                index.text_statistics(max_k);
            } else if (output_trie_nodes) {
                index.trie_nodes();
            } else if (greedy) {
				bit_vector factor_borders;
                size_type factors = index.greedy_parse(tmp_file_dir, factor_borders);
                std::cout << "factors = "<<factors << std::endl;
                std::cout << "avg factor length = "<<((double)index.size())/factors<< std::endl;
				std::cout << "k="<< index.k-1 << std::endl;
                int bpf = bit_magic::l1BP(index.k-1)+1;
                std::cout << "bit_magic::l1BP(index.k-1)+1="<< bpf << std::endl;
				std::cout << "lz_width="<<(int)index.lz_width<<std::endl;
                std::cout << "factorization size in MB:"<< ((double)factors*bpf)/(8*(1<<20)) << std::endl;
            } else if (reconstruct_text){
				index.reconstruct_text();
			}else if (factor_occ_freq){
				index.factor_frequency();
			}
        } else {
            std::cerr << "ERROR: could not open file "<< int_idx_file_name << endl;
            return 1;
        }
        return 0;
    }

    if (generate_patterns) {
        if (pattern_file_name == "") {
            std::string path = get_output_dir(input_file_name.c_str(), output_dir.c_str());
            pattern_file_name = path + "/" + util::basename(input_file_name) + "." + util::to_string(plen) +
                                "." + util::to_string(pno) +
                                "." + util::to_string(swapno) +
                                "." + util::to_string(pmin_occ) +
                                "." + util::to_string(pmax_occ) +
                                ".pattern";
        }
        cout << "-- start pattern generation -- " << endl;
        if (pmin_occ == pmax_occ and pmin_occ == 0) {  // if the occurrence restriction is not initialized
            pattern_file pf(pattern_file_name.c_str());
            pf.generate(input_file_name.c_str(), pno, plen, swapno);
        } else {
            std::string path = get_output_dir(input_file_name.c_str(), output_dir.c_str())+"/"+util::basename(input_file_name.c_str());
            string tmp_fwd_cst_file_name = path + "." +TMP_CST_SUFFIX;
            rosa<>::tCst cst;
            if (!util::load_from_file(cst, tmp_fwd_cst_file_name.c_str())) {
                std::cerr << "ERROR: could not open the compressed suffix tree file" << std::endl;
            } else {
                pattern_file pf(pattern_file_name.c_str());
                pf.generate_restricted(cst, input_file_name.c_str(), pno, plen,
                                       pmin_occ, pmax_occ);
            }
        }
        cout << "-- finished pattern generation -- " << endl;
        cout << "pattern stored in file: " <<pattern_file_name << endl;
        return 0;
    }

    if (run_benchmark or run_benchmark_int or run_benchmark_ext) {
        if (pattern_file_name == "") {
            std::string path = get_output_dir(input_file_name.c_str(), output_dir.c_str());
            pattern_file_name =  path + "." + util::to_string(plen) +
                                 "." + util::to_string(pno) +
                                 "." + util::to_string(swapno) + ".pattern";
        }
        std::cerr << "run benchmark" << std::endl;
        std::cerr << "pattern_file_name " << pattern_file_name << std::endl;
        benchmark<tIDX>(input_file_name.c_str(), b, pattern_file_name.c_str(), tmp_file_dir.c_str(),
                        output_dir.c_str(), repeated_in_memory_search,
                        run_benchmark_int, run_benchmark_ext);
    }
}
