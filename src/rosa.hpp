/*! This file contains the class for the external memory index rosa
 */
#ifndef ROSA_INCLUDED
#define ROSA_INCLUDED

#include "rosa_helping_structures.hpp" // for the construction
#include "rosa_helping_functions.hpp"  // for the construction and s
#include "rosa_tikz.hpp"			   // for tikz output

#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_io_wrappers.hpp>	// for encoding/decoding  
#include <sdsl/suffixtrees.hpp>		// for construction	
#include <sdsl/csa_wt.hpp>			// for default template parameters
#include <sdsl/rank_support.hpp>	// for default template parameters 
#include <sdsl/select_support.hpp>	// for default template parameters
#include <sdsl/wt_huff.hpp>   // for default template parameters
#include <sdsl/util.hpp>      // for load and serialize
#include <sdsl/testutils.hpp> // for file

#include <iostream>
#include <cstdlib> // for atoi
#include <cstring> // for strncmp
#include <string>  // for file_names
#include <vector>  // for construction
#include <utility> // for pair
#include <queue>   // for construction
#include <iomanip> // for setw
#include <map>	   // for a file map
#include <stack>

using namespace sdsl;
using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::vector;
using std::ifstream;
using std::map;
using std::ostream;

const char KEY_GREEDY_FACTOR[] = "greedy_factor";

#ifndef LCP_WRAP
#define LCP_WRAP 0 //VBYTE
#endif

// forward declaration of the rosa class with default template parameters
template<class BitVector = bit_vector // for bl and bf
,class RankSupport = typename BitVector::rank_1_type // for bl and bf
,class SelectSupport = typename BitVector::select_1_type // for bf
,class WaveletTree = wt_huff<bit_vector,
rank_support_v5<>,
select_support_mcl<1>,
select_support_mcl<0> > // for pruned BWT
#if LCP_WRAP == 0
,class LcpSerializeWrapper = int_vector_serialize_vbyte_wrapper<> // int_vector_serialize_wrapper<>
,class LcpLoadWrapper = int_vector_load_vbyte_wrapper<>  // int_vector_load_wrapper<>
#endif
#if LCP_WRAP == 1
,class LcpSerializeWrapper = int_vector_serialize_vlen_wrapper<> // int_vector_serialize_wrapper<>
,class LcpLoadWrapper = int_vector_load_vlen_wrapper<>  // int_vector_load_wrapper<>
#endif
>
class rosa;

#define TMP_CST_SUFFIX "tmp_fwd_cstX"
#define TMP_FWD_CSA_SUFFIX "tmp_fwd_csaX"
#define TMP_BWD_CSA_SUFFIX "tmp_bwd_csaX"


typedef bit_vector::size_type size_type;
typedef std::pair<size_type, size_type> tPII;
typedef vector< tPII > tVPII;
typedef std::queue< tPII > tQPII;
typedef std::queue< size_type > tQI;
typedef std::pair<unsigned char, size_type> tPCI;
typedef std::priority_queue<tPCI, vector<tPCI>, std::greater<tPCI> > tPQPCI;

template<class BitVector
,class RankSupport
,class SelectSupport
,class WaveletTree
,class LcpSerializeWrapper
,class LcpLoadWrapper
>
class rosa {
    public:
        typedef int_vector<>::size_type 									size_type;
        typedef BitVector 													bit_vector_type;
        typedef RankSupport 												rank_support_type;
        typedef SelectSupport 												select_support_type;
        typedef WaveletTree 												wavelet_tree_type;
        typedef csa_wt<wt_huff<rrr_vector<63> >,64000, 64000> 				tCsa;
//        typedef cst_sada<csa_wt<wt_rlmn<>,8, 256>, lcp_support_tree2<8> > 	tCst;
        typedef select_support_mcl<>  										bm_select_1_type;
        typedef select_support_mcl<0> 										bm_select_0_type;
        typedef rank_support_v5<10,2> 										bm_rank10_type;
    private:
        size_type				m_n;  // original text length
        size_type 				m_b;  // block size
        size_type				m_k;  // number of external blocks
		uint8_t					m_lz_width; // bit-width of LZ factors
		size_type				m_lz_size;  // #of factors in the LZ factorization 
        bit_vector_type 		m_bl; // indicates the BWT entries which corresponds to a representative suffix
        bit_vector_type 		m_bf; // indicates the representative suffixes in SA
        rank_support_type 		m_bl_rank; // rank support for m_bl
        rank_support_type 		m_bf_rank; // rank support for m_bf
        select_support_type 	m_bf_select; // select support for m_bf
        select_support_type     m_bl_select; // select support for m_bl
        wavelet_tree_type		m_wt; // wavelet tree for pruned BWT
        int_vector<64>			m_cC; // contains for each character c the m_rank_bf(C[c])
        int_vector<8>			m_char2comp;
        int_vector<8>			m_comp2char;
        bit_vector				m_bm; // bit sequence to map from backward intervals to blocks in the external part
        bm_select_1_type		m_bm_1_select; //
        bm_select_0_type		m_bm_0_select; //
        bm_rank10_type			m_bm_10_rank; //
        int_vector<>			m_min_depth;
        int_vector<>			m_pointer; // address of block [sp..ep] on disk,
        // if ep-sp > 0 or SA[sp] if sp=ep.

        string					m_file_name;  // file name of the supported text
        string					m_output_dir; // output directory

        mutable ifstream		m_text;	      // stream to the text
		mutable ifstream		m_glz_text;   // stream to the factored text 
        mutable ifstream		m_ext_idx;    // stream to the external memory part
        unsigned char*			m_buf;		  // buffer for the text read from disk to match against the pattern
		const uint64_t*      	m_buf_lz;
        size_type  				m_buf_size;   // buffer

#ifdef OUTPUT_STATS
        mutable size_type m_count_disk_access; 		// see corresponding public member
        mutable size_type m_count_gap_disk_access; 	// see corresponding public member
        mutable size_type m_count_int_steps;   		// see corresponding public member
        mutable size_type m_count_int_match;   		// see corresponding public member
        mutable size_type m_count_queries;	   		// see corresponding public member
        mutable size_type m_count_block_length;     // see corresponding public member
#endif

        //! Internal helper class for an external block
        /*  TODO: possible optimizations:
         *        * don't store the triple of the corresponding irreducible block
         *          and use delta_x=delta_d=0 if we don't find its backward id
         *        * subtract and add d from the lcp values
         */
        class disk_block {
                int_vector<> m_header;	// Each entry contains an encoded header triple.
                uint8_t m_width_bwd_id;	// Bit width of the maximal value of a backward id in the header.
                uint8_t m_width_delta_x;// Bit width of the maximal delta x in the header.
                int_vector<>  m_lcp; 	// Array for the relative bit-Lcp values.
                int_vector<>  m_sa;		// Array for the suffix array values.
                bit_vector    m_bp_ct;  // Bitvector for the balanced parentheses of the Cartesian Tree of the LCP Array
            public:

                const int_vector<>&  header;// const reference to the header
                const int_vector<>&  lcp;	// const reference to the lcp values
                const int_vector<>&  sa;	// const reference to the suffix array values
                const bit_vector&    bp_ct; // const reference to the balanced parentheses

                disk_block():header(m_header), lcp(m_lcp), sa(m_sa), bp_ct(m_bp_ct) {};

                size_type header_size_in_bytes() {
                    return sizeof(m_width_bwd_id) + sizeof(m_width_delta_x) +
                           util::get_size_in_bytes(m_header);
                }

                void set_header(const std::vector<header_item>& vh) {
                    size_type max_bwd_id = 0;
                    size_type max_delta_x = 0;
                    size_type max_delta_d = 0;
                    for (size_t i=0; i < vh.size(); ++i) {
                        max_bwd_id = std::max(max_bwd_id, vh[i].bwd_id);
                        max_delta_x = std::max(max_delta_x, vh[i].delta_x);
                        max_delta_d = std::max(max_delta_d, vh[i].delta_d);
                    }
                    // as vh is asserted to be sorted by bwd_id we can take the last element
                    m_width_bwd_id 			= bit_magic::l1BP(max_bwd_id)+1;
                    m_width_delta_x 		= bit_magic::l1BP(max_delta_x) + 1;
                    uint8_t width_delta_d 	= bit_magic::l1BP(max_delta_d) + 1;
                    uint8_t total_width		= m_width_bwd_id + m_width_delta_x + width_delta_d;
                    if (total_width > 64) {
                        std::cerr << "ERROR: header triple is to big for a 64bit integer / its " <<
                                  total_width << " bits " << std::endl;
                    }
                    m_header.set_int_width(total_width);
                    m_header.resize(vh.size());
                    for (size_t i=0; i < vh.size(); ++i) {
                        m_header[i] = vh[i].bwd_id;
                        m_header[i] = m_header[i] | (vh[i].delta_x << m_width_bwd_id);
                        m_header[i] = m_header[i] | (vh[i].delta_d << (m_width_bwd_id+m_width_delta_x));
                    }
                }

                // block boundaries [lb..rb] (both inclusive)
                /*!
                 *
                 * \pre rb>lb
                 */
                inline void set_content(const int_vector<64> &sa_buf, const int_vector<64> &lcp_buf, const unsigned char* text, size_type block_len) {
                    size_type max_sa = 0;
                    size_type max_lcp = 0;
                    m_lcp 	= int_vector<>(block_len);
                    m_sa  	= int_vector<>(block_len);
                    for (size_type i=1; i <= block_len; ++i) { // copy SA and LCP values of the block
                        m_sa[i-1] = sa_buf[i];
                    }
//					cout<<"entered sa values"<<endl;
                    calculate_bit_lcp(m_lcp, sa_buf, lcp_buf, text);
/*					
                    if (util::verbose) {
						std::cout<<"block_len="<<block_len<<std::endl;
                        std::cout << "LCP ="; for (size_type i=0; i<m_lcp.size(); ++i) {
                            std::cout << " " << m_lcp[i] << "("<< lcp_buf[i+1] <<") ";
                        }; std::cout << std::endl;
                    }
*/					
                    calculate_bp(m_bp_ct, m_lcp);
                    util::bit_compress(m_sa);			   // bit compress suffix array values
                    calculate_relative_bit_lcp(m_lcp, 0, m_lcp.size(), 0);
/*
                    if (util::verbose) {
                        std::cout << "[" << lb << "," << rb << "]" << std::endl;
                        std::cout << "SA ="; for (size_type i=0; i<m_sa.size(); ++i) {
                            std::cout << " " << m_sa[i];
                        }; std::cout << std::endl;
                        std::cout << "LCP ="; for (size_type i=0; i<m_lcp.size(); ++i) {
                            std::cout << " " << m_lcp[i];
                        }; std::cout << std::endl;
                    }
*/					
                }

				template<class tRank>
				void replace_pointers(const tRank& factor_borders_rank){
					for (size_type i=0; i<m_sa.size(); ++i){
						m_sa[i] = factor_borders_rank(m_sa[i]);
					}
					// TODO: use bit_compress and adjust pointers of in-memory structure
					// util::bit_compress(m_sa);
				}

                template<class rac>
                void calculate_bp(bit_vector& bp_ct, const rac& lcp) {
                    size_type N = lcp.size();
                    // initialize bitvectors
                    util::assign(bp_ct, bit_vector(2*N, 0));
                    // calculate content
                    stack<size_type> prev;     // stack for the previous seen LCP values
                    for (size_type i=0, idx=0, x; i< N; ++i) {
                        while (!prev.empty() and (x=lcp[prev.top()]) > lcp[i]) {
                            prev.pop();
                            idx++; // move forward in bp_ct == write a closing parenthesis
                        }
                        m_bp_ct[idx++] = 1; // write opening parenthesis for element i
                        prev.push(i);
                    }
                }

                template<class rac>
                void calculate_bit_lcp(rac& lcp, const int_vector<64>& sa_buf, const int_vector<64>& lcp_buf, const unsigned char* text) {
                    for (size_type i=0; i < lcp.size(); ++i) {
                        unsigned char c1 = text[sa_buf[i]   + lcp_buf[i+1]];
                        unsigned char c2 = text[sa_buf[i+1] + lcp_buf[i+1]];
                        lcp[i] = lcp_buf[i+1]*8 + 7-bit_magic::l1BP(c1^c2);
                        c1 = c2;
                    }
                }

                // TODO: till now, this is quadratic in the maximal block size
                // the complexity should be linear in a final implementation
                template<class rac>
                void calculate_relative_bit_lcp(rac& lcp, size_type lb, size_type rb, size_type min_lcp) {
                    if (lb < rb) {
                        size_type min_idx = lb;
                        for (size_type i = lb+1; i < rb; ++i) {
                            if (lcp[i] < lcp[min_idx]) {
                                min_idx = i;
                            }
                        }
                        calculate_relative_bit_lcp(lcp, lb, min_idx, lcp[min_idx]); // left sub-array
                        calculate_relative_bit_lcp(lcp, min_idx+1, rb, lcp[min_idx]);
                        lcp[min_idx] = lcp[min_idx] - min_lcp;
                    }
                }



                //! Decode the i-th triple of the header
                /*! \param i		Index of the header-triple.
                 *  \param bwd_id	Output reference for the backward id.
                 *  \param delta_x	Output reference for delta_x.
                 *	\param delta_d	Output reference for delta_d.
                 */
                void decode_header_triple(const size_type i, size_type& bwd_id, size_type& delta_x, size_type& delta_d)const {
                    bwd_id		= m_header[i] & bit_magic::Li1Mask[m_width_bwd_id];
                    uint64_t w 	= m_header[i] >> m_width_bwd_id;
                    delta_x 	= w & bit_magic::Li1Mask[m_width_delta_x];
                    w >>= m_width_delta_x;
                    delta_d = w;
                }

                //! Get the corresponding delta_x and delta_d values for a given backward id
                /*! \param bwd_id	Backward id of the block (information from the in-memory part).
                 *  \param delta_x	Output reference for delta_x.
                 *	\param delta_d	Output reference for delta_d.
                 */
                // TODO: using binary search may be and optimisation if header.size() is big
                void get_delta_x_and_d(const size_type bwd_id, size_type& delta_x, size_type& delta_d) {
                    for (size_type i=0; i < m_header.size(); ++i) {
                        if ((m_header[i] & bit_magic::Li1Mask[m_width_bwd_id]) == bwd_id) {
                            uint64_t w = m_header[i] >> m_width_bwd_id;
                            delta_x = w & bit_magic::Li1Mask[m_width_delta_x];
                            w >>= m_width_delta_x;
                            delta_d = w;
                            return;
                        }
                    }
                    throw std::logic_error("bwd id not found in header");
                }

                //! Get the number of the block
                size_type size()const {
                    return m_sa.size();
                }

                //! Write the block to the output stream
                size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
                    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
                    size_type written_bytes = 0;
                    written_bytes += util::write_member(m_width_bwd_id, out, child, "width_bwd_id");
                    written_bytes += util::write_member(m_width_delta_x, out, child, "width_delta_x");
                    written_bytes += m_header.serialize(out, child, "header");
                    written_bytes += m_bp_ct.serialize(out, child, "bp_ct");
                    written_bytes += LcpSerializeWrapper(m_lcp).serialize(out, child, "lcp");
                    written_bytes += m_sa.serialize(out, child, "sa");
                    structure_tree::add_size(child, written_bytes);
                    return written_bytes;
                }

                //! Load the block from the input stream
                void load(std::istream& in) {
                    util::read_member(m_width_bwd_id, in);
                    util::read_member(m_width_delta_x, in);
                    m_header.load(in);
                    m_bp_ct.load(in);
                    LcpLoadWrapper(m_lcp).load(in);
                    m_sa.load(in);
                }
        };

        // Private helper methods

        //! Open the stream of the text and the external index part
        /*! If the streams are already associated with a file they are closed
         *  before we open the new files.
         */
        void open_streams() {
			open_stream(m_text, m_file_name.c_str());
			open_stream(m_glz_text, get_factorization_filename().c_str());
			{
				int_vector_file_buffer<> glz_buffer(get_factorization_filename().c_str());
				m_lz_width = glz_buffer.int_width;
				m_lz_size  = glz_buffer.int_vector_size;
				if (util::verbose){
					cout<<"m_lz_width = "<<(int)m_lz_width<<" m_lz_size = "<<m_lz_size<<endl;
				}
			}
			open_stream(m_ext_idx, get_ext_idx_filename().c_str());
        }

		void open_stream(ifstream &stream, const char* file_name){
	        close_stream_if_open(stream);
            stream.open(file_name);
            if (!stream) {
                std::cerr << "Error: Could not open file: " << file_name << std::endl;
            } else {
                if (util::verbose) std::cerr << "Opened file " << file_name << "\n";
            }	
		}

        //! Wrapper for the seekg method of ifstream that increases the disk access counter
        void seekg(ifstream& in, size_type block_addr, bool count=true)const {
            in.seekg(block_addr, std::ios::beg);
            if (count) {
#ifdef OUTPUT_STATS
                ++m_count_disk_access;
#endif
            }
        }

    public:
        const bit_vector_type& bf;				//!< Bit vector which indicates intervals in the backward index.
        const bit_vector_type& bl; 				//!< Bit vector which contains the LF-permuted content of bf.
        const wavelet_tree_type& wt;			//!< Wavelet tree holding the condensed BWT in the backward index.
        const rank_support_type& bl_rank;		//!< Rank support for bl.
        const rank_support_type& bf_rank;		//!< Rank support for bf.
        const select_support_type& bf_select;	//!< Select support for bf.
        const select_support_type& bl_select;	//!< Select support for bl.
        const bit_vector& bm;					//!< Bit vector to map from intervals in the backward index to blocks in the external part.
        const bm_select_1_type& 	bm_1_select;//!< Select support for ones in bm.
        const bm_select_0_type& 	bm_0_select;//!< Select support for zeros in bm.
        const bm_rank10_type& bm_10_rank;		//!< Rank support for the bit-pattern 01 in bm.
        const int_vector<>& min_depth;			//!<
        const int_vector<>& pointer;			//!< Array of pointers into the external structures (blocks or suffix array).
        const string& file_name;				//!< File name of the original text string.
        const string& output_dir;				//!< Directory where the output is stored.
		const size_type& k;						//!< Number of block prefixes.
		const uint8_t& lz_width;				//!< Bit-width of the LZ factors
#ifdef OUTPUT_STATS
        const size_type& count_disk_access; 	//!< Counter for potential disk accesses.
        const size_type& count_gap_disk_access; //!< Counter for potential disk accesses caused by gaps in the fringe.
        const size_type& count_int_steps; 		//!< Counter for matches that can be answered with the in-memory part of the data structure.
        const size_type& count_int_match; 		//!< Counter for matches that can be answered with the in-memory part of the data structure.
        const size_type& count_queries;			//!< Counter for the queries.
        const size_type& count_block_length;	//!< The sum of the length of all fetched blocks in elements.
#endif

        ~rosa() {
            m_text.close();    // close stream to the text
			m_glz_text.close(); // close stream to the factorization
            m_ext_idx.close(); // close stream to the external part
            if ( NULL != m_buf ) delete [] m_buf;
			if ( NULL != m_buf_lz) delete [] m_buf_lz;
        }

        //! Get the name of the file where the external part of the index is stored.
        static string get_ext_idx_filename(const char* file_name, size_type b,  const char* output_dir=NULL) {
            return get_output_dir(file_name, output_dir) + "/" + util::basename(file_name)
                   + "." + util::to_string(b) + ".2." + SDSL_XSTR(LCP_WRAP)  + ".bt2.ext_idx";
        }

        //! Get the name of the file where the external part of the index is stored.
        string get_ext_idx_filename() {
            return get_ext_idx_filename(m_file_name.c_str(), m_b, m_output_dir.c_str());
        }

        //! Get the name of the file where the external part of the index is stored.
        string get_tmp_ext_idx_filename() {
            return get_ext_idx_filename(m_file_name.c_str(), m_b, m_output_dir.c_str())+".tmp";
		}
        

        //! Get the name of the file where the external part of the index is stored.
        string get_factorization_filename() {
			return get_output_dir(m_file_name.c_str(), m_output_dir.c_str()) + "/" + util::basename(m_file_name.c_str())
			       +"."+ util::to_string(m_b)	+ ".2.glz";
        }


        //! Get the name of the file where the in-memory part of the index is stored.
        static string get_int_idx_filename(const char* file_name, size_type b, const char* output_dir=NULL) {
            return get_output_dir(file_name, output_dir) + "/" + util::basename(file_name)
                   + "." + util::to_string(b) + ".2." + SDSL_XSTR(LCP_WRAP) + "." +INDEX_SUF+".int_idx";
        }

        //! Get the name of the file where the in-memory part of the index is stored.
        string get_int_idx_filename() {
            return get_int_idx_filename(m_file_name.c_str(), m_b, m_output_dir.c_str());
        }

        //! Remove the temporary files which are created during the construction.
        static void remove_tmp_files(const char* file_name, const char* output_dir) {
            string path = get_output_dir(file_name, output_dir) + "/" + util::basename(file_name);
            std::remove((path + "." + TMP_FWD_CSA_SUFFIX).c_str());
            std::remove((path + "." + TMP_BWD_CSA_SUFFIX).c_str());
        }

        /*!
         *  \param cst			CST which should be loaded or constructed.
         *	\param file_name	Location of the CST on disk.
         *	\param tmp_dir 		Location of the temporary directory for the construction.
         *  \param delete_tmp   Boolean flag, if the generated files uncompressed SA, LCP, and so on should be deleted.
         */
        void construct_or_load_fwd_csa(tCsa& csa, string file_name, string tmp_dir, bool delete_tmp) {
            ifstream tmp_fwd_csa_stream(file_name.c_str());
            if (!tmp_fwd_csa_stream) {
                if (util::verbose) {
                    cout<<"m_file_name="<<m_file_name<<endl;
                    cout<<"tmp_dir="<<tmp_dir<<endl;
                    cout<<"id="<<util::basename(m_file_name)<<endl;
                }
                cache_config config(delete_tmp, tmp_dir, util::basename(m_file_name));
                construct(csa, m_file_name.c_str(), config, 1);
                util::store_to_file(csa, file_name.c_str());
            } else {
                tmp_fwd_csa_stream.close();
                if (util::verbose) {
                    cout<< "load stored fwd_cst form file "<< file_name << endl;
                }
                util::load_from_file(csa, file_name.c_str());
            }
        }

        string get_rev_file_name(const string& file_name) {
            return file_name+"_rev";
        }

        /*!
         *  \param bwd_csa		CSA which should be loaded or constructed.
         *  \param file_name	Location of the CSA on disk.
         *	\param tmp_dir 		Location of the temporary directory for the construction.
         *  \param delete_tmp   Boolean flag, if the generated files uncompressed SA, LCP, and so on should be deleted.
         */
        void construct_or_load_bwd_csa(tCsa& bwd_csa, string file_name, string tmp_dir, bool delete_tmp) {
            ifstream tmp_bwd_csa_stream(file_name.c_str());
            if (!tmp_bwd_csa_stream) {
                string rev_file_name = get_rev_file_name(m_file_name);
                ifstream rev_text_stream(rev_file_name.c_str());
                if (!rev_text_stream) {
                    int_vector<8> text;
                    util::load_vector_from_file(text, m_file_name.c_str(), 1);
                    for (size_type i=0, j=text.size()-1; i < j and j < text.size(); ++i, --j) {
                        std::swap(text[i], text[j]);
                    }
                    util::store_to_plain_array<uint8_t>(text, rev_file_name.c_str());
                }
                cache_config config(delete_tmp, tmp_dir, (util::basename(rev_file_name)));
                construct(bwd_csa, rev_file_name.c_str(), config, 1);
                util::store_to_file(bwd_csa, file_name.c_str());
            } else {
                if (util::verbose) {
                    std::cout<< "load stored bwd_csa form file "<< file_name << std::endl;
                }
                util::load_from_file(bwd_csa, file_name.c_str());
            }
        }

        void calculate_bl_and_bf_and_mapping(const tCsa& csa, vector<block_info>& map_info) {
            bit_vector bl = bit_vector(csa.size());
            bit_vector bf = bit_vector(csa.size() + 1);   // add one bit for technical reasons
            tQPII q; q.push(tPII(0, csa.size()));     // init queue with root interval [0..csa.size()]
            tQI q_fp; q_fp.push(0); 				    // init queue for the forward position of the left bound of the interval
            tQI q_depth; q_depth.push(0); 			  // init queue for depth of the interval
            size_type k;                                  // number of different symbols in interval [lb..rb]
            vector<unsigned char> cs(256);   			  // cs[0..k-1] contains the different symbols of [lb..rb]
            vector<size_type> rank_lb(256), rank_rb(256); //rank_lb[i] = rank(lb, cs[i]) and rank_rb[i] = rank(rb, cs[i])

            while (!q.empty()) {
                size_type lb = q.front().first;
                size_type rb = q.front().second;
                size_type lb_forward = q_fp.front();
                q.pop(); q_fp.pop();
                size_type depth = q_depth.front(); q_depth.pop();
                if (!bf[lb]) {
                    bf[lb] = 1;
                    bl[ csa.psi[lb] ] = 1;
                }
                if (!bf[rb]) {
                    bf[rb] = 1;  // rb has to be marked
                    if (rb < csa.size()) {
                        bl[ csa.psi[rb] ] = 1;
                    }
                }
                if (rb-lb > m_b) { // if
                    // get all different symbols and ranks in the interval [lb..rb]
                    csa.wavelet_tree.interval_symbols(lb, rb, k, cs, rank_lb, rank_rb);
                    // sort the different symbols lexicographically
                    tPQPCI pq; // by using a priority queue
                    for (size_type i=0; i<k; ++i) { // push character and its index to pq to sort it
                        pq.push(tPCI(cs[i], i));    // this step would be not necessary if we use a sorted WT
                    }
                    for (size_type i=0, j, lex_smaller=0; i<k; ++i) {
                        j = pq.top().second; pq.pop();
                        size_type c_begin = csa.C[ csa.char2comp[cs[j]] ];
                        size_type lb_new = c_begin + rank_lb[j], rb_new =c_begin + rank_rb[j];
                        q.push(tPII(lb_new, rb_new));
                        q_fp.push(lb_forward + lex_smaller);
                        q_depth.push(depth+1);
                        lex_smaller += rank_rb[j]-rank_lb[j];
                    }
                } else {
                    map_info.push_back(block_info(lb, depth, lb_forward, rb-lb));
                }
            }
            sort(map_info.begin(),map_info.end());
            util::assign(m_bl, bl); if (util::verbose) std::cout << "m_bl was assigned.\n";
            util::assign(m_bf, bf); if (util::verbose) std::cout << "m_bf was assigned.\n";
            util::init_support(m_bf_rank, &m_bf); 		if (util::verbose) std::cout << "m_bf_rank was assigned.\n";
            util::init_support(m_bf_select, &m_bf);		if (util::verbose) std::cout << "m_bf_select was assigned.\n";
            util::init_support(m_bl_select, &m_bl);		if (util::verbose) std::cout << "m_bl_select was assigned.\n";
        }

        void calculate_bm_and_min_depth(const vector<block_info>& map_info, const size_type f_k) {
			if (util::verbose) cout<<"m_k="<<m_k<<endl;
			if (util::verbose) cout<<"f_k="<<f_k<<endl;
			if (util::verbose) cout<<"map_info.size()="<<map_info.size()<<endl;
            m_bm = bit_vector(f_k     + m_bf_rank(m_bf.size()), 0);
//                                 ^^^^^^^^ ^^^^^^^^^^^^^^^^^^^^^
//                                 # of 0s        # of 1s
            size_type max_min_depth = 0;
            // construct bm and get information about the max min_depth
            for (size_type i=1, j=0, p=0, end=m_bf_rank(m_bf.size()); i < end; ++i) {
                size_type lb = m_bf_select(i);
                if (j < map_info.size() and map_info[j].bwd_lb == lb) {
                    if (map_info[j].depth > max_min_depth)
                        max_min_depth = map_info[j].depth;
                }
                while (j < map_info.size() and map_info[j].bwd_lb == lb) {
                    ++p; // write zeros in m_map_block
                    ++j;
                }
                m_bm[p++] = 1; // mark the i-th entries in m_bm
            }
            if (util::verbose) cout << "m_mb was calculated.\n";
            util::init_support(m_bm_1_select, &m_bm); if (util::verbose) cout<<"m_bm_1_select was assigned.\n";
            util::init_support(m_bm_0_select, &m_bm); if (util::verbose) cout<<"m_bm_0_select was assigned.\n";
            util::init_support(m_bm_10_rank, &m_bm); if (util::verbose) cout<<"m_bm_10_rank was assigned.\n";
            m_min_depth.set_int_width(bit_magic::l1BP(max_min_depth)+1);
            m_min_depth.resize(m_bm_10_rank(m_bm.size())+1);
            // insert min_depth information
            for (size_type i=1, j=0, k=0, end=m_bf_rank(m_bf.size()); i < end; ++i) {
                size_type lb = m_bf_select(i);
                bool non_empty = false;
                size_type min_depth = 0;
                while (j < map_info.size() and map_info[j].bwd_lb == lb) {
                    if (!non_empty) {
                        min_depth = map_info[j].depth;
                        non_empty = true;
                    }
                    ++j;
                }
                if (non_empty) {
                    m_min_depth[k++] = min_depth;
                }
            }
            if (util::verbose) cout << "m_min_depth was calculated.\n";
        }

        void calculate_bwd_id_and_fill_singleton_pointers(const char* sa_file, const tCsa &csa,
                const vector<block_info>& map_info,
                const bit_vector& fwd_bf,
                vector<block_node>& v_block,
				bit_vector& is_singleton
				) {
            rank_support_v5<> fwd_bf_rank(&fwd_bf);
			int_vector<> block_sa(fwd_bf_rank(fwd_bf.size()), 0, bit_magic::l1BP(fwd_bf.size())+1);
			int_vector_file_buffer<> sa_buf(sa_file);
			for (size_type i=0,r=0,r_sum=0,idx=0; i < sa_buf.int_vector_size;) { 
				for (; i < r+r_sum; ++i) {
					if ( fwd_bf[i] ){
						block_sa[idx++] = sa_buf[i-r_sum];
					}
				}
				r_sum += r; r = sa_buf.load_next_block();
			}

            for (size_t bwd_id=0; bwd_id < map_info.size(); ++bwd_id) {
                size_type fwd_lb = map_info[bwd_id].fwd_lb;
                size_type fwd_id = fwd_bf_rank(fwd_lb+1)-1; // calculate corresponding forward id
                v_block[fwd_id].bwd_id = bwd_id;
                if (map_info[bwd_id].size == 1) {
					m_pointer[bwd_id] = block_sa[fwd_id];
					is_singleton[bwd_id] = 1;
                }
            }
            if (util::verbose) cout << "SA pointer for singleton block were inserted.\n";
            if (util::verbose) cout << "singleton block were marked in bit_vector is_singleton.\n";
            if (util::verbose) cout << "fwd_id <-> bwd_id mapping was calculated.\n";
        }

        void calculate_headers(const vector<block_node>& v_block,
                               const vector<block_info>& map_info,
                               vector<vector<header_item> >& header_of_external_block) {
            for (size_t fwd_id=0; fwd_id < v_block.size(); ++fwd_id) {
                size_type bwd_id = v_block[fwd_id].bwd_id;
                if (map_info[bwd_id].size > 1) {
                    size_type dest_block = v_block[fwd_id].dest_block;
                    size_type delta_x = v_block[fwd_id].delta_x;
                    size_type delta_d = v_block[fwd_id].delta_d;
                    header_of_external_block[dest_block].push_back(header_item(bwd_id,delta_x,delta_d));
                }
            }
            if (util::verbose) cout << "Header for external block calculated.\n";
        }

		uint64_t get_next_int(ifstream &stream, uint8_t width, uint8_t& offset, uint64_t& cached_word){
			uint64_t res = 0;
			if ( offset == 0 ){
				stream.read((char*)&cached_word, sizeof(uint64_t));
				res = bit_magic::read_int(&cached_word, offset, width);
				offset = (offset+width)&0x3F;
			}else if ( offset + width <= 64 ){
				res = bit_magic::read_int(&cached_word, offset, width);
				offset = (offset+width)&0x3F;
			}else{ // handle case where integers spans across borders
				res = (cached_word)>>offset;
				uint8_t  bit_of_old_word = 64-offset;
				stream.read((char*)&cached_word, sizeof(uint64_t));
				offset = (offset+width)&0x3F;
				res |= (((cached_word)&bit_magic::Li1Mask[offset])<<bit_of_old_word);
			}
			return res;
		}

		/*!
		 * \par Working set size
		 *      text size + constant size buffers + O(m_b) words
		 */
        void write_external_blocks_and_pointers(const char* sa_file,
												const char* lcp_file,
                                                const bit_vector& fwd_bf,
                                                const vector<block_node>& v_block,
                                                vector<vector<header_item> >& header_of_external_block) {
            size_type total_header_in_bytes = 0;
            size_type ext_idx_size_in_bytes = 0;
            std::ofstream ext_idx_out(get_tmp_ext_idx_filename().c_str(), std::ios_base::trunc);
            if (ext_idx_out) {
//TODO: test if the size of the output buffer chances the performance
                select_support_mcl<> fwd_bf_select(&fwd_bf);
                vector<size_type> block_addr(m_k+1, 0);
                char* text = NULL;
                file::read_text(m_file_name.c_str(), text);

				ifstream sa_stream(sa_file);
				ifstream lcp_stream(lcp_file);
				uint64_t sa_size, lcp_size, sa_word=0, lcp_word=0, sa_idx=0, lcp_idx=0;
				uint8_t sa_width, lcp_width, sa_off=0, lcp_off=0;
				sa_stream.read((char*)&sa_size,sizeof(sa_size)); lcp_stream.read((char*)&lcp_size,sizeof(lcp_size));
				sa_stream.read((char*)&sa_width,sizeof(sa_width)); lcp_stream.read((char*)&lcp_width,sizeof(lcp_width));
				int_vector<64> sa_buf(m_b+1, 0);
				int_vector<64> lcp_buf(m_b+1, 0);
//				int_vector<> lcp;
//				util::load_from_file(lcp, lcp_file);
				size_type sa_size_in_bytes = 0;
				size_type sa_size_in_bytes1 = 0;

                for (size_t fwd_id=0; fwd_id < header_of_external_block.size(); ++fwd_id) {
                    if (header_of_external_block[fwd_id].size() > 0) { // if it is an irreducible block
                        sort(header_of_external_block[fwd_id].begin(), header_of_external_block[fwd_id].end());
                        disk_block db;
                        db.set_header(header_of_external_block[fwd_id]);
						size_type lb = fwd_bf_select(fwd_id+1), rb = fwd_bf_select(fwd_id+2)-1;
						if ( lb==0 ){
							throw std::logic_error("lb==0: No sentinel character appended?");
						}
						while ( sa_idx < lb-1 ){
							get_next_int(sa_stream, sa_width, sa_off, sa_word);
							get_next_int(lcp_stream, lcp_width, lcp_off, lcp_word);
							++sa_idx; ++lcp_idx;
						}
						for (size_type i= (sa_idx==lb); sa_idx<=rb; ++i, ++sa_idx, ++lcp_idx){
							sa_buf[i]	=	get_next_int(sa_stream, sa_width, sa_off, sa_word);
							lcp_buf[i]	=	get_next_int(lcp_stream, lcp_width, lcp_off, lcp_word);
//							if ( lcp_buf[i] != lcp[lcp_idx] ){
//								cout<<"ERROR: i="<<i<<" lcp_idx="<<lcp_idx<<"   lcp_width="<<(int)lcp_width<<" lcp_off="<<(int)lcp_off<<endl;
//								cout<<"lcp_buf[i]="<<lcp_buf[i]<<" lcp[lcp_idx]="<<lcp[lcp_idx]<<endl;
//							}
						}
						size_type block_len = rb-lb+1;
//						cout<<"["<<lb<<","<<rb<<"]"<<endl;
						//             sa_buf and lcp_buf store block_len+1 integers 
                        db.set_content(sa_buf, lcp_buf, (const unsigned char*)text, block_len);
						sa_size_in_bytes1 += util::get_size_in_bytes(db.sa);

						sa_size_in_bytes += 9 + (((bit_magic::l1BP(block_len)+1)*block_len+63)/64)+8;
						sort(sa_buf.begin(), sa_buf.begin()+block_len);
						size_type sa_diff_bits = coder::elias_delta::encoding_length(sa_buf[0]);
						for (size_t i=1; i < block_len; ++i){ 
							sa_diff_bits += coder::elias_delta::encoding_length(sa_buf[i]-sa_buf[i-1]);
						}
						sa_size_in_bytes += ((sa_diff_bits+63)/64)*8;

                        block_addr[fwd_id] = ext_idx_size_in_bytes;
                        ext_idx_size_in_bytes += db.serialize(ext_idx_out);
                        total_header_in_bytes += db.header_size_in_bytes();
						sa_buf[0] = sa_buf[block_len];
						lcp_buf[0] = lcp_buf[block_len];
                    }
                }
                ext_idx_out.close();
                delete [] text;
                if (util::verbose)cout<<"# ext_idx_size_in_MB = "<<ext_idx_size_in_bytes/(1024.0*1024.0)<<"\n";
                if (util::verbose)cout<<"# ext_idx_header_size_in_MB = "<<total_header_in_bytes/(1024.0*1024)<<"\n";
				if (util::verbose)cout<<"# sa_size_in_MB = " << ((double)sa_size_in_bytes)/(1024*1024.0) << "\n";
				if (util::verbose)cout<<"# sa_size1_in_MB = " << ((double)sa_size_in_bytes1)/(1024*1024.0) << "\n";
                for (size_t fwd_idx = 0; fwd_idx < v_block.size(); ++fwd_idx) {
                    size_type lb = fwd_bf_select(fwd_idx+1);
                    if (!fwd_bf[lb+1]) {  // if not a singleton interval
                        size_type addr = block_addr[v_block[fwd_idx].dest_block];
                        m_pointer[ v_block[fwd_idx].bwd_id ] = addr;
                    }
                }
            } else {
                std::cerr << "ERROR: Could not open file " << get_tmp_ext_idx_filename() << endl;
            }
        }

        void init_bl_rank_and_cC(const tCsa& csa) {
            util::init_support(m_bl_rank, &m_bl);
            m_cC = int_vector<64>(csa.sigma+1, 0);  // condensed C
            for (size_type i=0; i <= csa.sigma; ++i) {
                m_cC[i] = m_bf_rank(csa.C[i]);
            }
        }

        void construct_condensed_bwt(const tCsa& csa, const string& tmp_dir) {
            string tmp_file_name =  tmp_dir+util::basename(m_file_name)+"_cBWT_"+util::to_string(util::get_pid())+"_"+util::to_string(util::get_id());
            size_type cn = m_bl_rank(m_n); // condensed n
            int_vector<8> temp_bwt(cn);
            for (size_type i=0,j=0; i<m_n; ++i) {
                if (m_bl[i]) {
                    temp_bwt[j++] = csa.bwt[i];
                }
            }
            util::store_to_file(temp_bwt, tmp_file_name.c_str());
            int_vector_file_buffer<8> temp_bwt_buf(tmp_file_name.c_str());
            util::assign(m_wt, wavelet_tree_type(temp_bwt_buf, cn));
            std::remove(tmp_file_name.c_str()); // remove file of BWT'
        }

		void output_tikz(){
            string base_name = util::basename(file_name)+"."+util::to_string(m_b);
            string bwd_csa_file_name = base_name + ".tikz." + TMP_BWD_CSA_SUFFIX;
			ofstream bwd_out((base_name+".bwd_idx.tex").c_str());
			tCsa bwd_csa;
            construct_or_load_bwd_csa(bwd_csa, bwd_csa_file_name, "./", true);
			std::vector<bwd_block_info> fwd_blocks_in_bwd;
			size_type max_bwd_id = m_k;
			for (size_type bwd_id=0; bwd_id<max_bwd_id; ++bwd_id){
				fwd_blocks_in_bwd.push_back(bwd_id_to_info(bwd_id));
			}
			write_tikz_output_bwd(bwd_out, bwd_csa, m_bf, m_bl, m_bm, m_b, fwd_blocks_in_bwd, m_min_depth);

			
			ofstream factor_out((base_name+".fac.tex").c_str());
				
			int_vector<> factorization;
			util::load_from_file(factorization, get_factorization_filename().c_str());
			write_tikz_array(factor_out, factorization, "facArray");
			int_vector<> factor_len(factorization.size(), 0);
			for (size_type i=0; i <factorization.size() ; ++i){
				string factor;
				extract_factor(factorization[i], factor);
				factor_len[i] = factor.size();
			}
			write_tikz_array(factor_out, factor_len, "facLen");
			string tmp = algorithm::extract(bwd_csa, 0, bwd_csa.size()-1);
			vector<string> text(tmp.size());
			for (int i=tmp.size()-2; i>=0; --i)
				text[tmp.size()-2-i] = util::to_latex_string((unsigned char)tmp[i]);
			text[tmp.size()-1] = util::to_latex_string((unsigned char)tmp[tmp.size()-1]);
			write_tikz_array(factor_out, text, "FwdText", true);
/*
			ofstream fwd_out((base_name+".fwd_idx.tex").c_str());
	        m_ext_idx.seekg(0, std::ios::end);
            std::streampos end = m_ext_idx.tellg(); // get the end position
            seekg(m_ext_idx, 0, false); // load the first block
            // iterate through all blocks
            while (m_ext_idx.tellg() < end) {
                disk_block db;
                db.load(m_ext_idx); // load block

                ++k_ir;             // increase number of irreducible blocks
                k_re += (db.header.size()-1); // increase the number of reducible blocks
                n_ir += db.sa.size();

                header_in_megabyte += util::get_size_in_mega_bytes(db.header);
                lcp_in_megabyte += util::get_size_in_mega_bytes(LcpSerializeWrapper(db.lcp));
                sa_in_megabyte += util::get_size_in_mega_bytes(db.sa);
// TODO: how to calculate fwd_id from bwd_id ??? 

                for (size_type i=0; i<db.header.size(); ++i) {
                    size_type bwd_id, delta_x, delta_d;
                    db.decode_header_triple(i, bwd_id, delta_x, delta_d);
					
                }
            }		
*/			
		}


        //! Constructor
        /*!
         *	\param file_name 		File name of the file which contains the text.
         *  \param b				Threshold parameter b for the maximum size of a block.
         *	\param output_tikz		Indicates if latex output should be produced for the internal and
         *							external part (so backward and forward part) of the index. The output
         *							is stored in the files file_name.fwd_idx.tex and file_name.bwd_idx.tex.
         *	\param delete_tmp		Indicates if the temporary files that are used during the construction should be removed after the construction.
         *  \param tmp_file_dir		Directory for the temporary files.
         *	\param output_dir		Directory for the output.
         *
         *	\Time complexity
         *		\f$ \Order{n \log\sigma} \f$, where n is the length of the text.
         */
        rosa(const char* file_name=NULL, size_type b=4096, bool output_tikz=false, bool delete_tmp=false,
             const char* tmp_file_dir="./", const char* output_dir=NULL):m_b(b), m_buf(NULL), m_buf_lz(NULL), m_buf_size(1024)
            ,bl(m_bl), bf(m_bf), wt(m_wt)
            ,bl_rank(m_bl_rank), bf_rank(m_bf_rank)
            ,bf_select(m_bf_select)
            ,bl_select(m_bl_select)
            ,bm(m_bm)
            ,bm_1_select(m_bm_1_select)
            ,bm_0_select(m_bm_0_select)
            ,bm_10_rank(m_bm_10_rank)
            ,min_depth(m_min_depth), pointer(m_pointer)
            ,file_name(m_file_name)
            ,output_dir(m_output_dir)
			,k(m_k)						 
			,lz_width(m_lz_width)
#ifdef OUTPUT_STATS
            ,count_disk_access(m_count_disk_access)
            ,count_gap_disk_access(m_count_gap_disk_access)
            ,count_int_steps(m_count_int_steps)
            ,count_int_match(m_count_int_match)
            ,count_queries(m_count_queries)
            ,count_block_length(m_count_block_length)
#endif
        {
            m_buf = new unsigned char[m_buf_size]; // initialise buffer for pattern
            m_buf_lz = new uint64_t[m_buf_size]; // initialise buffer for pattern
            if (NULL == file_name) {
                return; // if no file_name is specified do not construct the index
            }
//          (0) Initialise path for the external parts of the data structure
            m_file_name = string(file_name);
            m_output_dir = get_output_dir(file_name, output_dir);
            string path = m_output_dir + "/" + util::basename(file_name);
            string bwd_csa_file_name = path + "." + TMP_BWD_CSA_SUFFIX;
            string fwd_csa_file_name = path + "." + TMP_FWD_CSA_SUFFIX;
            string tmp_dir = (util::dirname(tmp_file_dir)+"/"+util::basename(tmp_file_dir)+"/");

            cache_config config(false, tmp_dir, util::basename(m_file_name));
			string lcp_file = util::cache_file_name(constants::KEY_LCP, config);
			string sa_file = util::cache_file_name(constants::KEY_SA, config);
			cout<<"sa_file = "<<sa_file<<endl;
			cout<<"lcp_file = "<<lcp_file<<endl;

//			(1) Load or construct the forward CST
            tCsa fwd_csa;
            write_R_output("fwd_csa","construct","begin");
            construct_or_load_fwd_csa(fwd_csa, fwd_csa_file_name, tmp_dir, delete_tmp);
            write_R_output("fwd_csa","construct","end");
			write_R_output("lcp","construct","begin");
			{
				ifstream lcp_in(lcp_file.c_str());
				if (!lcp_in){ // if the lcp file does not exist: construct it
					util::clear(fwd_csa);
					int_vector<> lcp;
					config.file_map[constants::KEY_SA] = sa_file;
					config.file_map[constants::KEY_TEXT] = util::cache_file_name(constants::KEY_TEXT, config);
					construct_lcp_kasai(lcp, config );
					util::store_to_file(lcp, lcp_file.c_str());
					util::clear(lcp);
            		construct_or_load_fwd_csa(fwd_csa, fwd_csa_file_name, tmp_dir, delete_tmp);
				}
			}
			write_R_output("lcp","construct","end");
			
/*
			std::string check_text = algorithm::extract(fwd_cst.csa, 0, fwd_cst.csa.size()-1);
			ofstream check_out((util::basename(file_name)+".check_text").c_str());
			check_out.write(check_text.c_str(), check_text.size());
			check_out.close();
*/

//          (2) Create the fwd_bf bit vector
            m_n = fwd_csa.size(); if (util::verbose) cout<<"m_n="<<m_n<<endl;
            bit_vector fwd_bf(m_n+1);
            fwd_bf[m_n] = 1; // mark last bit
            write_R_output("fwd_bf","construct","begin");
//            mark_blocks(fwd_cst, fwd_bf, b, m_k);

			mark_blocks(lcp_file.c_str(), fwd_bf, b, m_k); 
			cout<<"mark blocks: m_k="<<m_k<<endl;
            write_R_output("fwd_bf","construct","end");
            util::clear(fwd_csa);  if (util::verbose) cout<<"cleared fwd_csa"<<endl;
//          (3) Load or construct the backward CSA
            tCsa bwd_csa;
            write_R_output("bwd_csa","construct","begin");
            construct_or_load_bwd_csa(bwd_csa, bwd_csa_file_name, tmp_dir, delete_tmp);
            write_R_output("bwd_csa","construct","end");
            m_comp2char = bwd_csa.comp2char;
            m_char2comp = bwd_csa.char2comp;
//			(4) Create m_bl, m_bf and the mapping between fwd_ids and bwd_ids
            vector<block_info> map_info;
            write_R_output("bl,bf and bwd_id<->fwd_id mapping","construct","begin");
            calculate_bl_and_bf_and_mapping(bwd_csa, map_info);
            write_R_output("bl,bf and bwd_id<->fwd_id mapping","construct","end");
            util::clear(bwd_csa);  if (util::verbose) cout<<"cleared bwd_csa"<<endl;
            write_R_output("bwd_id<->fwd_id mapping","sort","begin");
            sort(map_info.begin(),map_info.end()); // sort according to bwd_lb, depth, fwd_lb and size
            write_R_output("bwd_id<->fwd_id mapping","sort","end");
//          (5) Create bm and min_depth
            write_R_output("bm and min_depth","construct","begin");
            calculate_bm_and_min_depth(map_info,  m_k);
            write_R_output("bm and min_depth","construct","end");
//          (6) Create the reducible graph
            write_R_output("fwd_csa","load","begin");
            construct_or_load_fwd_csa(fwd_csa, fwd_csa_file_name, tmp_dir, delete_tmp);
            write_R_output("fwd_csa","load","end");

            vector<block_node> v_block;  // block_node contains (delta_x, delta_d, dest_block, bwd_id)
            size_type red_blocks=0;
            size_type singleton_blocks = 0;
            size_type elements_in_irred_blocks = 0;
            write_R_output("reducible graph","construct","begin");
            calculate_reducible_graph(fwd_csa, fwd_bf, m_b, red_blocks, singleton_blocks, elements_in_irred_blocks, v_block);
            write_R_output("reducible graph","construct","end");
            // now the entries delta_x, delta_d and dest_block are known for each block in forward order.
            // still missing bwd_id
//          (7)
            util::assign(m_pointer, int_vector<>(m_k, 0, 64));
            write_R_output("bwd_id and singleton pointers","construct","begin");
			bit_vector is_singleton(m_k, 0);
            calculate_bwd_id_and_fill_singleton_pointers(sa_file.c_str(), fwd_csa, map_info, fwd_bf, v_block, is_singleton);
            util::clear(fwd_csa);   // cst not needed any more
            write_R_output("bwd_id and singleton pointers","construct","end");
//          (8) Calculate the headers of the external blocks
            vector<vector<header_item> > header_of_external_block(m_k);
            write_R_output("block_header","construct","begin");
            calculate_headers(v_block, map_info, header_of_external_block);
            write_R_output("block_header","construct","end");
//          (9) Write the external blocks
            write_R_output("external blocks","write","begin");

            write_external_blocks_and_pointers(sa_file.c_str(), lcp_file.c_str(), fwd_bf, v_block, header_of_external_block);
            write_R_output("external blocks","write","end");
//          (10) bit compress pointers
            util::bit_compress(m_pointer);
//			(11) init m_bl_rank_and_cC
            construct_or_load_bwd_csa(bwd_csa, bwd_csa_file_name, tmp_dir, delete_tmp);
            write_R_output("cC","construct","begin");
            init_bl_rank_and_cC(bwd_csa);
            write_R_output("cC","construct","end");
//          (12) construct condensed BWT (cBWT)
            write_R_output("cBWT","construct","begin");
            construct_condensed_bwt(bwd_csa, tmp_dir);
            write_R_output("cBWT","construct","end");
            util::clear(bwd_csa);
//			(13) greedy parse the text
			bit_vector factor_borders;
            write_R_output("parse","construct","begin");
			greedy_parse(tmp_dir, factor_borders);
            write_R_output("parse","construct","end");
//			(14) Replace SA text pointers by SA LZ-text pointers
            write_R_output("ext_idx","replace_pointers","begin");
			replace_pointers(factor_borders, is_singleton);
            write_R_output("ext_idx","replace_pointers","end");
			util::clear(factor_borders);
			util::clear(is_singleton);
//          (15) Open stream to text and external index for the matching
            open_streams();
        }

		//! Reconstruct text from the factorization and the condensed BWT
		void reconstruct_text(const std::string delimiter=""){
			string out_name = ("./"+util::basename(m_file_name)+".lz.txt");
			ofstream text_out(out_name.c_str());
			int_vector_file_buffer<> glz_buf(get_factorization_filename().c_str());
            for (size_type i=0,r=0,r_sum=0; i < glz_buf.int_vector_size;) { 
                for (; i < r+r_sum; ++i) {
					string factor_string;
					assert(glz_buf[i-r_sum]<m_k);
                   	extract_factor(glz_buf[i-r_sum],factor_string);
					if ( i+1 == glz_buf.int_vector_size ){
						factor_string = factor_string.substr(0, factor_string.size()-1);
					}
					text_out.write(factor_string.c_str(), factor_string.size());
					if ( delimiter.size() > 0 ){
						text_out.write(delimiter.c_str(), delimiter.size());
					}
                }
                r_sum += r; r = glz_buf.load_next_block();
            }
			text_out.close();	
			cout<<"The reconstructed text is stored in "<<out_name<<endl;
		}

		void factor_frequency(){
			string out_name = ("./"+util::basename(m_file_name)+".occ_freq");
			ofstream res_out(out_name.c_str());
			std::vector<size_type> freq(m_k, 0);
			size_type max_freq = 0;
			int_vector_file_buffer<> glz_buf(get_factorization_filename().c_str());
            for (size_type i=0,r=0,r_sum=0; i < glz_buf.int_vector_size;) { 
                for (; i < r+r_sum; ++i) {
                   	++freq[glz_buf[i-r_sum]];
					if ( freq[glz_buf[i-r_sum]] > max_freq ){
						max_freq = freq[glz_buf[i-r_sum]];
					}
                }
                r_sum += r; r = glz_buf.load_next_block();
            }
			std::vector<size_type> occ_freq(max_freq+1, 0);
		    for(size_type i=0; i<m_k; ++i){
				++occ_freq[freq[i]];
			}
			res_out<<"frequency factor_nr"<<endl;
			for(size_type i=0; i<occ_freq.size(); ++i){
				res_out << i << " " << occ_freq[i] << std::endl;
			}
			res_out.close();
		}

		/*!  
		 *  Adjust SA pointers in disk blocks to factorization
		 *  Adjust SA singleton pointers in condensed BWT
		 */
		void replace_pointers(const bit_vector& factor_borders, const bit_vector& is_singleton){
			rank_support_v<> factor_borders_rank(&factor_borders);
			ifstream tmp_ext_idx;
			open_stream(tmp_ext_idx, get_tmp_ext_idx_filename().c_str());
	        tmp_ext_idx.seekg(0, std::ios::end);
            std::streampos end = tmp_ext_idx.tellg(); // get the end position
            seekg(tmp_ext_idx, 0, false); // load the first block
			ofstream ext_idx(get_ext_idx_filename().c_str());
            // iterate through all blocks
            while (tmp_ext_idx.tellg() < end) {
                disk_block db;
                db.load(tmp_ext_idx); // load block
				db.replace_pointers(factor_borders_rank);
				db.serialize(ext_idx);
			}
			ext_idx.close();
			std::remove(get_tmp_ext_idx_filename().c_str());

			for (size_type i=0; i<m_k; ++i){
				if ( is_singleton[i] ) {
					m_pointer[i] = factor_borders_rank(m_pointer[i]);
				}
			}
		}

        //! Calculate the id of a block in the backward index using its left bound and depth.
        /*!
         *	\param lb		The left bound of the block in the backward index.
         *  \param depth	The depth of the block.
         *  \return	The backward id of the block (in [0..k-1], where k is the number of external blocks).
         *
         *	\par Time complexity
         *		 \f$ \Order{1} \f$
         */
        size_type get_bwd_id(size_type lb, size_type depth) const {
            /*            if (util::verbose) {
                            std::cout<<"get_bwd_id("<<lb<<","<<depth<<")\n";
                        }
            */            size_type run_nr = m_bf_rank(lb);
            size_type run_pos = 0;
            if (run_nr > 0) {
                run_pos = m_bm_1_select(run_nr)+1;
            }
            // optimization for the special case when the run is of size 1
            // m_mb[run_pos-1..run_pos] = 10
            if (m_bm[run_pos+1]) {  //i.e. m_mb[run_pos-1..run_pos+1] = 101
                // in the first (run_pos+1) position there are run_nr ones
                // => there are X=(run_pos+1)-run_nr zeros
                // => the id is X-1
                return run_pos - run_nr;
            } else { // i.e. m_mb[run_pos-1..run_pos+1] = 100...
                size_type min_depth = m_min_depth[m_bm_10_rank(run_pos+1)];
                return run_pos - run_nr + (depth-min_depth);
            }
        }

        //! List the occurrences of a pattern of length m
        /*!
         */
        size_type locate(const unsigned char* pattern, size_type m, vector<size_type>& res)const {
            res.resize(0);
            // (1) query in-memory data structure
            size_type lb, rb, d;
            if (get_interval(pattern, m, lb, rb, d)) {
                if (d == m) {  // the in-memory data structure answered the query
                    // we have to get all blocks between lb and rb and output the positions
                    if (util::verbose) {
                        std::cout<<"answered in-memory"<<std::endl;
                    }
                    return get_all_occurences(d, lb, rb, res);
                } else { // m > d
                    // (2) query external memory data structure
                    size_type bwd_id = get_bwd_id(lb, d);
                    if (rb+1-lb == 1) {
                        size_type sa = m_pointer[bwd_id];
                        if (match_pattern_lz(pattern, m, sa)) {
                            res.push_back(sa);
                            return 1;
                        }
                    } else {
                        size_type block_addr = m_pointer[bwd_id];
                        return search_block(pattern+d, m-d, d, rb+1-lb, bwd_id, block_addr, &res);
                    }
                }
            }
            return 0;

        }

        //! Count the number of occurrences of a pattern of length m
        /*!
         *  \param pattern 	A pointer to the start of the pattern.
         *  \param m		The length of the pattern.
         *  \param
         *  \return 		The number of occurrences of the pattern in the text.
         */
        size_type count(const unsigned char* pattern, size_type m, bool repeated_in_memory_search=false)const {
            if (util::verbose) {
                std::cout<<"count("<<pattern<<","<<m<<")\n";
            }
#ifdef OUTPUT_STATS
            ++m_count_queries;
#endif
            // (1) query in-memory data structure
            size_type lb, rb, d;
            if (get_interval(pattern, m, lb, rb, d)) {
                if (d == m) {  // the in-memory data structure answered the query
#ifdef OUTPUT_STATS
                    ++m_count_int_match;
#endif
                    return rb+1-lb;
                } else { // m > d
                    // (1a) try if the pattern P[d..m-1] and so on also appear in the text
                    if (repeated_in_memory_search) {
                        size_type td = d, tlb, trb;
                        const unsigned char* shorter_pattern = pattern + td;
                        while (m > td) {
                            size_type delta_d=0;
                            if (get_interval(shorter_pattern, m-td, tlb, trb, delta_d)) {
                                // delta_d should be > 0, since m-td > 0 and get_interval only returns true, if something is matched
                                td += delta_d;
                                shorter_pattern += delta_d;
//								if ( trb+1-tlb == 1 ){
//									// TODO: it is possible to check the text here with only one disk access
//								}
                            } else {
                                return 0;
                            }
                        }
                    }
                    // (2) query external memory data structure
                    size_type bwd_id = get_bwd_id(lb, d);
                    if (rb+1-lb == 1) {
                        size_type sa = m_pointer[bwd_id];
                        if (match_pattern_lz(pattern, m, sa))
                            return 1;
                    } else {
                        size_type block_addr = m_pointer[bwd_id];
                        return search_block(pattern+d, m-d, d, rb+1-lb, bwd_id, block_addr);
                    }
                }
            }
            return 0;
        }

        //! Match the pattern p reversed backwards against the pruned BWT until the interval <= b.
        /*! \param pattern  Pointer to the beginning of the pattern.
         *  \param m		The length of the pattern.
         *  \param lb		Left bound of the resulting interval in the backward index.
         *  \param rb		Right bound of the resulting interval in the backward index (inclusive).
         *	\param d		Number of matched characters.
         *	\return 		True if the pattern can occur in T and false otherwise.
         *
         *	\par Time complexity
         *		 \f$ \Order{m t_{rank\_bwt}} \f$
         */
        bool get_interval(const unsigned char* pattern, size_type m,
                          size_type& lb, size_type& rb, size_type& d)const {
            d  = 0;
            lb = 0; rb = m_n-1;
            if (util::verbose) {
                std::cout<<"pattern="<<pattern<<","<<m<<","<<lb<<","<<rb<<","<<d<<"\n";
            }
            const unsigned char* cp = pattern;//+m-1;
            while (d < m and rb-lb+1 > m_b) {
                if (util::verbose) {
                    std::cout<<d<<"-["<<lb<<","<<rb<<"]"<<" b="<<m_b<<std::endl;
                    std::cout<<"c="<< *cp <<" m_cC[m_char2comp["<<*cp<<"]]="<<m_cC[m_char2comp[*cp]]<<std::endl;
                }
                size_type lb1 = m_bl_rank(lb);
                size_type rb1 = m_bl_rank(rb+1);
                // Notice: rb1 > lb1, if bl is not empty
                unsigned char c = *cp++;
                size_type lb2 = m_wt.rank(lb1, c);
                size_type rb2 = m_wt.rank(rb1, c);
                if (lb2 == rb2) {// the pattern does not exist in the pruned BWT
                    return false;
                }
                ++d; // we have match another character
                lb = m_bf_select(m_cC[m_char2comp[c]] + lb2 + 1);
                rb = m_bf_select(m_cC[m_char2comp[c]] + rb2 + 1) - 1;
#ifdef OUTPUT_STATS
                ++m_count_int_steps;
#endif
            }
            return true;
        }

        void extract_factor(size_type bwd_id, string &factor_string)const {
			if(util::verbose){
				cout<<"extract_factor: bwd_id="<<bwd_id<<endl;
				cout<<"m_bm.size()="<<m_bm.size()<<" m_k="<<m_k<<endl;
			}
			assert ( bwd_id < m_k );
            size_type zero_pos	= m_bm_0_select(bwd_id+1);

            size_type ones 		= zero_pos - bwd_id; // ones in bm[0..zero_pos)
            size_type depth_idx = m_bm_10_rank(zero_pos+1);
			if(util::verbose){
				cout<<"extract_factor: depth_idx="<<depth_idx<<endl;
				cout<<"m_min_depth.size()="<<m_min_depth.size()<<endl;
			}
            size_type depth		= m_min_depth[depth_idx];
            if (ones) {
                depth += zero_pos - 1 - m_bm_1_select(ones);
            } else {
                depth += zero_pos;
            }
            size_type lb		= m_bf_select(ones + 1);
            unsigned char c = '\0';
			factor_string.resize(depth);
			if(util::verbose){
				cout<<"extract_factor: depth="<<depth<<endl;
			}
			size_type i=depth-1;
			do{
                c = first_row_character(lb);
				factor_string[i] = c;
				if ( i > 0 ){
					--i;
				}else
					break;
                size_type c_rank = m_bf_rank(lb)+m_bf[lb]-m_cC[m_char2comp[c]];
                size_type cpos   = m_wt.select(c_rank, c);
                lb = m_bl_select(cpos+1);
            }while(true);
        }
		
		bwd_block_info bwd_id_to_info(size_type bwd_id){
	            size_type zero_pos	= m_bm_0_select(bwd_id+1);
            size_type ones 		= zero_pos - bwd_id; // ones in bm[0..zero_pos)
            size_type depth_idx = m_bm_10_rank(zero_pos+1);
            size_type depth		= m_min_depth[depth_idx];
            if (ones) {
                depth += zero_pos - 1 - m_bm_1_select(ones);
            } else {
                depth += zero_pos;
            }
            size_type lb		= m_bf_select(ones + 1);
            unsigned char c = '\0';
            stack<unsigned char> factor;
            for (size_type i=0, _lb = lb; i < depth; ++i) {
                c = first_row_character(_lb);
                factor.push(c);
                size_type c_rank = m_bf_rank(_lb)+m_bf[_lb]-m_cC[m_char2comp[c]];
                size_type cpos   = m_wt.select(c_rank, c);
                _lb = m_bl_select(cpos+1);
            }
            size_type rb = m_n-1;
            while (!factor.empty()) {
                size_type rb1 = m_bl_rank(rb+1);
                c = factor.top(); factor.pop();
                size_type rb2 = m_wt.rank(rb1, c);
                rb = m_bf_select(m_cC[m_char2comp[c]] + rb2 + 1) - 1;
            }
			return bwd_block_info(lb, rb+1, depth, bwd_id);
		}

        unsigned char first_row_character(size_type i)const {
			// transform position in BWT into position in condensed BWT
            size_type ii = m_bf_rank(i) - (m_bf[i]==0);
            if (m_wt.sigma < 16) {
                size_type res = 1;
                while (m_cC[res] <= ii) {
                    ++res;
                }
                return m_comp2char[res-1];
            }else{
				size_type upper_c = m_wt.sigma, lower_c = 0; // lower_c inclusive, uppper_c exclusive
				size_type res = 0;
				do{
					res = (upper_c+lower_c)/2;
					if ( ii < m_cC[res] ){
						upper_c = res;
					}else{
						lower_c = res+1;
					}
				}while(ii < m_cC[res] or ii >= m_cC[res+1] );
				return m_comp2char[res];
			}
        }

        size_type size()const {
            return m_n;
        }

		void write_factor(uint64_t factor, ofstream &out, uint8_t num_bytes){
			if ( 8 == num_bytes ){
				out.write((char*)&factor, num_bytes);
			}else if( 4 == num_bytes){
				uint32_t x = factor;
				out.write((char*)&x, num_bytes);
			}
		}

		/*!
		 * Parse the text greedily into factors. Each factor (except the last one) corresponds to
		 * a block prefix of an external block. 
		 *
		 * \par Working set size
		 *      Constant size buffer for the text + constant size buffer for the output + n bits for factor borders 
		 */
        size_type greedy_parse(string tmp_dir, bit_vector &factor_borders) {
            if (m_b >= m_n) {
                throw std::logic_error("greedy_parse: m_b="+util::to_string(m_b)+" >= "+util::to_string(m_b)+"=m_n");
            }
            cache_config config(false, tmp_dir, util::basename(m_file_name));
            int_vector_file_buffer<8> text_buf(util::cache_file_name(constants::KEY_TEXT, config).c_str());
			// file name for temporary 
			string factor_file = tmp_dir+util::basename(m_file_name)+"_factors_"+util::to_string(util::get_pid())+"_"+util::to_string(util::get_id());
			uint8_t num_bytes = (bit_magic::l1BP(m_k-1)+1) > 32 ? 8 : 4;
//                              ^^^^^^^^^^^^^^^^^^^^^^^^^
//                              max #of bits for a factor 
			ofstream factor_stream(factor_file.c_str());
			if ( factor_stream ){
				util::assign(factor_borders, bit_vector(text_buf.int_vector_size,0));
				size_type factors = 0;
				size_type lb = 0, rb=m_n-1, d=0;

				for (size_type i=0,r=0,r_sum=0; i < text_buf.int_vector_size;) { 
					for (; i < r+r_sum; ++i) {
						size_type lb1 = m_bl_rank(lb);
						size_type rb1 = m_bl_rank(rb+1);
						uint8_t c = text_buf[i-r_sum];
						size_type lb2 = m_wt.rank(lb1, c);
						size_type rb2 = m_wt.rank(rb1, c);
						if (lb2 == rb2) {
							throw std::logic_error("greedy_parse: parse not possible with block prefixes");
						} else {
							++d; // we have matched another character
							lb = m_bf_select(m_cC[m_char2comp[c]] + lb2 + 1);
							rb = m_bf_select(m_cC[m_char2comp[c]] + rb2 + 1) - 1;
							if (rb-lb+1 <= m_b) { // reached 
								factor_borders[i] = 1;
								size_type bwd_id = get_bwd_id(lb, d);
								assert(bwd_id < m_k);
								write_factor(bwd_id, factor_stream, num_bytes);
								++factors;
								lb = 0; rb=m_n-1; d=0;
							}
						}
					}
					r_sum += r; r = text_buf.load_next_block();
				}
				factor_stream.close();
				int_vector<> factorization;
				util::load_vector_from_file(factorization, factor_file.c_str(), num_bytes);
				std::remove(factor_file.c_str()); // remove temp file
				util::bit_compress(factorization);
				m_lz_width = factorization.get_int_width();
            	util::store_to_file(factorization, get_factorization_filename().c_str());
            	return factors;
			}else{
				throw std::logic_error(("Greedy parse: Could not open temporary file: "+factor_file).c_str());	
				return 0;
			}
        }

        struct item {
            size_type d, lb, rb;
            item(size_type _d, size_type _lb, size_type _rb): d(_d), lb(_lb), rb(_rb) {};
        };
        struct block_item {
            size_type block_addr, bwd_id, size;
            block_item(size_type _block_addr, size_type _bwd_id, size_type _size):
                block_addr(_block_addr), bwd_id(_bwd_id), size(_size) {};
            bool operator<(const block_item& x)const {
                if (block_addr != x.block_addr) return block_addr < x.block_addr;
                if (bwd_id != x.bwd_id) return bwd_id < x.bwd_id;
                return size < x.size;
            }
        };


        //! Search for a pattern in the block
        /*!
         * \param d		 	 Number of characters matched so far.
         * \param lb		 Left bound of the lexicographic interval of the already matched prefix.
         * \param rb		 Right bound of the
         * \param bwd_id	 bwd_id of the matched prefix.
         * \param block_addr Start address of the block in the external index.
         * \param locations  Pointer to a vector, where the locations should be stored.
         * \param
         *
         * \par Time complexity
         *		2 disk accesses
         */
        size_type get_all_occurences(size_type d, size_type lb, size_type rb, std::vector<size_type>& loc)const {
            loc.resize(rb+1-lb);
            size_t loc_idx = 0;
            vector<block_item> res_block;
            queue<item> q;
            q.push(item(d, lb, rb));
            size_type block_size_sum=0;
            while (!q.empty()) {
                item x = q.front();
                q.pop();
                if (x.rb - x.lb + 1 > m_b) {
                    size_type k = 0;
                    std::vector<unsigned char> cs(m_wt.sigma);
                    std::vector<size_type> lb2(m_wt.sigma), rb2(m_wt.sigma);
                    size_type lb1 = m_bl_rank(x.lb);
                    size_type rb1 = m_bl_rank(x.rb+1);
                    m_wt.interval_symbols(lb1, rb1, k, cs, lb2, rb2);
                    for (size_type i=0; i<k; ++i) {
                        lb = m_bf_select(m_cC[m_char2comp[cs[i]]] + lb2[i] + 1);
                        rb = m_bf_select(m_cC[m_char2comp[cs[i]]] + rb2[i] + 1) - 1;
                        q.push(item(x.d+1, lb, rb));
                    }
                } else {
                    size_type bwd_id = get_bwd_id(x.lb, x.d);
                    if (x.rb+1-x.lb == 1) {
                        loc[loc_idx++] = m_pointer[bwd_id];
                    } else {
                        // push triple
                        res_block.push_back(block_item(m_pointer[bwd_id], bwd_id, x.rb+1-x.lb));
                        block_size_sum += x.rb+1-x.lb;
                    }
                }
            }
            if (util::verbose) {
                cout<<"queue done"<<endl;
            }
            sort(res_block.begin(), res_block.end()); // sort according to block_addr
            disk_block db;
            for (size_t i=0; i < res_block.size(); ++i) {
                if (0 == i or res_block[i-1].block_addr != res_block[i].block_addr) {
                    if (util::verbose) {
                        cout<<"i="<<i<<" res_block[i].block_addr="<<res_block[i].block_addr<<endl;
                    }
                    seekg(m_ext_idx, res_block[i].block_addr);
                    db.load(m_ext_idx);
                    if (util::verbose) {
                        cout<<"i="<<i<<" block_loaded loc.size()="<<loc.size()<<" loc_idx="<<loc_idx<<endl;
                    }
                }
                size_type delta_x = 0, delta_d = 0;
                db.get_delta_x_and_d(res_block[i].bwd_id, delta_x, delta_d);
                lb = delta_x;
                rb = delta_x + res_block[i].size - 1;
                for (size_type j=lb; j <= rb; ++j) {
                    loc[loc_idx++] = db.sa[j];
                }
            }
            if (util::verbose) {
                cout<<"ready"<<endl;
            }
            return loc.size();
        }


        //! This class represents a half-blind suffix tree for a disk block.
        /*!
         *  The half-blind suffix tree is constructed from a rosa disk-block in linear
         *  time. After that the get_interval operation efficiently determines
         *  if and how often a pattern can possible occur in the disk block.
         */
        class block_tree {
            public:
                typedef bit_vector::size_type  size_type;
                typedef bp_interval<size_type> node_type;
            private:
                const disk_block* m_db;             // pointer to the disk_block
                bp_support_sada<> m_bp_ct_support;  // balanced parentheses support for m_bp_ct

                // Get the first l index of a [i,j] interval.
                /* I.e. given an interval [i,j], the function returns the position of the smallest entry lcp[k] with \f$ i<k\leq j \f$
                 * \par Time complexity
                 * 	 \f$ \Order{1} \f$
                 */
                size_type get_first_l_index(const node_type& node, const bp_support_sada<>& bp_support, size_type& kpos, size_type& ckpos)const {
                    if (node.cipos > node.jp1pos) { // corresponds to m_lcp[i] <= m_lcp[j+1]
                        ckpos 	= node.jp1pos-1;
                    } else { // corresponds to m_lcp[i] > m_lcp[j+1]
                        ckpos	= node.cipos-1;
                    }
                    kpos	= bp_support.find_open(ckpos);
                    return bp_support.rank(kpos)-1;
                }

                // Get the next smaller value.
                /* \par Time complexity
                 *      \f$ \Order{1} \f$
                 */
                // Returns n if there is no next smaller value in [i+1..n-1]
                inline size_type nsv(size_type i, size_type ipos, const bp_support_sada<>& bp_support)const { // possible optimization: calculate also position of nsv, i.e. next ( following position cipos
                    size_type cipos = bp_support.find_close(ipos);
                    size_type result = bp_support.rank(cipos);
                    return result;
                }

                // Get the previous smaller value. Adapted for the case that there are no equal values in the array.
                /*
                 * \par Time complexity
                 *    \f$ \Order{\frac{\sigma}{w}} \f$, where w=64 is the word size, can be implemented in \f$\Order{1}\f$ with rank and select
                 */
                inline size_type psv(size_type i, size_type ipos, size_type cipos, size_type& psvpos, size_type& psvcpos, const bp_support_sada<>& bp_support)const {
                    if (i == 0) {  // if lcp[i]==0 => psv is the 0th index by definition
                        psvpos = 0;
                    } else {
                        psvpos = bp_support.enclose(ipos);
                    }
                    psvcpos = bp_support.find_close(psvpos);
                    return bp_support.rank(psvpos)-1;
                }

                //! Calculate the parent node of a node v.
                /*! \param v A valid node of the suffix tree.
                 *  \return The parent node of v or the root if v==root().
                 *  \par Time complexity
                 *       \f$ \Order{1}\f$
                 */
                node_type parent(const node_type& v, const bit_vector bp,const bp_support_sada<>& bp_support) const {
                    if (v.cipos > v.jp1pos) { // LCP[i] <= LCP[j+1]
                        size_type psv_pos, psv_cpos, psv_v, nsv_v, nsv_p1pos;
                        psv_v = psv(v.j+1, v.jp1pos, bp_support.find_close(v.jp1pos), psv_pos, psv_cpos, bp_support);
                        nsv_v = nsv(v.j+1, v.jp1pos, bp_support)-1;
                        if (nsv_v == size()-1) {
                            nsv_p1pos = bp.size();
                        } else { // nsv_v < size()-1
                            nsv_p1pos = bp_support.select(nsv_v+2);
                        }
                        return node_type(psv_v, nsv_v, psv_pos, psv_cpos, nsv_p1pos);
                    } else { // LCP[i] > LCP[j+1]
                        size_type psv_pos, psv_cpos, psv_v;
                        psv_v = psv(v.i, v.ipos, v.cipos, psv_pos, psv_cpos, bp_support);
                        return node_type(psv_v, v.j, psv_pos, psv_cpos, v.jp1pos);
                    }
                }

                //! Decide if a node is a leaf.
                bool is_leaf(const node_type& v)const {
                    return v.i==v.j;
                }

                bool is_root(const node_type& v)const {
                    return v.i==0 and v.j+1 == m_db->lcp.size();
                }
            public:
                block_tree(const disk_block& db):m_db(&db) {
                    util::init_support(m_bp_ct_support, &(m_db->bp_ct));
                    if (util::verbose) {
                        for (size_type i=0; i<m_db->lcp.size(); ++i) {
                            std::cout<<" "<<m_db->lcp[i];
                        } std::cout<<std::endl;
                        std::cout<<"m_bp_ct="<<m_db->bp_ct;
                    }
                }

                /*!
                 *  \param pattern  Pointer to the remaining pattern.
                 *  \param m		Length of the remaining pattern (in characters).
                 *	\param depth	Number of characters already matched (in characters).
                 *  \param lb		Left boundary of the search interval (inclusive).
                 *  \param rb   	Right boundary of the search interval (inclusive).
                 *	\return The number of potential occurrences of the pattern in the block.
                 *
                 *  Possible improvements: Determine in the search process, if there was
                 *                         a blind step. If not, we can report that and
                 *                         don't have to access disk for a check.
                 */
                size_type get_interval(const unsigned char* pattern, size_type m, size_type depth,
                                       size_type& lb, size_type& rb) {
                    // search process which can lead to a false positive result
                    size_type lb_pos = m_bp_ct_support.select(lb+1);
                    size_type rb_p1_pos = (rb+1 < size() ? m_bp_ct_support.select(rb+2) : 2*size());
                    node_type v(lb, rb, lb_pos, m_bp_ct_support.find_close(lb_pos), rb_p1_pos);
                    if (is_leaf(v)) {
                        return 1;
                    }

                    size_type v_bit_depth = 0;
                    node_type w = v;
                    while (!is_root(w)) {
                        node_type w_old = w;
                        if (util::verbose) {
                            std::cout<<w.i<<" "<<w.j<<" v_bit_depth="<<v_bit_depth<<std::endl;
                        }
                        w = parent(w, m_db->bp_ct, m_bp_ct_support);
                        if (w.j > w_old.j) {
                            v_bit_depth += m_db->lcp[w_old.j+1];
                        } else {
                            v_bit_depth += m_db->lcp[w_old.i];
                        }
                    }
                    v_bit_depth += m_db->lcp[0];
                    if (util::verbose) {
                        std::cout<<w.i<<" "<<w.j<<" v_bit_depth="<<v_bit_depth<<std::endl;
                    }

                    while (1) {
                        if (util::verbose) {
                            std::cout<<v.i<<" "<<v.j<<" v_bit_depth="<<v_bit_depth<<std::endl;
                        }
                        if (is_leaf(v)) {
                            return 1;
                        }
                        size_type lpos, clpos;
                        size_type l = get_first_l_index(v, m_bp_ct_support, lpos, clpos);
                        v_bit_depth += m_db->lcp[l];
                        if (util::verbose) {
                            std::cout<<"l-index="<<l<<" v_bit_depth="<<v_bit_depth<<std::endl;
                        }
                        if (v_bit_depth >= (depth+m)*8) {
                            return rb-lb+1;
                        }
                        unsigned char pc = *(pattern + (v_bit_depth/8 - depth));  // pattern char
                        if ((pc >> (7-(v_bit_depth%8)))&1) {   // 1-bit at depth v_bit_depth in the pattern
                            // => go to the right child in the tree
                            v = node_type(l ,v.j, lpos, clpos, v.jp1pos);
                        } else { // 0-bit at depth v_bit_depth in the pattern
                            // => go to the left child in the tree
                            v = node_type(v.i, l-1, v.ipos, v.cipos, lpos);
                        }
                        lb = v.i;
                        rb = v.j;
                    }
                    return 0;
                }

                size_type size()const {
                    return m_db->lcp.size();
                }
        };

        //! Search for a pattern in the block
        /*!
         * \param pattern    Pointer to the remaining pattern.
         * \param m			 Number of not yet matched characters. Equals the length of the remaining pattern.
         * \param depth 	 Number of characters matched so far.
         * \param size		 Size of the lexicographic interval of the already matched prefix.
         * \param bwd_id	 bwd_id of the matched prefix.
         * \param block_addr Start address of the block in the external index.
         * \param locations  Pointer to a vector, where the locations should be stored.
         *
         * \par Time complexity
         *		2 disk accesses
         */
        size_type search_block(const unsigned char* pattern, size_type m, size_type depth, size_type size,
                               size_type bwd_id, size_type block_addr, std::vector<size_type>* loc=NULL)const {
            if (util::verbose) {
                std::cout<<"serach_block("<<pattern<<","<<m<<","<<depth<<","<<size<<","<<bwd_id<<","<<block_addr<<")\n";
            }
            disk_block db;				                     // disk block object
            seekg(m_ext_idx, block_addr);                    // seek to the start address of the disk block
            db.load(m_ext_idx);			                     // fetch the disk block (load the header, LCP
#ifdef OUTPUT_STATS
            m_count_block_length += db.size(); // not db.size() >= size
#endif
            size_type delta_x = 0, delta_d = 0;
            db.get_delta_x_and_d(bwd_id, delta_x, delta_d);
            size_type lb = delta_x, rb = delta_x+size-1;     // determine left and right bound
#ifdef BENCHMARK_LOAD_ONLY
            // don't construct and match
#else
            block_tree tree(db);							 // create tree structure out of the block
#ifdef BENCHMARK_CREATE_ONLY
            // don't to the actual matching but something; but use the block_tree, so that
            // the construction is no optimized out
            size = std::min(size, tree.size());
#else
            size = tree.get_interval(pattern, m, depth+delta_d, lb, rb); // do the matching
#endif // END BENCHMARK_CREATE_ONLY
#endif // END BENCHMARK_LOAD_ONLY
            // if the search interval is not of size 0, and the check for the pattern
            // on disk is successful
#if defined BENCHMARK_LOAD_ONLY || defined BENCHMARK_CREATE_ONLY || defined BENCHMARK_SEARCH_BLOCK_ONLY
            return size;
#endif
            if (size > 0 and match_pattern_lz(pattern-depth, m+depth, db.sa[lb])) {
                if (NULL != loc) {
                    for (size_type i=lb; i<=rb; ++i) {
                        loc->push_back(db.sa[i]);
                    }
                }
                return size; // return the interval size
            } else {
                return 0;
            }

        }

        //! Check if pattern is a prefix of the suffix starting at text_offset in T.
        /*! \param pattern		A pointer to the start of the pattern.
         *  \param m			Length of the pattern.
         *  \param text_offset	Starting position in the text.
         *  \return				True if is a real pre
         *	\par Time complexity
         *		\f$  \Order{m} \f$ steps and one disk access
         */
         bool match_pattern(const unsigned char* pattern, size_type m, size_type text_offset)const {
            if (0==m or text_offset >= m_n-1) {  // avoid disk access if m==0
                return false;
            }
            // should be: text_offset+m <= m_n-1
            if (text_offset + m > m_n-1) {   // m_n-1 text size without the sentinel character
				return false;
            }
            seekg(m_text, text_offset);
            while (m > 0) {
                size_type len = std::min(m, m_buf_size);
                m_text.read((char*)m_buf, len);
                int res;
                if (0 != (res=memcmp(pattern, m_buf, len))) {
					return false;
                }
                m -= len;
                pattern += len;
            }
            return true;
        }
		
		 void calculate_kmp_table(const unsigned char* pattern, size_type m, int_vector<> &kmp_table)const {
			if(util::verbose) cout<<"KMP_TABLE="<<kmp_table[0];
			size_type i=1, j=0; 
			while ( i < m ){
				if ( pattern[i] == pattern[j] ){
					kmp_table[i++] = ++j;
if(util::verbose) cout<<" "<<kmp_table[i-1];
				}else{
					if(j>0){
						j = kmp_table[j-1];	
					}else{
						kmp_table[i] = 0;
if(util::verbose) cout<<" "<<kmp_table[i];
						++i;
					}
				}
			}
			if(util::verbose)cout<<endl;
		 } 

		//! Check if pattern is a prefix of a suffix starting in the factor lz_offset
		/* \param pattern		A pointer to the start of the pattern.
		 * \param m				Length of the pattern.
		 * \param lz_offset		Index to the factor 
		 */
		bool match_pattern_lz(const unsigned char* pattern, size_type m, size_type lz_offset)const{
			if ( util::verbose ) cout<<"match_pattern_lz("<<pattern<<","<<m<<","<<lz_offset<<")"<<endl;
			if ( m == 0 or lz_offset >= m_lz_size-1 ){ // if lz_offset == m_lz_size-1 then it is the terminal character
				return false;
			}
			size_type lz_bit_offset = m_lz_width*lz_offset;
			if ( util::verbose ) cout<<"lz_offset = "<<lz_offset<<"  /  lz_size = "<<m_lz_size<<endl;
			if ( util::verbose ) cout<<"lz_bit_offset = "<<lz_bit_offset<<endl;
			
			seekg(m_glz_text, ((lz_bit_offset)/64)*8+9);
			
			int_vector<> kmp_table(m, 0, bit_magic::l1BP(m)+1);
			calculate_kmp_table(pattern, m, kmp_table);
//TODO: handle the case where we load multiple blocks form disk				
			size_type len = ((m_lz_width*(lz_offset+m)-1)/64)-((m_lz_width*lz_offset)/64)+1; // len in 64 bit words
			len = std::min(len, m_buf_size);
			
			size_type cur_end_word = (m_lz_width*lz_offset)/64+len;
			size_type the_end_word = (m_lz_width*m_lz_size+63)/64;
			if ( cur_end_word > the_end_word ){ // check that we don't go beyond EOF
			    len = the_end_word - ((m_lz_width*lz_offset)/64);
			}

			if ( util::verbose ) cout<<"len="<<len<<endl;
			size_type matched = 0;
			
			/// TODO take care about failbit!!!

//			cout<<"m_glz_text.fail()="<<m_glz_text.fail()<<endl;
			m_glz_text.read((char*)m_buf_lz, 8*len);
				
			uint8_t bit_offset = lz_bit_offset&0x3F;
			if ( util::verbose ) cout<<"bit_offset="<<(int)bit_offset<<endl;
			const uint64_t *word = m_buf_lz;
			while ( word < m_buf_lz+len-1 or (word == m_buf_lz+len-1 and bit_offset+m_lz_width <= 64) ){
				uint64_t factor_id = bit_magic::read_int_and_move(word, bit_offset, m_lz_width);
				++lz_offset;
				if ( util::verbose ) cout << "factor id=" << factor_id << 
					" word="<<word <<
					" bit_offset="<<(int)bit_offset<< " m_lz_width="<<(int)m_lz_width << endl;
				string factor_string;
				extract_factor(factor_id, factor_string); 
				if ( util::verbose ) cout << "factor string=" << factor_string << endl;
				for(size_type i=0; i<factor_string.size(); ){
					if ( (unsigned char)factor_string[i] == pattern[matched] ){
						if (matched==m-1)
							return true;
						++i; ++matched;
					}else{
//						if( util::verbose  and factor_string.size()==67){
//							cout<<"factor_string["<<i<<"]="<<factor_string[i]<<"!="<<pattern[matched]<<"=pattern["<<matched<<"]"<<endl;
//							cout<<(int)factor_string[i]<<"!="<<(int)pattern[matched]<<endl;
//						}
						if ( matched > 0 ){
							matched = kmp_table[matched-1];
						}else{
							++i;
						}
					}
				}
				if( util::verbose ) cout << "matched so far = "<<matched<<endl;
				if (factor_id==0)
					break;
			}

			return false;
		}

		



        void reset_counters() {
#ifdef OUTPUT_STATS
            m_count_disk_access 	= 0;
            m_count_gap_disk_access = 0;
            m_count_int_steps   	= 0;
            m_count_int_match   	= 0;
            m_count_queries			= 0;
            m_count_block_length    = 0;
#endif
        }

        //! Get the size of the external part in megabytes.
        double get_ext_idx_size_in_mega_byte() {
            return ((double)util::get_file_size(get_ext_idx_filename().c_str()))/ (1024.0*1024.0);
        }

        //! Writes the in-memory part of the index into the output stream.
        /*! \param out Output stream to which the index should be written.
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += util::write_member(m_n, out, child, "n");
            written_bytes += util::write_member(m_b, out, child, "b");
            written_bytes += util::write_member(m_k, out, child, "k");
            written_bytes += m_bl.serialize(out, child, "bl");
            written_bytes += m_bf.serialize(out, child, "bf");
            written_bytes += m_bl_rank.serialize(out, child, "bl_rank");
            written_bytes += m_bf_rank.serialize(out, child, "bl_select");
            written_bytes += m_bf_select.serialize(out, child, "bf_select");
            written_bytes += m_bl_select.serialize(out, child, "bl_select");
            written_bytes += m_wt.serialize(out, child, "wt");
            written_bytes += m_cC.serialize(out, child, "cC");
            written_bytes += m_char2comp.serialize(out, child, "char2comp");
            written_bytes += m_comp2char.serialize(out, child, "comp2char");
            written_bytes += m_bm.serialize(out, child, "bm");
            written_bytes += m_bm_1_select.serialize(out, child, "bm_1_select");
            written_bytes += m_bm_0_select.serialize(out, child, "bm_0_select");
            written_bytes += m_bm_10_rank.serialize(out, child, "bp_01_rank");
            written_bytes += m_min_depth.serialize(out, child, "min_depth");
            written_bytes += m_pointer.serialize(out, child, "pointer");
            written_bytes += util::write_member(m_file_name, out, child, "file_name");
            written_bytes += util::write_member(m_output_dir, out, child, "output_dir");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the in-memory part of the index from the input stream.
        /*! \param in Input stream from which the index should be read.
         */
        void load(std::istream& in) {
            util::read_member(m_n, in);
            util::read_member(m_b, in);
            util::read_member(m_k, in);
            m_bl.load(in);
            m_bf.load(in);
            m_bl_rank.load(in, &m_bl);
            m_bf_rank.load(in, &m_bf);
            m_bf_select.load(in, &m_bf);
            m_bl_select.load(in, &m_bl);
            m_wt.load(in);
            m_cC.load(in);
            m_char2comp.load(in);
            m_comp2char.load(in);
            m_bm.load(in);
            m_bm_1_select.load(in, &m_bm);
            m_bm_0_select.load(in, &m_bm);
            m_bm_10_rank.load(in, &m_bm);
            m_min_depth.load(in);
            m_pointer.load(in);
            util::read_member(m_file_name, in);
            util::read_member(m_output_dir, in);
            open_streams();
        }

        void set_file_name(const string& file_name) {
            m_file_name = file_name;
            open_streams();
        }

        void set_output_dir(const string& output_dir) {
            m_output_dir = output_dir;
            open_streams();
        }

        //! Output statistics about the data structure
        // Maybe TODO: output the number of nodes of the corresponding trie
        void statistics()const {
            size_type k_ir = 0; // counter for irreducible blocks
            size_type k_re = 0; // counter for reducible blocks
            size_type k_s  = 0; // counter for singleton blocks

            size_type n_ir = 0; // number of elements in irreducible blocks

            double header_in_megabyte = 0.0;
            double lcp_in_megabyte    = 0.0;
            double sa_in_megabyte     = 0.0;

            size_type max_bwd_id = 0;
            size_type max_delta_x = 0;
            size_type max_delta_d = 0;

            m_ext_idx.seekg(0, std::ios::end);
            std::streampos end = m_ext_idx.tellg(); // get the end position
            seekg(m_ext_idx, 0, false); // load the first block
            // iterate through all blocks
            while (m_ext_idx.tellg() < end) {
                disk_block db;
                db.load(m_ext_idx); // load block
#ifdef OUTPUT_BLOCK_BOUNDS
                std::cout << m_ext_idx.tellg() << std::endl;
#endif
                ++k_ir;             // increase number of irreducible blocks
                k_re += (db.header.size()-1); // increase the number of reducible blocks
                n_ir += db.sa.size();

                header_in_megabyte += util::get_size_in_mega_bytes(db.header);
                lcp_in_megabyte += util::get_size_in_mega_bytes(LcpSerializeWrapper(db.lcp));
                sa_in_megabyte += util::get_size_in_mega_bytes(db.sa);


                for (size_type i=0; i<db.header.size(); ++i) {
                    size_type bwd_id, delta_x, delta_d;
                    db.decode_header_triple(i, bwd_id, delta_x, delta_d);

                    max_bwd_id = std::max(max_bwd_id, bwd_id);
                    max_delta_x = std::max(max_delta_x, delta_x);
                    max_delta_d = std::max(max_delta_d, delta_d);
                }
            }
            k_s = m_k-k_ir-k_re; // calculate the number of singleton blocks

            std::cout << "# file_name = " <<m_file_name << std::endl;
            std::cout << "# index_suf = " << INDEX_SUF << std::endl;
            std::cout << "# sigma = " << (int)m_wt.sigma << std::endl;
            std::cout << "# n = " << m_n << std::endl;
            std::cout << "# b = " << m_b << std::endl;
            std::cout << "# k = " << m_k << std::endl;
            std::cout << "# k_ir = " << k_ir << std::endl;
            std::cout << "# k_re = " << k_re << std::endl;
            std::cout << "# k_s = " << k_s << std::endl;
            std::cout << "# n_ir = " << n_ir << std::endl;
            std::cout << "# k_bf = " << m_bf_rank(m_bf.size()) << std::endl;
            std::cout << "# k_bl = " << m_bl_rank(m_bl.size()) << std::endl;
            std::cout << "# bm = " << m_bm.size() << std::endl;
            std::cout << "# max_bwd_id = " << max_bwd_id << std::endl;
            std::cout << "# max_delta_x = " << max_delta_x << std::endl;
            std::cout << "# max_delta_d = " << max_delta_d << std::endl;
            std::cout << "# header_in_megabyte = " << header_in_megabyte << std::endl;
            std::cout << "# lcp_in_megabyte = " << lcp_in_megabyte << std::endl;
            std::cout << "# sa_in_megabyte = " << sa_in_megabyte << std::endl;
            std::cout << "# label_in_megabyte = " << 0 << std::endl;
        }

};

#endif // end of include guard
