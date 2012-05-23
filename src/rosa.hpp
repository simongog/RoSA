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

using namespace sdsl;
using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::vector;
using std::ifstream;
using std::map;
using std::ostream;

// forward declaration of the rosa class with default template parameters
template<class BitVector = bit_vector // for bl and bf
,class RankSupport = typename BitVector::rank_1_type // for bl and bf
,class SelectSupport = typename BitVector::select_1_type // for bf
,class WaveletTree = wt_huff<bit_vector,
rank_support_v5<>,
select_support_dummy,
select_support_dummy> // for pruned BWT
,class LcpSerializeWrapper = int_vector_serialize_vbyte_wrapper<> // int_vector_serialize_wrapper<>
,class LcpLoadWrapper = int_vector_load_vbyte_wrapper<>  // int_vector_load_wrapper<>
>
class rosa;

#define TMP_CST_SUFFIX "tmp_fwd_cst4"
#define TMP_CSA_SUFFIX "tmp_bwd_csa4"


typedef bit_vector::size_type size_type;
typedef std::pair<size_type, size_type> tPII;
typedef vector< tPII > tVPII;
typedef std::queue< tPII > tQPII;
typedef std::queue< size_type > tQI;
typedef std::pair<unsigned char, size_type> tPCI;
typedef std::priority_queue<tPCI, vector<tPCI>, std::greater<tPCI> > tPQPCI;

#ifdef OUTPUT_STATS
size_type gl_X;
#endif


template<class BitVector
,class RankSupport
,class SelectSupport
,class WaveletTree
,class LcpSerializeWrapper
,class LcpLoadWrapper
>
class rosa
{
    public:
        typedef int_vector<>::size_type size_type;
        typedef BitVector bit_vector_type;
        typedef RankSupport rank_support_type;
        typedef SelectSupport select_support_type;
        typedef WaveletTree wavelet_tree_type;
        typedef csa_wt<wt_huff<>,64000, 64000> tCsa;
        typedef cst_sada<csa_wt<wt_huff<>,4, 64000>, lcp_dac<> > tCst;
        typedef tCst::node_type node_type;
        typedef select_support_mcl<> bm_select_type;
        typedef rank_support_v5<10,2> bm_rank01_type;
    private:
        size_type				m_n;  // original text length
        size_type 				m_b;  // block size
        size_type				m_k;  // number of external blocks
        bit_vector_type 		m_bl; // indicates the BWT entries which corresponds to a representative suffix
        bit_vector_type 		m_bf; // indicates the representative suffixes in SA
        rank_support_type 		m_bl_rank; // rank support for m_bl
        rank_support_type 		m_bf_rank; // rank support for m_bf
        select_support_type 	m_bf_select; // select support for m_bf
        wavelet_tree_type		m_wt; // wavelet tree for pruned BWT
        int_vector<64>			m_C1; // contains for each character c the m_rank_bf(C[c])
        bit_vector				m_bm; // bit sequence to map from backward intervals to blocks in the external part
        bm_select_type			m_bm_select; //
        bm_rank01_type			m_bm_rank01; //
        int_vector<>			m_min_depth;
        int_vector<>			m_pointer; // address of block [sp..ep] on disk,
        // if ep-sp > 0 or SA[sp] if sp=ep.

        string					m_file_name;  // file name of the supported text
        string					m_output_dir; // output directory

        mutable ifstream	m_text;	      // stream to the text, is needed to check
        // the
        mutable ifstream	m_ext_idx;    // stream to the external memory part
        const size_type  		m_buf_size;
        unsigned char*			m_buf;

#ifdef OUTPUT_STATS
        mutable size_type m_count_disk_access; 		// see corresponding public member
        mutable size_type m_count_gap_disk_access; 	// see corresponding public member
        mutable size_type m_count_int_steps;   		// see corresponding public member
        mutable size_type m_count_int_match;   		// see corresponding public member
        mutable size_type m_count_queries;	   		// see corresponding public member
        mutable size_type m_X;
#endif

        //! Internal helper class for an external block
        /*  TODO: possible optimizations:
         *        * don't store the triple of the corresponding irreducible block
         *          and use delta_x=delta_d=0 if we don't find its backward id
         *        * subtract and add d from the lcp values
         */
        class disk_block
        {
                int_vector<> m_header;	// Each entry contains an encoded header triple.
                uint8_t m_width_bwd_id;	// Bit width of the maximal value of a backward id in the header.
                uint8_t m_width_delta_x;// Bit width of the maximal delta x in the header.
                int_vector<>  m_lcp; 	// Array for the lcp values.
                int_vector<>  m_sa;		// Array for the suffix array values.

            public:

                const int_vector<>& header;	// const reference to the header
                const int_vector<>& lcp;	// const reference to the lcp values
                const int_vector<>& sa;		// const reference to the suffix array values

                disk_block():header(m_header), lcp(m_lcp), sa(m_sa) {};


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
                 * \pre rb>lb
                 */
                template<class tCst>
                void set_content(const tCst& cst, const unsigned char* text, size_type lb, size_type rb) {
                    size_type max_sa = 0;
                    size_type max_lcp = 0;
                    int_vector<64> t_csa(rb-lb+1); // buffer for csa values
                    int_vector<64> t_lcp(rb-lb+1); // buffer for lcp values
//				size_type d = get_block_depth(cst, lb, rb);
                    // determine maximum of the elements
                    for (size_type i=lb; i <= rb ; ++i) {
                        t_csa[i-lb] = cst.csa[i];
                        t_lcp[i-lb] = cst.lcp[i];
                        max_sa = std::max(max_sa, t_csa[i-lb]);
                        max_lcp = std::max(max_lcp, t_lcp[i-lb]);
                    }
                    m_sa = int_vector<>(rb-lb+1, 0, bit_magic::l1BP(max_sa)+1);
                    m_lcp = int_vector<>(rb-lb+1, 0, bit_magic::l1BP(max_lcp/*-d*/)+1);
//				std::cout << "set_context m_sa.size()=" << m_sa.size() << endl;
                    // copy elements
                    for (size_type i=lb; i <= rb; ++i) {
                        m_sa[i-lb] 		= t_csa[i-lb];
                        m_lcp[i-lb] 	= t_lcp[i-lb]/*-d*/;
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
                size_type serialize(std::ostream& out)const {
                    size_type written_bytes = 0;
                    written_bytes += util::write_member(m_width_bwd_id, out);
                    written_bytes += util::write_member(m_width_delta_x, out);
                    written_bytes += m_header.serialize(out);

                    written_bytes += LcpSerializeWrapper(m_lcp).serialize(out);
                    written_bytes += m_sa.serialize(out);
                    return written_bytes;
                }

                //! Load the block from the input stream
                void load(std::istream& in) {
#ifdef DEBUG_DISK_ACCESS
                    cout << "load block" << endl;
#endif
                    util::read_member(m_width_bwd_id, in);
                    util::read_member(m_width_delta_x, in);
#ifdef DEBUG_DISK_ACCESS
                    cout << "loaded first two header items" << endl;
#endif
                    m_header.load(in);
#ifdef DEBUG_DISK_ACCESS
                    cout << "m_header.size() = "<<m_header.size() << endl;
#endif
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
            close_stream_if_open(m_text);
            m_text.open(m_file_name.c_str());
            if (!m_text) {
                std::cerr << "Error: Could not open file: " << m_file_name << std::endl;
            } else {
                std::cerr << "Opened file " << m_file_name << std::endl;
            }
            close_stream_if_open(m_ext_idx);
            m_ext_idx.open(get_ext_idx_filename().c_str());
            if (!m_ext_idx) {
                std::cerr << "Error: Could not open file: " << get_ext_idx_filename() << std::endl;
            } else {
                std::cerr << "Opened file "<< get_ext_idx_filename() << std::endl;
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
        const bit_vector& bm;					//!< Bit vector to map from intervals in the backward index to blocks in the external part.
        const select_support_mcl<>& bm_select;  //!< Select support for bm.
        const bm_rank01_type& bm_rank01;		//!< Rank support for the bit-pattern 01 in bm.
        const int_vector<>& min_depth;			//!<
        const int_vector<>& pointer;			//!< Array of pointers into the external structures (blocks or suffix array).
        const string& file_name;				//!< File name of the original text string.
        const string& output_dir;				//!< Directory where the output is stored.
        const string tmp_cst_suffix;			//!< File name of the temporary generated suffix tree.
#ifdef OUTPUT_STATS
        const size_type& count_disk_access; 	//!< Counter for potential disk accesses.
        const size_type& count_gap_disk_access; //!< Counter for potential disk accesses caused by gaps in the fringe.
        const size_type& count_int_steps; 		//!< Counter for matches that can be answered with the in-memory part of the data structure.
        const size_type& count_int_match; 		//!< Counter for matches that can be answered with the in-memory part of the data structure.
        const size_type& count_queries;			//!< Counter for the queries.
#endif

        ~rosa() {
            m_text.close();    // close stream to the text
            m_ext_idx.close(); // close stream to the external part
            delete [] m_buf;
        }

        //! Get the name of the file where the external part of the index is stored.
        static string get_ext_idx_filename(const char* file_name, size_type b,  const char* output_dir=NULL) {
            return get_output_dir(file_name, output_dir) + "/" + util::basename(file_name)
                   + "." + util::to_string(b) + ".ext_idx";
        }

        //! Get the name of the file where the external part of the index is stored.
        string get_ext_idx_filename() {
            return get_ext_idx_filename(m_file_name.c_str(), m_b, m_output_dir.c_str());
        }


        //! Get the name of the file where the in-memory part of the index is stored.
        static string get_int_idx_filename(const char* file_name, size_type b, const char* output_dir=NULL) {

            return get_output_dir(file_name, output_dir) + "/" + util::basename(file_name)
                   + "." + util::to_string(b) + "."+INDEX_SUF+".int_idx";
        }

        //! Get the name of the file where the in-memory part of the index is stored.
        string get_int_idx_filename() {
            return get_int_idx_filename(m_file_name.c_str(), m_b, m_output_dir.c_str());
        }

        //! Remove the temporary files which are created during the construction.
        static void remove_tmp_files(const char* file_name, const char* output_dir) {
            std::string path = get_output_dir(file_name, output_dir) + "/" + util::basename(file_name);
            std::remove((path + "." + TMP_CST_SUFFIX).c_str());
            std::remove((path + "." + TMP_CSA_SUFFIX).c_str());
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
        rosa(const char* file_name=NULL, size_type b=4000, bool output_tikz=false, bool delete_tmp=false,
             const char* tmp_file_dir="./", const char* output_dir=NULL):m_b(b), m_buf_size(1024)
            ,bl(m_bl), bf(m_bf), wt(m_wt)
            ,bl_rank(m_bl_rank), bf_rank(m_bf_rank)
            ,bf_select(m_bf_select), bm(m_bm)
            ,bm_select(m_bm_select), bm_rank01(m_bm_rank01)
            ,min_depth(m_min_depth), pointer(m_pointer)
            ,file_name(m_file_name)
            ,output_dir(m_output_dir)
#ifdef OUTPUT_STATS
            ,count_disk_access(m_count_disk_access)
            ,count_gap_disk_access(m_count_gap_disk_access)
            ,count_int_steps(m_count_int_steps)
            ,count_int_match(m_count_int_match)
            ,count_queries(m_count_queries)
#endif
        {
            m_buf = new unsigned char[m_buf_size];
            if (NULL == file_name) {  // if no file_name is specified do not construct the index
                return;
            }
            tCsa bwd_csa;
            tCst fwd_cst;
            m_file_name = string(file_name);
            m_output_dir = get_output_dir(file_name, output_dir);
            std::cout<<"m_output_dir="<<m_output_dir<<std::endl;
            std::string path = m_output_dir + "/" + util::basename(file_name) ;
            string tmp_bwd_csa_file_name = path + "." + TMP_CSA_SUFFIX;
            string tmp_fwd_cst_file_name = path + "." + TMP_CST_SUFFIX;

            std::string tmp_file_dir2 = (util::dirname(tmp_file_dir)+"/"+util::basename(tmp_file_dir)+"/");
            std::cout<<"tmp_file_dir2="<<tmp_file_dir2<<std::endl;

            ifstream tmp_bwd_csa_stream(tmp_bwd_csa_file_name.c_str());
            if (!tmp_bwd_csa_stream) {
                tMSS file_map;
                construct_csa_of_reversed_text(m_file_name, bwd_csa, file_map, delete_tmp, tmp_file_dir2, ("reversed_"+util::basename(m_file_name)).c_str());
                util::store_to_file(bwd_csa, tmp_bwd_csa_file_name.c_str());
            } else {
                std::cout<< "load stored bwd_csa form file "<< tmp_bwd_csa_file_name << std::endl;
                util::load_from_file(bwd_csa, tmp_bwd_csa_file_name.c_str());
            }
            std::cout<<"bwd_csa.sigma="<<(int)bwd_csa.sigma<<std::endl;
            ifstream tmp_fwd_cst_stream(tmp_fwd_cst_file_name.c_str());
            if (!tmp_fwd_cst_stream) {
                tMSS file_map;
                std::cout<<"m_file_name="<<m_file_name<<std::endl;
                std::cout<<"tmp_file_dir2="<<tmp_file_dir2<<std::endl;
                std::cout<<"id="<<util::basename(m_file_name)<<std::endl;
                construct_cst(m_file_name, fwd_cst, file_map, delete_tmp, tmp_file_dir2, false, util::basename(m_file_name));
                util::store_to_file(fwd_cst, tmp_fwd_cst_file_name.c_str());
            } else {
                tmp_fwd_cst_stream.close();
                std::cout<< "load stored fwd_cst form file "<< tmp_fwd_cst_file_name << std::endl;
                util::load_from_file(fwd_cst, tmp_fwd_cst_file_name.c_str());
            }

            m_n = bwd_csa.size();

            bit_vector fwd_bf(m_n+1);
            fwd_bf[m_n] = 1;
            size_type trie_nodes = 1;
            size_type red_blocks=0, red_blocks_delta0 = 0;
            size_type singleton_blocks = 0;
            // guarantee that only one additional disk access is required
            size_type elements_in_irred_blocks = 0;
            get_block_info_and_mark_blocks(fwd_cst, fwd_bf, b,trie_nodes, m_k);
            rank_support_v5<> fwd_bf_rank(&fwd_bf);

#ifdef OUTPUT_STATS
            cout << "# text_size = " << fwd_cst.csa.size() << endl;
            cout << "# k = "<< m_k << endl;
            cout << "# trie_nodes = "<< trie_nodes << endl;
#endif
            vector<block_node> v_block;
            get_block_info(fwd_cst, fwd_bf, b,red_blocks, red_blocks_delta0,
                           singleton_blocks, elements_in_irred_blocks, v_block);

#ifdef OUTPUT_STATS
            cout << "# red_blocks = " << red_blocks << endl;
            cout << "ratio_red: " << ((double)red_blocks)/m_k<<endl;
            cout << "# red_blocks_delta0 = " << red_blocks_delta0 << endl;
            cout << "ratio_red_delta0: " << ((double)red_blocks_delta0)/m_k<<endl;
            cout << "# singleton_blocks = " << singleton_blocks << endl;
            cout << "ratio_singleton_blocks: " << ((double)singleton_blocks)/m_k<<endl;
            cout << "# elements_in_irred_blocks = " << elements_in_irred_blocks << endl;
            cout << "ratio_irred_elements: " << elements_in_irred_blocks/(double)fwd_cst.csa.size() << endl;
#endif
            output_block_info_statistics(v_block);

            vector<bwd_block_info> fwd_blocks_in_bwd;
            vector<vector<header_item> > external_block(m_k);
            {
                bit_vector bl = bit_vector(m_n);
                bit_vector bf = bit_vector(m_n +1);   // add one bit for technical reasons
                tQPII q;
                q.push(tPII(0, bwd_csa.size()));   // init queue with root interval
                tQI q_fp; // queue for the forward position
                q_fp.push(0);

                tQI q_depth; q_depth.push(0); // for statistics

                size_type k;
                vector<unsigned char> cs(256);
                vector<size_type> rank_c_lb(256), rank_c_rb(256);
                std::vector<block_info> map_info;

                tPQPCI pq; // min priority queue for calculating the forward position
                while (!q.empty()) {
                    size_type lb = q.front().first;
                    size_type rb = q.front().second;
                    size_type lb_forward = q_fp.front();
                    q.pop(); q_fp.pop();
                    size_type depth = q_depth.front(); q_depth.pop();
                    if (! bf[lb]) {
                        bf[lb] = 1;
                        bl[ bwd_csa.psi[lb] ] = 1;
                    }
                    if (!bf[rb]) {
                        bf[rb] = 1;  // rb has to be marked for the correctness
                        if (rb < bwd_csa.size()) {
                            bl[ bwd_csa.psi[rb] ] = 1;
                        }
                    }
//					cout << "range=[" << lb <<", " << rb << "]" << endl;
                    if (rb-lb > b) {
                        bwd_csa.wavelet_tree.interval_symbols(lb, rb, k, cs, rank_c_lb, rank_c_rb);
                        for (size_type i=0; i<k; ++i)
                            pq.push(tPCI(cs[i], i));    // push character and its index to pq to sort it
                        // this step would be not necessary if we use a sorted WT

                        for (size_type i=0, j, lex_smaller=0; i<k; ++i) {
                            j = pq.top().second; pq.pop();
                            size_type c_begin = bwd_csa.C[ bwd_csa.char2comp[cs[j]] ];
                            size_type lb_new = c_begin + rank_c_lb[j], rb_new =c_begin + rank_c_rb[j];
//							if( !( bf[lb_new] and bf[rb_new] ) ){ // not a good idea, since it is possible that
                            // one interval could be reached more than one times on different depths
                            q.push(tPII(lb_new, rb_new));
                            q_fp.push(lb_forward + lex_smaller);
                            q_depth.push(depth+1);
//							}

                            lex_smaller += rank_c_rb[j]-rank_c_lb[j];
                        }
                    } else {
                        if (m_n < 100) {
//							cout << "backward interval: " << depth << "-["<< lb << "," << rb << "]" << endl;
//							cout << " forward interval: " << depth << "-["<< lb_forward << ","
//								 << lb_forward + rb-lb  << "] equals block " << fwd_bf_rank(lb_forward+1)-1
//								 << endl;
                            fwd_blocks_in_bwd.push_back(bwd_block_info(lb, rb, depth, fwd_bf_rank(lb_forward+1)-1));
                        }
                        map_info.push_back(block_info(lb, depth, lb_forward, rb-lb));
                    }
                }
                sort(map_info.begin(),map_info.end());
                cout << "finished "<< endl;
                cout << "assign1" << endl;
                util::assign(m_bl, bl);
                cout << "assign2" << endl;
                util::assign(m_bf, bf);
                m_bf_rank.init(&m_bf);
                m_bf_select.init(&m_bf);
                m_bm = bit_vector(m_k     + m_bf_rank(m_bf.size()), 0);
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
                m_bm_select.init(&m_bm);
                m_bm_rank01.init(&m_bm);
                m_min_depth.set_int_width(bit_magic::l1BP(max_min_depth)+1);
                m_min_depth.resize(m_bm_rank01(m_bm.size()));
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

                // allocate pointers; each pointer with 64 bits
                m_pointer.set_int_width(64);
                m_pointer.resize(m_k);

                // calculate external information
                for (size_t bwd_id=0; bwd_id <map_info.size(); ++bwd_id) {
                    size_type fwd_lb = map_info[bwd_id].fwd_lb;
                    size_type fwd_id = fwd_bf_rank(fwd_lb+1)-1;
                    v_block[fwd_id].bwd_id = bwd_id;
                    if (map_info[bwd_id].size == 1) {
                        m_pointer[bwd_id] =  fwd_cst.csa[fwd_lb];
                    }
                }

                for (size_t fwd_id=0; fwd_id < v_block.size(); ++fwd_id) {
                    size_type bwd_id = v_block[fwd_id].bwd_id;
                    size_type dest_block = v_block[fwd_id].dest_block;
                    size_type delta_x = v_block[fwd_id].delta_x;
                    size_type delta_d = v_block[fwd_id].delta_d;
                    if (map_info[bwd_id].size > 1) {
                        external_block[dest_block].push_back(header_item(bwd_id,delta_x,delta_d));
                    }
                }
            }

            {
                size_type total_header_in_bytes = 0;
                size_type ext_idx_size_in_bytes = 0;
                std::ofstream ext_idx_out(get_ext_idx_filename().c_str(), std::ios_base::trunc);
                if (ext_idx_out) {
//TODO: test if the size of the output buffer chances the performance
#ifdef USE_CUSTOM_BUFFER
                    const size_type output_buf_size = 20 * (1ULL<<20);
                    char* output_buf = new char[output_buf_size];
                    ext_idx_out.rdbuf()->pubsetbuf(output_buf, output_buf_size);
#endif
                    select_support_mcl<> fwd_bf_select(&fwd_bf);
                    vector<size_type> block_addr(m_k+1, 0);
                    char* text = NULL;
                    file::read_text(m_file_name.c_str(), text);
                    for (size_t fwd_id=0; fwd_id < external_block.size(); ++fwd_id) {
                        if (external_block[fwd_id].size() > 0) {
                            sort(external_block[fwd_id].begin(), external_block[fwd_id].end());
                            disk_block db;
                            db.set_header(external_block[fwd_id]);
                            db.set_content(fwd_cst, (const unsigned char*)text, fwd_bf_select(fwd_id+1), fwd_bf_select(fwd_id+2)-1);
                            block_addr[fwd_id] = ext_idx_size_in_bytes;
                            if (m_n < 50) {
                                cout << "fwd_id="<<fwd_id<<" block_addr=" <<block_addr[fwd_id] << endl;
                            }
                            ext_idx_size_in_bytes += db.serialize(ext_idx_out);
                            total_header_in_bytes += db.header_size_in_bytes();
                        }
                    }
                    ext_idx_out.close();
                    delete [] text;
                    cout << "# ext_idx_size_in_MB = " << ext_idx_size_in_bytes/(1024.0*1024.0) << endl;
                    cout << "# ext_idx_header_size_in_MB = " << total_header_in_bytes / (1024.0*1024) << endl;
                    for (size_t fwd_idx = 0; fwd_idx < v_block.size(); ++fwd_idx) {
                        size_type lb = fwd_bf_select(fwd_idx+1);
                        if (!fwd_bf[lb+1]) {  // if not a singleton interval
                            size_type addr = block_addr[v_block[fwd_idx].dest_block];
//	cout << "fwd_id="<<fwd_idx<<" bwd_id = "<<v_block[fwd_idx].bwd_id<<" addr="<<addr<<endl;
                            m_pointer[ v_block[fwd_idx].bwd_id ] = addr;
                        }
                    }
#ifdef USE_CUSTOM_BUFFER
                    delete [] output_buf;
#endif
                } else {
                    std::cerr << "ERROR: Could not open file " << get_ext_idx_filename() << endl;
                }
            }
            util::bit_compress(m_pointer);
            if (output_tikz) {
                if (check_size(m_n)) {
                    std::ofstream tikzout((path+"."+util::to_string(b)+".fwd_idx.tex").c_str());
                    if (tikzout) {
                        write_tikz_output_fwd(tikzout, fwd_cst, fwd_bf, b, v_block, external_block);
                        tikzout.close();
                    }
                }
            }

            if (output_tikz) {
                if (check_size(m_n)) {
                    std::ofstream tikzout((path+"."+util::to_string(b)+".bwd_idx.tex").c_str());
                    if (tikzout) {
                        sort(fwd_blocks_in_bwd.begin(), fwd_blocks_in_bwd.end());
                        write_tikz_output_bwd(tikzout, bwd_csa, m_bf, m_bl,m_bm, b, fwd_blocks_in_bwd, m_min_depth);
                        tikzout.close();
                    }
                }
            }
            m_bl_rank.init(&m_bl);
            size_type n1 = m_bl_rank(m_n);
#ifdef OUTPUT_STATS
            cout << "#1_in_bf = " << m_bf_rank(m_bf.size()) << endl;
            cout << "#1_in_bl = " << n1 << endl;
#endif
            m_C1 = int_vector<64>(256,0);
            for (size_type i=0; i<bwd_csa.sigma; ++i) {
                m_C1[bwd_csa.comp2char[i]] = m_bf_rank(bwd_csa.C[i]);
                std::cerr<<"m_C1["<< bwd_csa.comp2char[i] << "]="	<< m_C1[bwd_csa.comp2char[i]]  << std::endl;
            }
            // construct condensed BWT (cBWT)
            string tmp_file_name =  tmp_file_dir2+util::basename(m_file_name)+"_cBWT_"+util::to_string(util::get_pid())+"_"+util::to_string(util::get_id());
            {
                int_vector<8> temp_bwt(n1);
                for (size_type i=0,j=0; i<m_n; ++i) {
                    if (m_bl[i]) {
                        temp_bwt[j++] = bwd_csa.bwt[i];
                    }
                }
                util::store_to_file(temp_bwt, tmp_file_name.c_str());
            }
            int_vector_file_buffer<8, size_type> temp_bwt_buf(tmp_file_name.c_str());
            m_wt.construct(temp_bwt_buf, n1);
            std::remove(tmp_file_name.c_str()); // remove file of BWT'

            open_streams();
        }

        //! Calculate the id of a block in the backward index using its left bound and depth.
        /*!
         *	\param lb		The left bound of the block in the backward index.
         *  \param depth	The depth of the block.
         *  \return	The backward id of the block (in [0..k-1]).
         *
         *	\par Time complexity
         *		 \f$ \Order{1} \f$
         */
        size_type get_bwd_id(size_type lb, size_type depth) const {
            size_type run_nr = m_bf_rank(lb);
            size_type run_pos = 0;
            if (run_nr > 0) {
                run_pos = m_bm_select(run_nr)+1;
            }
            // optimization for the special case when the run is of size 1
            // m_mb[run_pos-1..run_pos] = 10
            if (m_bm[run_pos+1]) {  //i.e. m_mb[run_pos-1..run_pos+1] = 101
                // in the first (run_pos+1) position there are run_nr ones
                // => there are X=(run_pos+1)-run_nr zeros
                // => the id is X-1
                return run_pos - run_nr;
            } else { // i.e. m_mb[run_pos-1..run_pos+1] = 100...
                size_type min_depth = m_min_depth[m_bm_rank01(run_pos+1)];
                return run_pos - run_nr + (depth-min_depth);
            }
        }

        //! Count the number of occurrences of a pattern of length m
        /*!
         *  \param pattern 	A pointer to the start of the pattern.
         *  \param m		The length of the pattern.
         *  \param
         *  \return 		The number of occurrences of the pattern in the text.
         */
        size_type count(const unsigned char* pattern, size_type m, bool repeated_in_memory_search=false)const {
#ifdef OUTPUT_STATS
            ++m_count_queries;
            m_X = gl_X;
            if (m_count_queries == m_X) {
                cout<<"match pattern "<<pattern<<" m="<<m<<" m_count_queries="<< m_count_queries <<endl;
            }
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
                        if (0 == match_pattern(pattern+d, m-d, sa+d))
                            return 1;
                    } else {
                        size_type block_addr = m_pointer[bwd_id];
                        return search_block(pattern+d, m-d, d, rb+1-lb, bwd_id, block_addr);
                    }
                }
            }
            return 0;
        }

        // TODO: implement locate, take care of the case where the pattern does not reach the
        // external index
        //void locate()

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
            const unsigned char* cp = pattern;//+m-1;
            while (d < m and rb-lb+1 > m_b) {
#ifdef ROSA_DEBUG
                std::cout<<d<<"-["<<lb<<","<<rb<<"]"<<" b="<<m_b<<std::endl;
                std::cout<<"c="<<c<<" m_C1["<<*cp<<"]="<<m_C1[*cp]<<std::endl;
#endif
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
                lb = m_bf_select(m_C1[c] + lb2 + 1);
                rb = m_bf_select(m_C1[c] + rb2 + 1) - 1;
#ifdef OUTPUT_STATS
                ++m_count_int_steps;
#endif
            }
            return true;
        }

        //! Search for a pattern by binary search in a block
        /*!
         * \param pattern    Pointer to the remaining pattern.
         * \param m			 Number of not yet matched characters. Equals the length of the remaining pattern.
         * \param depth 	 Number of characters matched so far.
         * \param size		 Size of the lexicographic interval of the already matched prefix.
         * \param bwd_id	 bwd_id of the matched prefix.
         * \param block_addr Start address of the block in the external index.
         *
         * \par Time complexity
         *		\f$ \Order{\log b} \f$ disk accesses
         */
        size_type search_block(const unsigned char* pattern, size_type m, size_type depth, size_type size,
                               size_type bwd_id, size_type block_addr)const {
            disk_block db;
            seekg(m_ext_idx, block_addr);
            db.load(m_ext_idx);
            size_type delta_x = 0;
            size_type delta_d = 0;
            db.get_delta_x_and_d(bwd_id, delta_x, delta_d);
            // do a binary search on the block for the pattern
            // possible TODO: if match_pattern also returns the length of the length X of the longest matching prefix
            //       of pattern with the text, then we can further speed up the computation
            //       by limit the interval to the range where LCP[i] >= X
            size_type lb = delta_x, rb = delta_x+size; // [lb..rb)
            while (lb < rb) {
                size_type mid = (lb+rb)/2;
                int cmp =  match_pattern(pattern, m, db.sa[mid] + depth + delta_d);
                if (0 == cmp) {
                    size_type l=mid, r=mid;
                    while (l > lb and db.lcp[l] >= (depth+m)+delta_d) {
                        --l;
                    }
                    while (r+1 < rb and db.lcp[r+1] >= (depth+m)+delta_d) {
                        ++r;
                    }
                    return r-l+1;
                } else if (cmp < 0) { // pattern is smaller than current suffix
                    rb = mid;
                } else { // pattern is greater than current suffix
                    lb = mid+1;
                }
            }
            return 0;
        }

        //! Check if pattern is a prefix of the suffix starting at position text_offset in T.
        /*!
         *  \param pattern		A pointer to the start of the pattern.
         *	\param m			Length of the pattern.
         *  \param text_offset	Starting position in the text.
         *  \sa match_pattern
         *	\par Time complexity
         *		\f$  \Order{m} \f$ steps and one disk access
         */
        int match_pattern(const unsigned char* pattern, size_type m, size_type text_offset)const {
            size_type matched = 0;
            return match_pattern(pattern, m, text_offset, matched);
        }

        //! Check if pattern is a prefix of the suffix starting at text_offset in T.
        /*! \param pattern		A pointer to the start of the pattern.
         *  \param m			Length of the pattern.
         *  \param text_offset	Starting position in the text.
         *  \return				0 if the pattern is a prefix of the suffix starting
         *                      at text_offset in T, a value smaller 0 if the pattern is lex.
         *                      smaller and greater zero if the pattern is lex. greater
         *                      than the suffix starting at text_offset in T.
         *	\par Time complexity
         *		\f$  \Order{m} \f$ steps and one disk access
         */
        int match_pattern(const unsigned char* pattern, size_type m, size_type text_offset, size_type& matched)const {
            matched = 0;
            if (0==m or text_offset >= m_n-1) {  // avoid disk access if m==0
                return 1;
            }
            bool truncated = false;
            // should be: text_offset+m <= m_n-1
            if (text_offset + m > m_n-1) {   // m_n-1 text size without the sentinel character
                m = (m_n-1) - text_offset;
                truncated = true;
            }
            seekg(m_text, text_offset);
            while (m > 0) {
                size_type len = std::min(m, m_buf_size);
                m_text.read((char*)m_buf, len);
                int res;
                if (0 != (res=memcmp(pattern, m_buf, len))) {
                    const unsigned char* p = pattern, *t = m_buf;
                    while (*(p++)==(*t++)) {
                        ++matched;
                    }
                    return res;
                }
                m -= len;
                pattern += len;
                matched += len;
            }
            if (truncated) {
                return 1; // we assume that the character behind the text is smaller than all other
            } else {
                return 0;
            }
        }



        void reset_counters() {
#ifdef OUTPUT_STATS
            m_count_disk_access 	= 0;
            m_count_gap_disk_access = 0;
            m_count_int_steps   	= 0;
            m_count_int_match   	= 0;
            m_count_queries			= 0;
#endif
        }

        //! Get the size of the external part in megabytes.
        double get_ext_idx_size_in_mega_byte() {
            return ((double)get_file_size(get_ext_idx_filename().c_str()))/ (1024.0*1024.0);
        }

        //! Writes the in-memory part of the index into the output stream.
        /*! \param out Output stream to which the index should be written.
         */
        size_type serialize(std::ostream& out)const {
            size_type written_bytes = 0;
            written_bytes += util::write_member(m_n, out);
            written_bytes += util::write_member(m_b, out);
            written_bytes += util::write_member(m_k, out);
            written_bytes += m_bl.serialize(out);
            written_bytes += m_bf.serialize(out);
            written_bytes += m_bl_rank.serialize(out);
            written_bytes += m_bf_rank.serialize(out);
            written_bytes += m_bf_select.serialize(out);
            written_bytes += m_wt.serialize(out);
            written_bytes += m_C1.serialize(out);
            written_bytes += m_bm.serialize(out);
            written_bytes += m_bm_select.serialize(out);
            written_bytes += m_bm_rank01.serialize(out);
            written_bytes += m_min_depth.serialize(out);
            written_bytes += m_pointer.serialize(out);
            written_bytes += util::write_member(m_file_name, out);
            written_bytes += util::write_member(m_output_dir, out);
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
            m_wt.load(in);
            m_C1.load(in);
            m_bm.load(in);
            m_bm_select.load(in, &m_bm);
            m_bm_rank01.load(in, &m_bm);
            m_min_depth.load(in);
            m_pointer.load(in);
            util::read_member(m_file_name, in);
            util::read_member(m_output_dir, in);
            open_streams();
        }

        void set_file_name(const std::string& file_name) {
            m_file_name = file_name;
            open_streams();
        }

        void set_output_dir(const string& output_dir) {
            m_output_dir = output_dir;
            open_streams();
        }

#ifdef MEM_INFO
        void mem_info(std::string label="")const {
            if (label=="") {
                label = "rosa";
            }
            double megabytes = util::get_size_in_mega_bytes(*this);
            std::cout << "list(label=\""<<label<<"\", size="<< megabytes << "\n,";
            m_bl.mem_info("bl"); std::cout<<",";
            m_bf.mem_info("bf"); std::cout<<",";
            m_bl_rank.mem_info("bl_rank"); std::cout<<",";
            m_bf_rank.mem_info("bf_rank"); std::cout<<",";
            m_bf_select.mem_info("bf_select"); std::cout<<",";
            m_wt.mem_info("wt"); std::cout<<",";
            m_C1.mem_info("C1"); std::cout<<",";
            m_bm.mem_info("bm"); std::cout<<",";
            m_bm_select.mem_info("bm_select"); std::cout<<",";
            m_bm_rank01.mem_info("bm_rank01"); std::cout<<",";
            m_min_depth.mem_info("min_depth"); std::cout<<",";
            m_pointer.mem_info("pointer");
            std::cout << ")\n";
        }
#endif

        //! Output statistics about the data structure
        // Maybe TODO: output the number of nodes of the corresponding trie
        void statistics()const {
            size_type k_ir = 0; // counter for irreducible blocks
            size_type k_re = 0; // counter for reducible blocks
            size_type k_s  = 0; // counter for singleton blocks

            size_type n_ir = 0; // number of elements in irreducible blocks

            double header_in_megabyte = 0.0;
            double lcp_in_megabyte = 0.0;
            double sa_in_megabyte  = 0.0;

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
        }

        //! Print statistics about the input text.
        void text_statistics(size_type max_k=0) {
            std::string path = m_output_dir + "/" + util::basename(file_name) ;
            string tmp_fwd_cst_file_name = path + "." + TMP_CST_SUFFIX;
            tCst cst;
            if (!util::load_from_file(cst, tmp_fwd_cst_file_name.c_str())) {
                std::cerr << "ERROR: could not open the compressed suffix tree file" << std::endl;
            } else {
                std::cout << "k   H_k         contexts" << std::endl;
                for (size_type k=0, contexts=0; k<=max_k; ++k) {
                    double hk = Hk(cst, k, contexts);
                    std::cout << std::setw(4) << k <<" "<< std::setw(10) << hk << " "<<contexts << std::endl;
                }
            }
        }

        //! Output the number of nodes in the in-memory trie of LOF-SA
        void trie_nodes() const {
            std::string path = m_output_dir + "/" + util::basename(file_name) ;
            string tmp_fwd_cst_file_name = path + "." + TMP_CST_SUFFIX;
            tCst cst;
            if (!util::load_from_file(cst, tmp_fwd_cst_file_name.c_str())) {
                std::cerr << "ERROR: could not open the compressed suffix tree file" << std::endl;
            } else {
                size_type trie_nodes;
                bit_vector dummy_bf(m_bf.size());
                size_type dummy_cnt;
                get_block_info_and_mark_blocks(cst, dummy_bf, m_b, trie_nodes, dummy_cnt);
                std::cout << "# trie_nodes = " << trie_nodes << std::endl;
            }
        }
};

#endif // end of include guard
