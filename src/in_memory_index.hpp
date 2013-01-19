/*! \file in_memory_index.hpp
    \author Simon Gog (simon.gog@unimelb.edu.au)
 */
#ifndef SDSL_IN_MEMORY_INDEX
#define SDSL_IN_MEMORY_INDEX

#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp> // for basename
#include <sdsl/algorithms_for_string_matching.hpp> // for backward search
#include <sdsl/construct.hpp>
#include <cstring> // remove? 
#include <string>
#include <fstream>

#include "rosa_helping_structures.hpp"
#include "rosa_helping_functions.hpp"

using std::string;
using std::ifstream;

#ifndef INDEX_SUF
#define INDEX_SUF "_"
#endif

namespace sdsl
{

//! A wrapper class for a compressed suffix array that provides the same interface as the rosa class
/*
 * The template parameter specifies the compressed suffix array
 */
template<class Csa>
class in_memory_index
{
    public:
        typedef bit_vector::size_type size_type;
        typedef void* tCst;
    private:
        Csa m_csa;
        size_type m_b;
        string m_file_name;
		size_type m_k;
		size_type m_lz_width;
        mutable size_type m_count_disk_access;
        mutable size_type m_count_gap_disk_access;
        mutable size_type m_count_int_steps;
        mutable size_type m_count_int_match;
        mutable size_type m_count_queries;
		mutable size_type m_count_block_length;

    public:
        const string& file_name;
        const size_type& count_disk_access;
        const size_type& count_gap_disk_access;
        const size_type& count_int_steps;
        const size_type& count_int_match;
        const size_type& count_queries;
        const size_type& count_block_length;
        const string tmp_cst_suffix;
		const size_type &k;
		const uint8_t &lz_width;

        //! Constructor
        /*!
         *	\param file_name 		File name of the file which contains the text.
         *  \param b				Threshold parameter b for the maximum size of a block.
         *	\param output_tikz		This parameter is not used.
         *	\param delete_tmp		Indicates if the temporary files that are used during the construction should be removed after the construction.
         *  \param tmp_file_dir		Directory for the temporary files.
         *	\param output_dir		Directory for the output.
         *
         *	\Time complexity
         *		\f$ \Order{n \log\sigma} \f$, where n is the length of the text.
         */
        in_memory_index(const char* file_name=NULL, size_type b=4096, size_type f_fac_dens=1, 
                        bool output_tikz=false, bool delete_tmp=false, const char* tmp_file_dir="./", const char* output_dir=NULL):
			k(m_k),
			lz_width(m_lz_width),
            file_name(m_file_name),
            count_disk_access(m_count_disk_access),
            count_gap_disk_access(m_count_gap_disk_access),
            count_int_steps(m_count_int_steps),
            count_int_match(m_count_int_match),
            count_queries(m_count_queries),
            count_block_length(m_count_block_length),
            tmp_cst_suffix("") {
			m_k=0; m_lz_width = 0;
            if (file_name == NULL)
                return;
            m_b = b;
            m_file_name = string(file_name);
            std::string tmp_dir = (util::dirname(tmp_file_dir)+"/"+util::basename(tmp_file_dir)+"/");
            std::ifstream in(get_int_idx_filename().c_str());
            if (!in) {
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
                construct(m_csa, rev_file_name.c_str(), config, 1);
            } else {
                load(in);
                in.close();
            }
        }

		void factor_frequency(){}

        string get_rev_file_name(const string& file_name) {
            return file_name+"_rev";
        }

        static string get_int_idx_filename(const char* file_name, size_type b, size_type fac_dens=0, const char* output_dir=NULL) {
            return get_output_dir(file_name, output_dir) + "/" + util::basename(file_name)
                   +"." + INDEX_SUF +".int_idx";
        }

        string get_ext_idx_filename() {
            return get_int_idx_filename();
        }

        static string get_ext_idx_filename(const char* file_name, size_type b, size_type fac_dens=0, const char* output_dir=NULL) {
            return get_int_idx_filename(file_name, b, fac_dens, output_dir);
        }

        double get_ext_idx_size_in_mega_byte() {
            return 0.0;
        }

        string get_int_idx_filename()const {
            return get_int_idx_filename(m_file_name.c_str(), m_b);
        }

        static void remove_tmp_files(const char* file_name, const char* output_dir) {

        }

        //! Get the name of the file where the factorization of the index is stored.
        string get_factorization_filename() const {
			return "";
        }

	    //! Get the name of the file where the factorization of the index is stored.
        static string get_factorization_filename(const char* file_name, size_type b, const char* output_dir=NULL) {
			return "";
        }	



		size_type size()const{
			return m_csa.size();
		}

		void reconstruct_text(const std::string delimiter=""){
		}

		void output_factors(){}

		void output_tikz(){}

        size_type greedy_parse(string tmp_dir, bit_vector &factor_borders, size_type fac_dens) {
			return 0;
		}

        //! Count the number of occurrences of a pattern of length m
        /*!
         *  \param pattern 	A pointer to the start of the pattern.
         *  \param m		The length of the pattern.
         *  \return 		The number of occurrences of the pattern in the text.
         */
        size_type count(const unsigned char* pattern, size_type m, bool dummy)const {
            // count(): search for the reversed pattern, since we use the bwt of the reversed string
            size_type lb = 0, rb = m_csa.size()-1, d = 0;
            const unsigned char* cp = pattern;//+m-1;
            while (d < m and rb-lb+1  > 0) {
                size_type lb2, rb2;
                algorithm::backward_search(m_csa, lb, rb, *cp, lb2, rb2);
                ++d;
                ++cp;
                lb = lb2; rb = rb2;
#ifdef OUTPUT_STATS
                ++m_count_disk_access; // a potential disk access for every step
#endif
            }
            return rb-lb+1;
        }

        size_type locate(const unsigned char* pattern, size_type m, vector<size_type>& res)const {
	        size_type lb = 0, rb = m_csa.size()-1, d = 0;
            const unsigned char* cp = pattern;//+m-1;
            while (d < m and rb-lb+1  > 0) {
                size_type lb2, rb2;
                algorithm::backward_search(m_csa, lb, rb, *cp, lb2, rb2);
                ++d;
                ++cp;
                lb = lb2; rb = rb2;
#ifdef OUTPUT_STATS
                ++m_count_disk_access; // a potential disk access for every step
#endif
            }
			if ( rb > lb+1 ){
				res.resize(rb+1-lb);
				for(size_type i=lb; i<=rb; ++i){
					res[i-lb] = m_csa[i];
#ifdef OUTPUT_STATS
                m_count_disk_access += m_csa.sa_sample_dens; // a potential disk access for every step
#endif
				}
			}
            return rb-lb+1;
		}

        //! Match the pattern p reversed backwards against the pruned BWT until the interval <= b.
        /*! \param pattern  Pointer to the beginning of the pattern.
         *  \param m		The length of the pattern.
         *  \param lb		Left bound of the resulting interval in the backward index.
         *  \param rb		Right bound of the resulting interval in the backward index (inclusive).
         *	\param d		Number of matched characters.
         *	\return 		True if the pattern can occur in T and false otherwise.
         */
        bool get_interval(const unsigned char* pattern, size_type m,
                          size_type& lb, size_type& rb, size_type& d)const {
            d  = 0;
            lb = 0; rb = m_csa.size()-1;
            const unsigned char* cp = pattern;//+m-1;
            while (d < m and rb-lb+1 > m_b) {
#ifdef IN_MEMORY_INDEX_DEBUG
                std::cout<<d<<"-["<<lb<<","<<rb<<"] b="<<m_b<<std::endl;
                std::cout<<"c="<<*cp<<std::endl;
#endif
                size_type lb2, rb2;
                algorithm::backward_search(m_csa, lb, rb, *cp, lb2, rb2);
                if (lb2 == rb2+1) {// the pattern does not exist in the pruned BWT
                    return false;
                }
                ++d; // we have match another character
                ++cp;
                lb = lb2;
                rb = rb2;
            }
            return true;
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_csa.serialize(out, child, "csa");
            written_bytes += util::write_member(m_b, out, child, "b");
            written_bytes += util::write_member(m_file_name, out, child, "file_name");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in) {
            m_csa.load(in);
            util::read_member(m_b, in);
            util::read_member(m_file_name, in);
        }

        void set_file_name(const std::string& file_name) {
            m_file_name = file_name;
        }

        void set_output_dir(const string& output_dir) {
        }

#ifdef MEM_INFO
        void mem_info(std::string label="")const {
            if (label=="") {
                label = "in_memory_index";
            }
            double megabytes = util::get_size_in_mega_bytes(*this);
            std::cout << "list(label=\""<<label<<"\", size="<< megabytes << "\n,";
            m_csa.mem_info("csa");
            std::cout << ")\n";
        }
#endif

        void reset_counters() {
#ifdef OUTPUT_STATS
            m_count_disk_access = 0;
#endif
        }

        void statistics()const {
            std::cout << "# sigma = " << (int)m_csa.sigma << std::endl;
        }

        void text_statistics(size_type max_k=0)const {
        }

        void trie_nodes()const {
        }
};

}// end namespace
#endif
