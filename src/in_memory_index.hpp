/*! \file in_memory_index.hpp
    \author Simon Gog (simon.gog@unimelb.edu.au)
 */
#ifndef SDSL_IN_MEMORY_INDEX
#define SDSL_IN_MEMORY_INDEX

#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp> // for basename
#include <sdsl/algorithms_for_string_matching.hpp> // for backward search
#include <sdsl/csa_construct.hpp>
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

namespace sdsl{

//! A wrapper class for a compressed suffix array that provides the same interface as the rosa class
/*
 * The template parameter specifies the compressed suffix array
 */
template<class Csa> 
class in_memory_index{
	public:
		typedef bit_vector::size_type size_type;
		typedef void* tCst;
	private:
		Csa m_csa;	
		size_type m_b;
		string m_file_name;
		mutable size_type m_count_disk_access;
		mutable size_type m_count_gap_disk_access;
		mutable size_type m_count_int_steps;   		
		mutable size_type m_count_int_match;   	
		mutable size_type m_count_queries;	   

	public:
		const string &file_name;
		const int fringe;
		const size_type fringe_len;
		const size_type &count_disk_access;
		const size_type &count_gap_disk_access;
		const size_type &count_int_steps;
		const size_type &count_int_match;
		const size_type &count_queries;	
		const string tmp_cst_suffix;

		//! Constructor
		/*!
		 *	\param file_name 		File name of the file which contains the text.
		 *  \param b				Threshold parameter b for the maximum size of a block. 
		 *	\param fringe_type		This parameter is not used.
		 *	\param fringe_length	This parameter is not used.
		 *	\param output_tikz		This parameter is not used.
		 *	\param delete_tmp		Indicates if the temporary files that are used during the construction should be removed after the construction.
		 *  \param tmp_file_dir		Directory for the temporary files.
		 *	\param output_dir		Directory for the output.
		 *
		 *	\Time complexity
		 *		\f$ \Order{n \log\sigma} \f$, where n is the length of the text.
		 */
		in_memory_index(const char *file_name=NULL, size_type b=4000, fringe_type x=NO_FRINGE, size_type fringe_len=0, 
				        bool output_tikz=false, bool delete_tmp=false, const char *tmp_file_dir="./", const char *output_dir=NULL):
			file_name(m_file_name),
			fringe(0),
			fringe_len(0),
			count_disk_access(m_count_disk_access),
			count_gap_disk_access(m_count_gap_disk_access),
			count_int_steps(m_count_int_steps),
			count_int_match(m_count_int_match),
			count_queries(m_count_queries),
			tmp_cst_suffix("")
		{
			if ( file_name == NULL )
				return;
			m_b = b;
			m_file_name = string(file_name);
			std::string tmp_file_dir2 = (util::dirname(tmp_file_dir)+"/"+util::basename(tmp_file_dir)+"/");
			std::ifstream in( get_int_idx_filename().c_str() );
			if ( !in ){
				tMSS file_map;
				construct_csa_of_reversed_text(m_file_name.c_str(), m_csa, file_map, delete_tmp, tmp_file_dir2, ("reversed_"+util::basename(string(file_name))).c_str() );
			}else{
				load(in);
				in.close();
			}
		}

		static string get_int_idx_filename(const char *file_name, size_type b, size_type fringe_len=0, 
										   fringe_type fringe=NO_FRINGE, const char *output_dir=NULL){
			return get_output_dir(file_name, output_dir) + "/" + util::basename(file_name) 
					+ "."+util::to_string(b)+"." + INDEX_SUF +".int_idx";
		}

		string get_ext_idx_filename(){
			return get_int_idx_filename();
		}

		static string get_ext_idx_filename(const char *file_name, size_type b, size_type fringe_len, fringe_type fringe, const char *output_dir=NULL){
			return get_int_idx_filename(file_name, b, fringe_len, fringe, output_dir);
		}

		double get_ext_idx_size_in_mega_byte(){
			return 0.0;
		}

		string get_int_idx_filename()const{
			return get_int_idx_filename(m_file_name.c_str(), m_b);
		}

		static void remove_tmp_files(const char *file_name, const char *output_dir){
			
		}

		//! Count the number of occurrences of a pattern of length m
		/*!
		 *  \param pattern 	A pointer to the start of the pattern.
		 *  \param m		The length of the pattern.
		 *  \return 		The number of occurrences of the pattern in the text.
		 */
		size_type count(const unsigned char *pattern, size_type m, bool dummy)const{
			// count(): search for the reversed pattern, since we use the bwt of the reversed string	
			size_type lb = 0, rb = m_csa.size()-1, d = 0;
			const unsigned char *cp = pattern;//+m-1;
			while( d < m and rb-lb+1  > 0 ){
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

		//! Match the pattern p reversed backwards against the pruned BWT until the interval <= b.
		/*! \param pattern  Pointer to the beginning of the pattern.
		 *  \param m		The length of the pattern.
		 *  \param lb		Left bound of the resulting interval in the backward index.
		 *  \param rb		Right bound of the resulting interval in the backward index (inclusive).
		 *	\param d		Number of matched characters.
		 *	\return 		True if the pattern can occur in T and false otherwise.
		 */
		bool get_interval(const unsigned char *pattern, size_type m, 
				               size_type &lb, size_type &rb, size_type &d)const{
			d  = 0;
			lb = 0; rb = m_csa.size()-1;
			const unsigned char *cp = pattern;//+m-1;
			while( d < m and rb-lb+1 > m_b ){
#ifdef IN_MEMORY_INDEX_DEBUG				
std::cout<<d<<"-["<<lb<<","<<rb<<"] b="<<m_b<<std::endl;				
std::cout<<"c="<<*cp<<std::endl;	
#endif
				size_type lb2, rb2;
				algorithm::backward_search(m_csa, lb, rb, *cp, lb2, rb2);
				if( lb2 == rb2+1 ){// the pattern does not exist in the pruned BWT
					return false;
				}
				++d; // we have match another character
				++cp;
				lb = lb2;
				rb = rb2;
			}
			return true;
		}

		size_type serialize(std::ostream &out)const{
			size_type written_bytes = 0;
			written_bytes += m_csa.serialize(out);
			written_bytes += util::write_member(m_b, out);
			written_bytes += util::write_member(m_file_name, out);
			return written_bytes;
		}	

		void load(std::istream &in){
			m_csa.load(in);
			util::read_member(m_b, in);
			util::read_member(m_file_name, in);
		}

		void set_file_name(const std::string &file_name){
			m_file_name = file_name;
		}

		void set_output_dir(const string &output_dir){
		}

#ifdef MEM_INFO
		void mem_info(std::string label="")const{
			if ( label=="" ){
				label = "in_memory_index";
			}
			double megabytes = util::get_size_in_mega_bytes(*this);
			std::cout << "list(label=\""<<label<<"\", size="<< megabytes << "\n,";
			m_csa.mem_info("csa"); 
			std::cout << ")\n";
		}
#endif

		void reset_counters(){
#ifdef OUTPUT_STATS			
			m_count_disk_access = 0;
#endif			
		}

		void statistics()const{
			std::cout << "# sigma = " << (int)m_csa.sigma << std::endl;
		}

		void text_statistics(size_type max_k=0)const{
		}

		void trie_nodes()const{
		}
};

}// end namespace
#endif
