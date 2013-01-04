#ifndef ROSA_PATTERN_FILE
#define ROSA_PATTERN_FILE

#include "bu_interval.hpp"
#include <sdsl/int_vector.hpp> // for bit_vector
#include <sdsl/util.hpp> // for read_member
#include <cstdlib> // for rand()
#include <iostream>
#include <fstream>
#include <vector>

using namespace sdsl;
using namespace std;

/*! A helper class for producing and reading patterns
 *  
 */
class pattern_file{
	public:
		typedef bit_vector::size_type size_type;
	private:
		ifstream m_pattern_stream;
		size_type m_pattern_cnt;
		size_type m_pattern_len;
		size_type m_swaps;
		char *m_buf;
		string m_pattern_file_name;
	public:
		size_type &pattern_cnt;
		size_type &pattern_len;
		size_type &swaps;

	//! Constructor takes the file name where the pattern should be written to or read from	
	pattern_file(const char *pattern_file_name);

	//! Opens the pattern file; after the call to rest we get the next pattern be calling get_next_pattern
	/*!
	 * \sa get_next_pattern
	 */
	void reset();

	//! Generates a number of patterns from a given text and store them in a file 
	/* 
	 * \param text_file_name		File name of the text.
	 * \param pattern_cnt			Number of pattern that should be generated.
	 * \param pattern_len			Length of each pattern that should be generated.
	 * \param swaps					The number of random neighbour character swaps applied to the pattern.
	 *
	 * \par Details of the generation
	 *      We generate a pattern by randomly jump to a position \f$i\f$ in the text T and then
	 *      we extract the substring \f$P=T[i..i+pattern_len-1]\f$. After the extraction
	 *      we chose \f$swaps\f$ times a random position \f$r\f$ in [0..pattern_len-2] and swap characters
	 *      \f$P[r]\f$ and \$fP[r+1]\f$.
	 */	
	void generate(const char *text_file_name, size_type pattern_cnt, size_type pattern_len, size_type swaps=0);

	//! Generates a number of patterns with some restrictions (length and occurrences) 
	/*
	 * \param cst 					Compressed suffix tree of the text.
	 * \param text_file_name 		File name of the text.
	 * \param pattern_cnt			Number of pattern that should be generated.
	 * \param pattern_len			Length of each pattern that should be generated.
	 * \param min_occ				Minimum number of occurrences of the pattern.
	 * \param max_occ				Maximum number of occurrences of the pattern.
	 */
	template<class tCst>
	void generate_restricted(const tCst &cst, const char *text_file_name, size_type pattern_cnt, 
							 size_type pattern_len, size_type min_occ, size_type max_occ);

	void generate_restricted(const char* lcp_file, const char *sa_file, const char *text_file_name, size_type pattern_cnt, 
							 size_type pattern_len, size_type min_occ, size_type max_occ);
	//! Destructor
	~pattern_file();

	//! Get the next pattern from the pattern file
	/*! Before a call to this function we have to call the init function.
	 * \sa init 
	 */
	const char* get_next_pattern();

	//! Pointer to current pattern.
	const char* get_pattern();

	//! Processing function for  the bottom-up traversal. 
	void get_candidate(bu_interval *v, size_t pattern_len, size_t min_occ, size_t max_occ, std::vector<size_t> &candidates);

	//! Remove the associated pattern file from disk.
	void remove();
};

template<class tCst>
void pattern_file::generate_restricted(const tCst &cst, const char *text_file_name, size_type pattern_cnt, 
						 size_type pattern_len, size_type min_occ, size_type max_occ){
	typedef typename tCst::node_type node_type;
	m_pattern_cnt = pattern_cnt;
	m_pattern_len = pattern_len;
	m_swaps		  = 0;

	vector<size_type> candidates;
	for (typename tCst::const_iterator it = cst.begin(), end=cst.end(); it!=end; ++it){
		if ( it.visit() == 1 ){
			node_type v = *it;
			if ( cst.depth(v) >= pattern_len ){ // if the depth >= pattern_len we determine
										// if min_occ <= occ <= max_occ
				size_type occ = cst.leaves_in_the_subtree( v );
				if ( min_occ <= occ and occ <= max_occ ){
					candidates.push_back( cst.lb(v) );
				}
				it.skip_subtree();
			}else{ // if the depth < pattern_len we process the subtree
				// if the subtree size is to small we don't process the subtree
				if ( min_occ > cst.leaves_in_the_subtree(v) ){
					it.skip_subtree();
				}
			}
		}
	}
	bit_vector already_used(candidates.size(), 0);
	std::ofstream pattern_stream( m_pattern_file_name.c_str() );
	if ( pattern_stream ){
		if ( candidates.size() == 0 ){
			m_pattern_cnt = 0;
			std::cerr << "Found no candidates for pattern_len = "<< pattern_len << std::endl;
			std::cerr << "min_occ = "<< min_occ << " and max_occ = " << max_occ << std::endl;
		}else if( candidates.size() < 100 ){
			std::cerr << "Less then 100 candidates for pattern_len = " << pattern_len << std::endl;
			std::cerr << "min_occ = "<< min_occ << " and max_occ = " << max_occ << std::endl;
		}else{
			std::cerr << candidates.size() << " candidates for pattern_len = " << pattern_len << std::endl;
			std::cerr << "min_occ = "<< min_occ << " and max_occ = " << max_occ << std::endl;
		}
		util::write_member(m_pattern_cnt, pattern_stream); 
		util::write_member(m_pattern_len, pattern_stream); 
		util::write_member(m_swaps, pattern_stream); 
		if ( candidates.size() > 0 )
		{
			char *pattern_buf = new char[pattern_len+1];
			pattern_buf[pattern_len] = '\0';
			ifstream text_stream(text_file_name);	
			if ( text_stream ){
				for (size_type i=0, j, k; i < pattern_cnt; ++i){
					j=rand()%candidates.size();
					if ( already_used[j] ) {
						// find next free location 
						size_type j1=j+1;
						while ( j1 < already_used.size() and already_used[j1] ){
							++j1;
						}
						if ( j1 >= already_used.size() ){
							j1 = j;
						}
						j = j1;
					}
					k = cst.csa[ candidates[j] ];
					// extract the pattern
					text_stream.seekg(k ,std::ios::beg);
					text_stream.read(pattern_buf, pattern_len);
// cout << i << " " << pattern_buf <<endl;
					pattern_stream.write( pattern_buf, pattern_len );
					already_used[j] = 1;
				}
				text_stream.close();
			}else{
				std::cerr << "ERROR: Could not open text file " << text_file_name << endl;	
			}
			pattern_stream.close();
		}
	}else{
		std::cerr << "ERROR: Could not open pattern file " << m_pattern_file_name << endl;
	}
}

#endif
