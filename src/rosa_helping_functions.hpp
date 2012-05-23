#ifndef ROSA_HELPING_FUNCTIONS
#define ROSA_HELPING_FUNCTIONS

#include "rosa_helping_structures.hpp"

#include <sdsl/rank_support_v.hpp>
#include <sdsl/select_support_mcl.hpp>
#include <sdsl/util.hpp>

#include <stack>
#include <string>
#include <vector>

using std::stack;
using std::vector;

using namespace sdsl;

typedef bit_vector::size_type size_type;

//! Calculates for a block how many fringe symbols are needed to full cover it
/* 
 * \param cst The compressed suffix tree of the input text.  
 * \param v The block as suffix tree node.
 * \param p The parent node of the block.
 * \param fringe_len The length of the fringe
 * \return The number of fringe symbols that a needed to full cover the block,
 *         i.e. that only one additional disk access is needed to answer a
 *         existence query
 */
template<class tCst>
size_type get_size_of_full_fringe(const tCst &cst, typename tCst::node_type v, typename tCst::node_type p, size_type fringe_len=1){
	size_type size = 0;
	// after the search in the in-memory data structure we
	// know that the prefixes in the block have d=cst.depth(p)+1
	// characters in common
	// Note that there can be a gap between the minimal 
	// value in LCP[lb..rb] and d!
	size_type known_prefix_len = cst.depth(p)+1;
	size_type lb = cst.lb(v), rb = cst.rb(v);
	for(size_type i=lb+1; i <= rb; ++i){
		if( cst.lcp[i] > known_prefix_len ){
			size += (cst.lcp[i] - known_prefix_len);
		}
		known_prefix_len = cst.lcp[i] + fringe_len;
		size += fringe_len;	
	}
	return size;
}

//! Mark each block [lb..rb] with a one at positions lb and rb+1 in bf and calculate statistics
/* \param cst A reference to the compressed suffix tree of the input text.
 * \param bf  A reference to the bit_vector which is marked with the block boundaries.
 * \param b   The maximal block size / block threshold.
 * \param trie_nodes A reference which will contain the number of nodes of the in-memory trie.
 * \param blocks A reference which will contain the number of blocks.
 * \param size_of_full_fringe Returns the number of fringe characters which have to be stored
 *                            in order to answer all queries with only one additional disk access
 */
template<class tCst>
void get_block_info_and_mark_blocks(const tCst &cst, bit_vector &bf, size_type b, size_type &trie_nodes,
									size_type &blocks, size_type &size_of_full_fringe){
	typedef typename tCst::node_type node_type;
	trie_nodes = 1; // initialize with the root
	blocks     = 0;
	size_of_full_fringe = 0;
	for (typename tCst::const_iterator it = cst.begin(), end=cst.end(); it!=end; ++it){
		if ( it.visit() == 1 ){
			node_type v = *it;
			node_type p = cst.parent(v);
			if ( cst.leaves_in_the_subtree(v) <= b ){
				++trie_nodes; // count the leaf node of the trie
				++blocks; // count the block
				bf[cst.lb(v)] = 1; // mark the beginning of the block
				size_of_full_fringe += get_size_of_full_fringe(cst, v, p, 1);
				it.skip_subtree();
			}
			else{ // count the inner nodes of the trie 
				// each character on the edge from p to v gets a node
				trie_nodes += cst.depth(v) - cst.depth(p);	
			}
		}
	}
}

bool check_size(size_type size);

template<class tCst>
size_type get_block_depth(const tCst &cst, size_type lb, size_type rb){
	size_type result = cst.lcp[lb]+1;
	if ( rb+1 < cst.csa.size() and cst.lcp[rb+1]+1 > result ){
		result = cst.lcp[rb+1]+1;
	}
	return result;
}

//! Get more info about the block. See above. get_block_info_and_mark_blocks
/*
 *  \param cst The compressed suffix tree of the input text.
 *  \param bf  A bit vector which marks each block [lb..rb] at positions lb and rb+1. 
 *  \param b   Threshold for the block size.  
 *  \param red_blocks 			The number of reducible blocks.
 *  \param red_blocks_delta0	The number of reducible blocks which point
 *                              
 */
template<class tCst>
void get_block_info(const tCst &cst, const bit_vector &bf, size_type b, 
		            size_type &red_blocks, size_type &red_blocks_delta0, 
		            size_type &singleton_blocks, size_type &elements_in_irred_blocks, 
					size_type &size_of_full_fringe, vector<block_node> &v_block){
	typedef typename tCst::node_type node_type;
	red_blocks = 0;
	red_blocks_delta0 = 0;
	singleton_blocks = 0;
	elements_in_irred_blocks = 0;
	size_of_full_fringe = 0;

	rank_support_v<> bf_rank(&bf);
	select_support_mcl<> bf_select(&bf);
	v_block = vector<block_node>( bf_rank( cst.csa.size() ) );

	size_type block_id = 0;
	for (typename tCst::const_iterator it = cst.begin(), end=cst.end(); it!=end; ++it){
		if ( it.visit() == 1 ){
			node_type v = *it;
			node_type p = cst.parent(v);
			if ( cst.leaves_in_the_subtree(v) <= b ){
				it.skip_subtree();
				unsigned char c1;
				size_type lb = cst.lb(v), rb = cst.rb(v) + 1;
				size_type rank1 = cst.csa.wavelet_tree.rank_ith_symbol(lb, c1);
				size_type rank2 = cst.csa.rank_bwt(rb, c1);
				if( rank2-rank1 == rb-lb ){
					if( rb-lb == 1 ){
						++singleton_blocks;
						v_block[block_id] = block_node(0, 0, block_id);	
					}
					else{
						++red_blocks;
						// determine who many times we can reduce the block until
						// we reach an irreducible block
					
						stack<size_type> cur_block_id; 
						stack<size_type> cur_delta_x; 
						size_type dest_block_id = block_id;

						do {
							cur_block_id.push( dest_block_id );
							dest_block_id = bf_rank(cst.csa.psi(lb)+1)-1;
							size_type delta_x = cst.csa.psi(lb) - bf_select(dest_block_id+1);
							cur_delta_x.push(delta_x);

							lb = bf_select(dest_block_id+1);
							rb = bf_select(dest_block_id+2);
							rank1 = cst.csa.wavelet_tree.rank_ith_symbol(lb, c1);
							rank2 = cst.csa.rank_bwt(rb, c1);	
						} //   v_block[dest_block] not already defined and block is reducible
						while ( v_block[dest_block_id].dest_block == 0 and rank2-rank1 == rb-lb );
						
						if( rank2 - rank1 != rb-lb ) { // if 
							v_block[dest_block_id].dest_block = dest_block_id;
						}
						size_type delta_d = v_block[dest_block_id].delta_d;
						size_type delta_x = v_block[dest_block_id].delta_x;
						dest_block_id   = v_block[dest_block_id].dest_block;
						while( !cur_block_id.empty() ){
							size_type x = cur_block_id.top(); cur_block_id.pop();
							++delta_d;
							delta_x += cur_delta_x.top(); cur_delta_x.pop();
							v_block[ x ] = block_node(delta_x, delta_d, dest_block_id);
						}
					}
				}else{
					size_of_full_fringe += get_size_of_full_fringe(cst, v, p, 1);
					elements_in_irred_blocks += (rb-lb);
					v_block[block_id] = block_node(0, 0, block_id);
				}
				++block_id;
			}
		}
	}
}

//! Outputs the maximal delta_x and delta_d and the average delta_x and delta_d 
/*
 * \param v_block A vector of block_node items
 */
void output_block_info_statistics(const vector<block_node> &v_block);

//! Close a stream if it is already opened
void close_stream_if_open(std::ifstream &in);

//! Returns the output directory
/*
 *	\param file_name	File name of the original text file
 *	\param output_dir	Char pointer to the output directory 
 *	\return If output_dir==NULL the dirname(file_name) will be returned and
 *          output_dir otherwise.
 */
std::string get_output_dir(const char *file_name, const char *output_dir=NULL);

#endif
