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


template<class tCst>
void get_trie_nodes_and_blocks(const tCst& cst, size_type b, size_type& trie_nodes, size_type& blocks)
{
    typedef typename tCst::node_type node_type;
    trie_nodes = 1; // initialize with the root
    blocks     = 0; // initialize blocks
    for (typename tCst::const_iterator it = cst.begin(), end=cst.end(); it!=end; ++it) {
        if (it.visit() == 1) {
            node_type v = *it;
            node_type p = cst.parent(v);
            if (cst.leaves_in_the_subtree(v) <= b) {
                ++trie_nodes; // count the leaf node of the trie
                ++blocks; // count the block
            } else { // count the inner nodes of the trie
                // each character on the edge from p to v gets a node
                trie_nodes += cst.depth(v) - cst.depth(p);
            }
        }
    }
}

//! Mark each block [lb..rb] with a one at positions lb and rb+1 in bf.
/* \param cst A reference to the compressed suffix tree of the input text.
 * \param bf  A reference to the bit_vector which is marked with the block boundaries.
 * \param b   The maximal block size / block threshold.
 * \param blocks A reference which will contain the number of blocks.
 */
template<class tCst>
void mark_blocks(const tCst& cst, bit_vector& bf, size_type b, size_type& blocks) {
    typedef typename tCst::node_type node_type;
    blocks     = 0;
    for (typename tCst::const_iterator it = cst.begin(), end=cst.end(); it!=end; ++it) {
        if (it.visit() == 1) {
            node_type v = *it;
            if (cst.leaves_in_the_subtree(v) <= b) {
                ++blocks; // count the block
                bf[cst.lb(v)] = 1; // mark the beginning of the block
                it.skip_subtree();
            }
        }
    }
}

// bottom up interval
struct bu_interval{
	typedef int_vector<>::size_type size_type;
	size_type lcp, lb, rb;
	std::vector<bu_interval*> children;
	bu_interval(size_type _lb, size_type _rb, size_type _lcp) : lb(_lb), rb(_rb), lcp(_lcp){};
	~bu_interval(){ delete_children(); }
	void delete_children(){
		for (size_t i=0; i<children.size(); ++i){ //std::cout<<"delete "<<children[i]->lcp<<"-["<<children[i]->lb<<","<<children[i]->rb<<"]"<<std::endl;
			delete children[i]; 
		}
		children.resize(0);
	}
	size_type size()const{
		return rb-lb+1;
	}
};

void mark_interval(bu_interval *v, bit_vector& bf, size_type b);
void mark_blocks(const char* lcp_file, bit_vector& bf, size_type b, size_type& blocks);

bool check_size(size_type size);

//! Get more info about the block. See above. get_block_info_and_mark_blocks
/*
 *  \param cst The compressed suffix tree of the input text.
 *  \param bf  A bit vector which marks each block [lb..rb] at positions lb and rb+1.
 *  \param b   Threshold for the block size.
 *  \param red_blocks 			The number of reducible blocks.
 *
 *	\par Time complexity
 *		 \f$  \f$
 *	\par Space complexity
 *
 */
template<class tCsa>
void calculate_reducible_graph(const tCsa& csa, const bit_vector& bf, size_type b,
                               size_type& red_blocks, size_type& singleton_blocks, size_type& elements_in_irred_blocks,
                               vector<block_node>& v_block)
{
    red_blocks = 0;
    singleton_blocks = 0;
    elements_in_irred_blocks = 0;

    rank_support_v<> bf_rank(&bf);
    select_support_mcl<> bf_select(&bf);
    v_block = vector<block_node>(bf_rank(bf.size()-1));

    size_type block_id = 0;
	while ( block_id < v_block.size() ){
		unsigned char c1;
		size_type lb = bf_select(block_id+1), rb = bf_select(block_id+2);

		size_type rank1 = csa.wavelet_tree.inverse_select(lb, c1);
		size_type rank2 = csa.rank_bwt(rb, c1);
		if (rank2-rank1 == rb-lb) {
			if (rb-lb == 1) {
				++singleton_blocks;
				v_block[block_id] = block_node(0, 0, block_id);
			} else {
				++red_blocks;
				// determine who many times we can reduce the block until
				// we reach an irreducible block
				stack<size_type> cur_block_id;
				stack<size_type> cur_delta_x;
				size_type dest_block_id = block_id;
				// follow LF until we reach an already handled block or an
				// irreducible block
				do {
					cur_block_id.push(dest_block_id);
					dest_block_id = bf_rank(csa.psi(lb)+1)-1;
					size_type delta_x = csa.psi(lb) - bf_select(dest_block_id+1);
					cur_delta_x.push(delta_x);

					lb = bf_select(dest_block_id+1);
					rb = bf_select(dest_block_id+2);
					rank1 = csa.wavelet_tree.inverse_select(lb, c1);
					rank2 = csa.rank_bwt(rb, c1);
				} //   v_block[dest_block] not already defined and block is reducible
				while (v_block[dest_block_id].dest_block == 0 and rank2-rank1 == rb-lb);

				if (rank2 - rank1 != rb-lb) {  // if
					v_block[dest_block_id].dest_block = dest_block_id;
				}
				size_type delta_d = v_block[dest_block_id].delta_d;
				size_type delta_x = v_block[dest_block_id].delta_x;
				dest_block_id   = v_block[dest_block_id].dest_block;
				while (!cur_block_id.empty()) {
					size_type x = cur_block_id.top(); cur_block_id.pop();
					++delta_d;
					delta_x += cur_delta_x.top(); cur_delta_x.pop();
					v_block[ x ] = block_node(delta_x, delta_d, dest_block_id);
				}
			}
		} else { // irreducible block
			elements_in_irred_blocks += (rb-lb);
			v_block[block_id] = block_node(0, 0, block_id);
		}
		++block_id;
	}
}

//! Outputs the maximal delta_x and delta_d and the average delta_x and delta_d
/*
 * \param v_block A vector of block_node items
 */
void output_block_info_statistics(const vector<block_node>& v_block);

//! Close a stream if it is already opened
void close_stream_if_open(std::ifstream& in);

//! Returns the output directory
/*
 *	\param file_name	File name of the original text file
 *	\param output_dir	Char pointer to the output directory
 *	\return If output_dir==NULL the dirname(file_name) will be returned and
 *          output_dir otherwise.
 */
std::string get_output_dir(const char* file_name, const char* output_dir=NULL);


#endif
