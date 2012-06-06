/* This file contains a bunch of struct which are used to construct and analyse
 * the rosa.
 */
#ifndef ROSA_HELPING_STRUCTURES
#define ROSA_HELPING_STRUCTURES

#include <iostream>
#include <sdsl/int_vector.hpp>

typedef sdsl::bit_vector::size_type size_type;

enum  fringe_type	{ NO_FRINGE=0, RIGHT_FRINGE=1, LEFT_FRINGE=2, FULL_FRINGE=3 };

//! block_info holds information about the in-memory part of a rosa block
struct block_info {
    size_type bwd_lb; // Left bound of the block in the backward index.
    size_type depth;  // Depth of the block.
    size_type fwd_lb; // Left bound of the block in the forward index.
    size_type size;   // Size of the block.
    block_info(size_type bwd_lb, size_type depth, size_type fwd_lb, size_type size);
    bool operator<(const struct block_info& b)const; // sort according to bwd_lb, depth, fwd_lb and size
};

//! block_node holds information about the external-memory part of a rosa block
struct block_node {
    size_type delta_x; // Start position in in the corresponding irreducible block.
    size_type delta_d; // Number of iterations to reach the corresponding irreducible block.
    size_type dest_block; // Pointer to the irreducible block.
    size_type bwd_id; // id of the corresponding block in the backward index.

    block_node(size_type delta_x=0, size_type delta_d=0, size_type dest_block=0, size_type bwd_id=0);
};

//! bwd_block_info holds information about the  Stores information of blocks in the backward index
struct bwd_block_info {
    size_type bwd_lb; // Left bound of the block in the backward index.
    size_type bwd_rb; // Right bound of the block in the backward index
    size_type depth;  // Depth of the block.
    size_type block_id; // Id in the forward index.
    bwd_block_info(size_type bwd_lb=0, size_type bwd_rb=0, size_type depth=0, size_type block_id=0);
    bool operator<(const struct bwd_block_info& b)const;
};

std::ostream& operator<<(std::ostream& out, const bwd_block_info& x);

//! header_item holds the information for on header entry of an external block
struct header_item {
    size_type bwd_id;	// Backward index of the internal block
    size_type delta_x;	// delta_x of the internal block
    size_type delta_d;	// delta_d of the internal block
    header_item(size_type bwd_id=0, size_type delta_x=0, size_type delta_d=0);
    bool operator<(const struct header_item& it)const;
};


std::ostream& operator<<(std::ostream& out, const header_item& x);

#endif
