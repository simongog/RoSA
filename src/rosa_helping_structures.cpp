#include "rosa_helping_structures.hpp"

block_info::block_info(size_type bwd_lb, size_type depth, size_type fwd_lb, size_type size):
				       bwd_lb(bwd_lb), depth(depth), fwd_lb(fwd_lb), size(size){}

bool block_info::operator<(const struct block_info &b)const{
	if( bwd_lb != b.bwd_lb )
		return bwd_lb < b.bwd_lb;
	else if(depth != b.depth)
		return depth < b.depth;
	else if(fwd_lb != b.fwd_lb)
		return fwd_lb < b.fwd_lb;
	else
		return size < b.size;
}

block_node::block_node(size_type delta_x, size_type delta_d, size_type dest_block,
			   size_type bwd_id): delta_x(delta_x), delta_d(delta_d), 
								  dest_block(dest_block), bwd_id(bwd_id) {}

bwd_block_info::bwd_block_info(size_type bwd_lb, size_type bwd_rb, size_type depth, size_type block_id):
		       bwd_lb(bwd_lb), bwd_rb(bwd_rb), depth(depth), block_id(block_id){}

bool bwd_block_info::operator<(const struct bwd_block_info &b)const{
	if ( bwd_lb != b.bwd_lb )
		return bwd_lb < b.bwd_lb;
	else if (depth != b.depth)
		return depth < b.depth;
	return block_id < b.block_id;
}

std::ostream& operator<<(std::ostream &out, const bwd_block_info &x){
	out << "{" << x.bwd_lb << "," << x.bwd_rb << "," << x.depth << ",";
	out << x.block_id << "}" << std::endl;
	return out;
}

header_item::header_item(size_type bwd_id, size_type delta_x, size_type delta_d):
				bwd_id(bwd_id),delta_x(delta_x),delta_d(delta_d){}

bool header_item::operator<(const struct header_item &it)const{
	if ( bwd_id != it.bwd_id )
		return bwd_id < it.bwd_id;
	else if ( delta_x != it.delta_x)
		return delta_x < it.delta_x;
	else
		return delta_d < it.delta_d;
}


std::ostream& operator<<(std::ostream &out, const header_item &x){
	out << "{" << x.bwd_id << "," << x.delta_x << "," << x.delta_d << "}" << std::endl;
	return out;
}
