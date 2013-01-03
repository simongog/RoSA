#include "rosa_helping_functions.hpp"

bool check_size(size_type size){
	const uint32_t max_size = 100;
	if ( size > max_size ){
		std::cerr << "ERROR: The tikz_output size is limited to examples of size " 
			 << max_size << ". We are sorry for that." << std::endl; 
		return false;
	}
	return true;
}

void output_block_info_statistics(const vector<block_node> &v_block){
#ifdef OUTPUT_STATS
	size_type max_delta_x = 0;
	size_type max_delta_d = 0;
	long double avg_delta_x = 0;
	long double avg_delta_d = 0;
	for (size_t i=0; i < v_block.size(); ++i){
		if ( v_block[i].delta_x > max_delta_x ){
			max_delta_x = v_block[i].delta_x;
		}	
		avg_delta_x += v_block[i].delta_x;
		if ( v_block[i].delta_d > max_delta_d ){
			max_delta_d = v_block[i].delta_d;
		}
		avg_delta_d += v_block[i].delta_d;
	}
	avg_delta_x /= v_block.size();
	avg_delta_d /= v_block.size();
	std::cout << "# max_delta_x = " << max_delta_x << std::endl;
	std::cout << "# max_delta_d = " << max_delta_d << std::endl;
	std::cout << "# avg_delta_x = " << avg_delta_x << std::endl;
	std::cout << "# avg_delta_d = " << avg_delta_d << std::endl;
#endif	
}

void close_stream_if_open(std::ifstream &in){
	if ( in.is_open() ){
#ifdef DEBUG_DISK_ACCESS		
		std::cout << "close stream" << std::endl;
#endif			
		in.close();
	}
}


std::string get_output_dir(const char *file_name, const char *output_dir){
	if ( NULL==output_dir ){
		return sdsl::util::dirname(file_name);
	}
	return std::string(output_dir);
}

void mark_interval(bu_interval *v, bit_vector& bf, size_type b, size_type& blocks){
	if ( v->size() > b ){
		for (size_t i=0; i< v->children.size(); ++i){
			if (v->children[i]->size() <= b ){
				size_type lb = v->children[i]->lb;
				if ( 0 == bf[lb] ){
					bf[lb] = 1;
					++blocks;
				}
			}
		}
	}
	v->delete_children();	
}


void mark_blocks(const char* lcp_file, bit_vector& bf, size_type b, size_type& blocks){
	blocks = 0;
	typedef int_vector<>::size_type size_type;
	int_vector_file_buffer<> lcp_buf(lcp_file);
	bu_interval root(0,0,0);
	bu_interval *last_interval = NULL;
	std::stack<bu_interval*> stk;
	stk.push(&root);
	for (size_type i=0,r=0,r_sum=0; i < lcp_buf.int_vector_size;) { 
		for (; i < r+r_sum; ++i) {
			if (i > 0 ){
				uint64_t lcp = lcp_buf[i-r_sum];
				size_type lb = i-1;
				while ( lcp < stk.top()->lcp ){
					stk.top()->rb = i-1;
					last_interval = stk.top(); stk.pop();
					// process node
					mark_interval(last_interval, bf, b, blocks);
					lb = last_interval->lb;
					if ( lcp <= stk.top()->lcp ){
						stk.top()->children.push_back(last_interval);
						last_interval = NULL;
					}
				}
				if ( lcp > stk.top()->lcp ){
					bu_interval *v = new bu_interval(lb, 0, lcp);
					if ( last_interval != NULL ){
						v->children.push_back(last_interval);
						last_interval = NULL;
					}
					stk.push(v);
				}
			}
			bu_interval *leaf = new bu_interval(i, i, lcp_buf.int_vector_size);
			stk.top()->children.push_back(leaf);
		}
		r_sum += r; r = lcp_buf.load_next_block();
	}
	while ( !stk.empty() ){
		stk.top()->rb = lcp_buf.int_vector_size-1;
		last_interval = stk.top(); stk.pop();
		// process node
		mark_interval(last_interval, bf, b, blocks);
	}
}

