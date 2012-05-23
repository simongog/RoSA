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


