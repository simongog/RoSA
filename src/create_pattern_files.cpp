#include <vector>
#include <utility>
#include <string>
#include <stack>
#include <algorithm>
#include "bu_interval.hpp"
#include "sdsl/int_vector.hpp"
#include <sdsl/util.hpp>

using namespace std;
using namespace sdsl;

typedef vector<size_t> tVI;

void get_candidate(bu_interval* v, const tVI &pattern_len, 
		                           const tVI &min_occ, 
								   const tVI &max_occ, 
								   vector<tVI> &candidates,
								   tVI &candidates_nr){
	
	for (size_t p=0; p<pattern_len.size(); ++p){
		if ( v->lcp < pattern_len[p] ){ // if depth of node is smaller than the desired pattern_len
			for (size_t i=0; i<v->children.size(); ++i){   // iterate over all children
				if ( v->children[i]->lcp >= pattern_len[p] ){ // and check, if the depth is greater or equal than pattern_len
					size_t occ = v->children[i]->size();
					if ( min_occ[p] <= occ and occ <= max_occ[p] ){ // if frequency constraints are meet
						++candidates_nr[p];
						if ( candidates_nr[p] < 1000000 ){ // limit candidates to ten million
							candidates[p].push_back( v->children[i]->lb ); // add candidate
						}
					}
				}
			}
		}
	}
	
	v->delete_children();
}

int main(int argc, char *argv[]){
	string lcp_file, sa_file, text_file, basename;
	size_t pattern_cnt=1000;
	if ( argc < 4  ){
		cout << "Usage: " << argv[0]<<" lcp_file sa_file text_file [pattern_cnt]" << endl;
		cout <<" output 25 files basename(text_file).<pattern_len>.<pattern_cnt>.0.<min_occ>.<max_occ>.pattern"<<endl;
		return 1;
	}
	lcp_file = argv[1];
	cout<<"lcp_file = "<<lcp_file<<endl;
	sa_file  = argv[2];
	cout<<"sa_file = "<<sa_file<<endl;
	text_file= argv[3]; 
	cout<<"text_file = "<<text_file<<endl;
	basename = util::basename(argv[3]);
	cout<<"basename = "<<basename<<endl;

	tVI len;
	
	len.push_back(4);
	len.push_back(10);
	len.push_back(20);
	len.push_back(50);
	len.push_back(100);
	tVI occ1, occ2;
	occ1.push_back(1); occ2.push_back(1);
	occ1.push_back(8); occ2.push_back(12);
	occ1.push_back(75); occ2.push_back(125);
	occ1.push_back(750); occ2.push_back(1250);
	occ1.push_back(7500); occ2.push_back(12500);
	tVI pattern_len, min_occ, max_occ;
	for(size_t i=0; i<len.size(); ++i){
		for(size_t j=0; j<occ1.size(); ++j){
			pattern_len.push_back(len[i]);
			min_occ.push_back(occ1[j]);
			max_occ.push_back(occ2[j]);
		}
	}
	int_vector_file_buffer<> lcp_buf(lcp_file.c_str());
	bu_interval root(0,0,0);
	bu_interval *last_interval = NULL;
	stack<bu_interval*> stk;
	vector<tVI> candidates(pattern_len.size(), tVI());
	tVI candidates_nr(pattern_len.size());
	stk.push(&root);
	for (size_t i=0,r=0,r_sum=0; i < lcp_buf.int_vector_size;) { 
		for (; i < r+r_sum; ++i) {
			if (i > 0 ){
				uint64_t lcp = lcp_buf[i-r_sum];
				size_t lb = i-1;
				while ( lcp < stk.top()->lcp ){
					stk.top()->rb = i-1;
					last_interval = stk.top(); stk.pop();
					// process node
					get_candidate(last_interval, pattern_len, min_occ, max_occ, candidates, candidates_nr);
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
			bu_interval *leaf = new bu_interval(i,i, lcp_buf.int_vector_size);
			stk.top()->children.push_back(leaf);
		}
		r_sum += r; r = lcp_buf.load_next_block();
	}
	while ( !stk.empty() ){
		stk.top()->rb = lcp_buf.int_vector_size-1;
		last_interval = stk.top(); stk.pop();
		// process node
		get_candidate(last_interval, pattern_len, min_occ, max_occ, candidates, candidates_nr);
	}

	for(size_t rr=0; rr<pattern_len.size(); ++rr){
		tVI &cand = candidates[rr];
		bit_vector already_used(cand.size(), 0);
		string pattern_file = basename+"."+util::to_string(pattern_len[rr])+"."
											       +util::to_string(pattern_cnt)+".0."
												   +util::to_string(min_occ[rr])+"."
												   +util::to_string(max_occ[rr])+".pattern";
		std::ofstream pattern_stream( pattern_file.c_str() );
		cout<<"pattern_file="<<pattern_file<<endl;
		if ( pattern_stream ){
			cout<<"# pattern_len = " <<pattern_len[rr]<<endl;
			cout<<"# min_occ = "<<min_occ[rr]<<endl;
			cout<<"# max_occ = "<<max_occ[rr]<<endl;
			cout<<"# candidatess = "<<candidates_nr[rr]<<endl;
			util::write_member(pattern_cnt, pattern_stream); 
			util::write_member(pattern_len[rr], pattern_stream); 
			uint64_t swaps=0;
			util::write_member(swaps, pattern_stream); 
			if ( cand.size() > 0 ) {

				ifstream text_stream(text_file.c_str());	
				if ( text_stream ){
					std::vector<size_t> sa_indexes;
					for (size_t i=0, j, k; i < pattern_cnt; ++i){
						j=rand()%cand.size();
						if ( already_used[j] ) {
							// find next free location 
							size_t j1=j+1;
							while ( j1 < already_used.size() and already_used[j1] ){
								++j1;
							}
							if ( j1 >= already_used.size() ){
								j1 = j;
							}
							j = j1;
						}
						sa_indexes.push_back(cand[j]);

						already_used[j] = 1;
					}

					int_vector_file_buffer<> sa_buf(sa_file.c_str());
					size_t idx_sa_indexes = 0;
					std::cerr<<"--sort sa_indexes--"<<std::endl;
					std::sort(sa_indexes.begin(), sa_indexes.end());
					std::cerr<<"--stream sa--"<<std::endl;
					char *pattern_buf = new char[pattern_len[rr]+1];
					pattern_buf[pattern_len[rr]] = '\0';
					for (size_t i=0,r=0,r_sum=0; i < sa_buf.int_vector_size and idx_sa_indexes < sa_indexes.size();) { 
						for (; i < r+r_sum and idx_sa_indexes < sa_indexes.size(); ++i) {
							if ( i == sa_indexes[idx_sa_indexes] ){
	//							std::cout<<"idx_sa_indexes="<<idx_sa_indexes<<" sa_indexes["<<idx_sa_indexes<<"]="<<sa_indexes[idx_sa_indexes]<<std::endl;
								++idx_sa_indexes;
	//							std::cout<<idx_sa_indexes<<" sa_indexes.size()"<<sa_indexes.size()<<endl;
								size_t sa_val = sa_buf[i-r_sum];
	//							std::cout<<"sa_val="<<sa_val<<std::endl;
								// extract the pattern
								text_stream.seekg(sa_val ,std::ios::beg);
								text_stream.read(pattern_buf, pattern_len[rr]);
	//							cout<<pattern_buf<<endl;
								pattern_stream.write( pattern_buf, pattern_len[rr] );
							}
						}
						if ( idx_sa_indexes < sa_indexes.size() ){
							r_sum += r; r = sa_buf.load_next_block();
						}
					}
					cout<<"--end stream sa"<<endl;
					delete [] pattern_buf;
					cout<<"--deleted pattern_buf"<<endl;
					text_stream.close();
					cout<<"text_stream.closed"<<endl;
				}else{
					std::cerr << "ERROR: Could not open text file " << text_file << endl;	
				}
			}
			pattern_stream.close();
			cout<<"pattern_stream.closed"<<endl;
		}else{
			std::cerr << "ERROR: Could not open pattern file " << pattern_file << endl;
		}
	}



}


