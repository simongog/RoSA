#include "pattern_file.hpp"
#include <sdsl/testutils.hpp>
#include <stack>

pattern_file::pattern_file(const char* pattern_file_name):pattern_cnt(m_pattern_cnt), pattern_len(m_pattern_len), swaps(m_swaps) {
    m_buf = NULL;
    m_pattern_file_name = string(pattern_file_name);
}

void pattern_file::reset() {
    if (m_pattern_stream.is_open()) {
        m_pattern_stream.close();
    }
    if (m_buf != NULL)
        delete [] m_buf;
    m_pattern_stream.open(m_pattern_file_name.c_str());
    if (m_pattern_stream.is_open()) {
        util::read_member(m_pattern_cnt, m_pattern_stream);
        util::read_member(m_pattern_len, m_pattern_stream);
        util::read_member(m_swaps, m_pattern_stream);
        m_buf = new char[m_pattern_len];
    } else {
        std::cerr << "ERROR: Could not open pattern file " << m_pattern_file_name << endl;
    }
}

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
void pattern_file::generate(const char* text_file_name, size_type pattern_cnt, size_type pattern_len, size_type swaps) {
    m_pattern_cnt = pattern_cnt;
    m_pattern_len = pattern_len;
    m_swaps		  = swaps;
    std::ofstream pattern_stream(m_pattern_file_name.c_str());
    if (pattern_stream) {
        util::write_member(m_pattern_cnt, pattern_stream);
        util::write_member(m_pattern_len, pattern_stream);
        util::write_member(m_swaps, pattern_stream);
        char* pattern_buf = new char[pattern_len+1];
        pattern_buf[pattern_len] = '\0';
        size_type n = 0;
        ifstream text_stream(text_file_name);
        n = util::get_file_size(text_file_name);
        if (text_stream) {
            for (size_type i=0, j, k; i < pattern_cnt; ++i) {
                j=rand()%n;
                if (j+pattern_len > n)
                    j = 0;
                // extract the pattern
                text_stream.seekg(j ,std::ios::beg);
                text_stream.read(pattern_buf, k = std::min(pattern_len, n-j));
                for (; k < pattern_len; ++k)
                    pattern_buf[k] = (char)255;
                // apply the swaps
                for (size_type s=0; s < swaps and pattern_len > 1; ++s) {
                    size_type r = rand()%(pattern_len-1);
                    std::swap(pattern_buf[r], pattern_buf[r+1]);
                }
// cout << i << " " << pattern_buf <<endl;
                pattern_stream.write(pattern_buf, pattern_len);
            }
            text_stream.close();
        } else {
            std::cerr << "ERROR: Could not open text file " << text_file_name << endl;
        }
        pattern_stream.close();
    } else {
        std::cerr << "ERROR: Could not open pattern file " << m_pattern_file_name << endl;
    }
}

pattern_file::~pattern_file() {
    if (m_buf != NULL) {
        delete [] m_buf;
    }
    if (m_pattern_stream.is_open()) {
        m_pattern_stream.close();
    }
}

const char* pattern_file::get_next_pattern() {
    m_pattern_stream.read(m_buf, m_pattern_len);
    return m_buf;
}

void pattern_file::remove() {
    std::remove(m_pattern_file_name.c_str());
}

void pattern_file::generate_restricted(const char* lcp_file, const char *sa_file, const char *text_file_name, size_type pattern_cnt, 
							 size_type pattern_len, size_type min_occ, size_type max_occ){
	typedef int_vector<>::size_type size_type;
	int_vector_file_buffer<> lcp_buf(lcp_file);
	bu_interval root(0,0,0);
	bu_interval *last_interval = NULL;
	m_pattern_cnt = pattern_cnt;
	m_pattern_len = pattern_len;
	std::stack<bu_interval*> stk;
	std::vector<size_t> candidates;
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
					get_candidate(last_interval, pattern_len, min_occ, max_occ, candidates);
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
		get_candidate(last_interval, pattern_len, min_occ, max_occ, candidates);
	}
	
}

void pattern_file::get_candidate(bu_interval* v, size_t pattern_len, size_t min_occ, size_t max_occ, std::vector<size_t> &candidates){
	if ( v->lcp < pattern_len ){ // if depth of node is smaller than the desired pattern_len
		for (size_t i=0; i<v->children.size(); ++i){   // iterate over all children
			if ( v->children[i]->lcp >= pattern_len ){ // and check, if the depth is greater or equal than pattern_len
				size_type occ = v->children[i]->size();
				if ( min_occ <= occ and occ <= max_occ ){ // if frequency constraints are meet
					candidates.push_back( v->children[i]->lb ); // add candidate
				}
			}
		}
	}
}
