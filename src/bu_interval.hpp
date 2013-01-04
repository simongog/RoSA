// A class for the intervals of a bottom-up traversal of the suffix tree.

#ifndef BU_INTERVAL
#define BU_INTERVAL

#include <vector>
#include <sdsl/int_vector.hpp>

// bottom up interval
struct bu_interval{
	typedef sdsl::int_vector<>::size_type size_type;
	size_type lcp, lb, rb; // lcp = depth of the interval, lb/rb = lef/right border of the interval
	std::vector<bu_interval*> children;
	bu_interval(size_type _lb, size_type _rb, size_type _lcp);
	~bu_interval(); 
	void delete_children();
	size_type size()const;
};

#endif
