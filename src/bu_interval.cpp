#include "bu_interval.hpp"


bu_interval::bu_interval(size_type _lb, size_type _rb, size_type _lcp) : lb(_lb), rb(_rb), lcp(_lcp){};

bu_interval::~bu_interval(){ delete_children(); }

void bu_interval::delete_children(){
		for (size_t i=0; i<children.size(); ++i){ //std::cout<<"delete "<<children[i]->lcp<<"-["<<children[i]->lb<<","<<children[i]->rb<<"]"<<std::endl;
			delete children[i]; 
		}
		children.resize(0);
}

bu_interval::size_type bu_interval::size()const{
	return rb-lb+1;
}

