#ifndef FULLYMIXED_H
#define FULLYMIXED_H 
#include"network.h"

class CFullymixed : public CNetwork{
	public:
	CFullymixed(int pop);
 	private:
};

CFullymixed::CFullymixed(int pop){

	for(int i=0; i<pop; i++){
		CNode *p_node=new CNode(i);
		nodes.push_back(p_node);
		p_node->degree=pop-1;
	}
}

#endif /* FULLYMIXED_H */
