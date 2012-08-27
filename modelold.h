#ifndef MODEL_H
#define MODEL_H 
#include"network.h"

double beta=2.5;
double delta=1;
double gama=0.05;

enum{INF=0, REC=1, SUS=2};


class CModel
	{
	public:
	CModel(CNetwork *net, int ns){
		network=net;
		nstates=ns;
		}

	void Initial_Conditions();
	void Fill_Network();

	int Choose_Transition();
	void Execute_Transition(int);
	CNode *Transition(int from, int to);
	void Update_State();
	
	int Count_State(int start, int end);

	CNetwork *network;
	double elapsed;
	int nstates;
	vector< vector<CNode*> > system_state;
 	private:
	};

void CModel::Initial_Conditions(){

	elapsed=0;
	int inf=20748;
	int sus=498038;
	int rec=network->get_N()-inf-sus;
	ERROR(rec<0, "The sum of infected and susceptibles is larger than total.");

	std::tr1::uniform_real<double> unif(0, 1);

	double s, i, rand;
	int j=network->get_N();
	list<CNode *>::iterator it;
	//Distribute the state randomly
	for (it=network->nodes.begin(); it!=network->nodes.end(); it++) {
		s=(double)sus/(double)j; i=(double)inf/(double)j;
		rand=unif(eng);
		if (rand < s ) {(*it)->state=SUS; sus--;}
		else if (rand <s+i) {(*it)->state=INF; inf--;}
		else {(*it)->state=REC; rec--;}
		--j;
		}
	
	for(int i=0; i<nstates; i++){
		system_state.push_back(vector<CNode*>() );
		}	
	for (it=network->nodes.begin(); it!=network->nodes.end(); it++) {
		if( (*it)->state!=SUS)system_state.at((*it)->state).push_back((*it));
		else{
			int nc=(*it)->count_neighbours_state(INF);
			system_state.at(nc+SUS).push_back((*it));
			}
		}

	for(int i=0; i<nstates; i++){
		cerr<<system_state.at(i).size()<<endl;
		}
}

int CModel::Count_State(int start, int end){
	int count=0;
	for(int i=start; i<=end; i++){
		count+=system_state.at(i).size();
	}
	return count;
}

int CModel::Choose_Transition(){
	double prob[nstates];
	double p_sum=0;
	prob[0]=p_sum+=delta*system_state.at(INF).size();
	prob[1]=p_sum+=gama*system_state.at(REC).size();

	prob[2]=0;
	for(int i=3; i<nstates; i++){
		prob[i]=p_sum+=beta*system_state.at(i).size()*(i-SUS);
	}
	
	std::tr1::uniform_real<double> unif(0, 1);
	elapsed-=log(1.0-unif(eng))/p_sum;
	p_sum*=unif(eng);
	int k;
	for(k=0; p_sum>=prob[k]; k++);
	return k;
}

void UpdateNode(CNode* node){

	system_state.at(from).at(pos)=system_state.at(from).back();
	system_state.at(from).pop_back();
	node->state=to;
	system_state.at(to).push_back(node);
}

CNode *CModel::Transition(int from, int to){
	assert(from!=to);

	int max=system_state.at(from).size();
	std::tr1::uniform_int<long int> unif(0, max-1);
	int pos=unif(eng);
	CNode* node=system_state.at(from).at(pos);

	if(from==INF)node.UpdateNeighbourState(-1);	
	if(to==INF)node.UpdateNeighbourState(1);	
	return node;
}

void CModel::Execute_Transition(int c){
}


void CModel::Update_State(){

}

#endif /* MODEL_H */
