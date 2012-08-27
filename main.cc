#include<iostream>
#include"strain.h"
#include"node.h"
#include"rrgraph.h"
#include"lattice.h"
#include"model.h"
using namespace std;
int seed=0;
std::tr1::ranlux64_base_01 eng(seed);
std::tr1::uniform_real<double> unif(0, 1);

//To to be used for assigning ID's
int stotal;
//List of all alive strains
vector<CStrain*> allstrains;
vector<CStrain*> strains;
CStrain *top=NULL;
int Nd = 5;
double rec_rate = 0.02;
double r0 = 2.;
double avconnect=4.;
double mig_rate = r0*rec_rate/(double)(Nd*avconnect);
double mutat_rate = 0.0001;
int nt=100;
int t=1, tmax=1000000, tstep=1;
//time

//CNetwork *contacts=new CRRGraph(5000, 4);
CNetwork *contacts=new CLattice(30, 30);
CModel model(contacts);

int distance(CStrain *s1, CStrain *s2){

	if (s1==s2) return 0;
	if(s1->gen > s2->gen) return 1+distance(s1->father(), s2);
	if(s2->gen > s1->gen) return 1+distance(s1, s2->father());
	if(s2->gen == s1->gen) return 2+distance(s1->father(), s2->father());

	return -1;

}

double GeneticDiversityGlobal(){

	vector<CStrain*> sample;

	std::tr1::uniform_int<int> unif4(0,model.system_state.at(INF).size()-1);
	std::tr1::uniform_int<int> unif5(0,Nd-1);

	for (int i=0; i<nt; i++){
		int chosen=unif4(eng);
		CNode* p_node=model.system_state.at(INF).at(chosen);
		chosen=unif5(eng);
		sample.push_back(p_node->pathogens.at(chosen));
	}

	double dist=0.;
	int k=0;

	for(int i=0; i<nt; i++){
		for(int j=i+1; j<nt; j++){
			k=distance( sample.at(i), sample.at(j) );
			assert(k!=-1);
			dist+=k;
		}
	}

	dist=dist/(nt*(nt-1)/2.);

	return dist;
}

void InitialConditions(){

	stotal=0;
	top = new CStrain(stotal,NULL);
	//top->color=1;
	allstrains.push_back(top);
	stotal++;

	for (unsigned int i=0; i < model.system_state.at(INF).size(); i++){
		for(int j=0; j < Nd; j++){
			model.system_state.at(INF).at(i)->pathogens.push_back(top);
		}
	}
}

void Recovery(){

	int j=0;

	vector<CNode*> infectives;

	infectives=model.system_state.at(INF);

	for (unsigned int i=0; i < infectives.size(); i++){
		if( unif(eng) < rec_rate ){
			inf--; sus++;
			j++;
			CNode* p_node=infectives.at(i);
			model.UpdateSystemState( p_node, INF, SUS);
			p_node->pathogens.clear();
		}
	}
	

	//cerr << "number of recoveries  " << j << endl;
}

void Migration(){

	for (unsigned int i=0; i < model.system_state.at(INF).size(); i++){
		CNode* p_node=model.system_state.at(INF).at(i);
		std::tr1::poisson_distribution<double> poisson( Nd*mig_rate*p_node->degree );
		unsigned int nmigrants = (unsigned int)poisson(eng);
		assert( nmigrants <= (unsigned int) Nd);//{
			//cerr<<"Number of migrants <= Nd"<<endl;
		//};
		unsigned int np=p_node->pathogens.size();
		std::tr1::uniform_int<int> unif1(0,np-1);
		while(p_node->migrants.size()<nmigrants){
			int chosen=unif1(eng);
			p_node->migrants.push_back(p_node->pathogens.at(chosen));
		}
	}

	
	for (unsigned int i=0; i < model.system_state.at(INF).size(); i++){
		CNode* p_node=model.system_state.at(INF).at(i);
		unsigned int nrecipients=p_node->count_neighbours_state(SUS);

		if(nrecipients==0) {
			p_node->migrants.clear(); 
			continue;
		}

		std::tr1::uniform_int<int> unif2(1,nrecipients);

		for(unsigned int j=0; j < p_node->migrants.size(); j++){
			list<CNode*>::iterator it; int count=0;
			int chosen=unif2(eng);
			for(it=p_node->neighbours.begin(); count!=chosen; it++) {
				if( (*it)->state == SUS ){
					count++;
					if(count==chosen){
						break;
					}
				}
			}
			(*it)->pathogens.push_back(p_node->migrants.at(j));
		}
		p_node->migrants.clear();

	}

	int j=0;

	vector<CNode *>::iterator it;
	for (it=model.network->nodes.begin(); it!=model.network->nodes.end(); it++) {
		if( (*it)->pathogens.size()>0 && (*it)->state==SUS ) {
			sus--; inf++;
			model.UpdateSystemState( (*it) , SUS, INF);
			j++;
			
		}
		assert( (*it)->migrants.size()==0 );
		if ( (*it)->pathogens.size()>0) { assert ( (*it)->state==INF ); }
	}

	//cerr << "number of migrations  " << j << endl;

}

void Reproduction(){

	int k=0;

	for (unsigned int i=0; i < model.system_state.at(INF).size(); i++){
		CNode* p_node=model.system_state.at(INF).at(i);
		vector<CStrain*> newgeneration;
		unsigned int npath=p_node->pathogens.size();
		std::tr1::uniform_int<int> unif3(0,npath-1);
		for(int j=0; j < Nd; j++){
			int chosen=unif3(eng);
			CStrain *s = p_node->pathogens.at(chosen);
			if( unif(eng)< mutat_rate){
				k++;
				CStrain *s2 = new CStrain(stotal,s);
				stotal++;
				allstrains.push_back(s2);
				s=s2;
			}
			newgeneration.push_back(s);
		}
		p_node->pathogens.clear();
		p_node->pathogens=newgeneration;

		assert(p_node->pathogens.size()==(unsigned int)Nd);
	}

	//cerr << "number of mutations " << k << endl;

}

void Update(){
	Recovery();
	Migration();
	Reproduction();
}

void Iterate(){
	while(t<=tmax and inf>0 and sus>=0 and rec>=0){
		if(t%100==0){
			cout<< t <<"\t"<< sus/(double)model.network->get_N() <<"\t"<< inf/(double)model.network->get_N() <<"\t"<< rec/(double)model.network->get_N() <<"\t"<< GeneticDiversityGlobal()<<endl;
		}
	
		Update();
		//cerr << network->average_degree() <<"\t"<< model.state_degree(SUS) <<"\t"<< model.state_degree(INF) <<"\t" << model.state_degree(REC) <<endl;
		t+=tstep;
		
		assert( (unsigned int)(sus+inf+rec)==model.network->get_N());
		assert((unsigned int)sus==model.system_state.at(SUS).size() && (unsigned int)inf==model.system_state.at(INF).size() && (unsigned int)rec==model.system_state.at(REC).size() );
	}

}


int main(int argc, char **argv){

	if(argc>1)seed=atoi(argv[1]);
	eng.seed(seed);

	model.Initial_Conditions();
	InitialConditions();
	Iterate();

return 0;

}
