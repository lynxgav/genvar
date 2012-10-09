#include<iostream>
#include"strain.h"
#include"node.h"
#include"rrgraph.h"
#include"lattice.h"
#include"fullymixed.h"
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
double rec_rate = 0.05;
double r0;
int pop=900;
double avconnect=4;//pop-1;//4.;
double mig_rate;// = r0*rec_rate/(double)(Nd*avconnect);
double mutat_rate = 0.0004;
int nt=50;
int nd=5;
int t=1;
int tmax=1000000;
int tstep=1;
int tprint=100;
//time
bool fullymixed=false;

// change avconnect, CNetwork and poisson if fullymixed=true;

//CNetwork *contacts=new CRRGraph(5000, 4);
CNetwork *contacts=new CLattice(32, 32);
//CNetwork *contacts=new CFullymixed(pop);
CModel model(contacts);

unsigned int visited=0;
unsigned int iteration=0;

int distance(CStrain *s1, CStrain *s2){

	if(s1==s2) {return 0;}
	if(s1->gen > s2->gen) swap(s1,s2);
	if(s2->gen > s1->gen) return 1+distance(s1, s2->father());
	if(s2->gen == s1->gen) return 2+distance(s1->father(), s2->father());

	
	return -1;
}

CStrain* CommonFatherForTwoNodes(CStrain *s1, CStrain *s2){

	if(s1==s2) return s1;
	if(s1->gen > s2->gen) return CommonFatherForTwoNodes(s1->father(), s2);
	if(s2->gen > s1->gen) return CommonFatherForTwoNodes(s1, s2->father());
	if(s2->gen == s1->gen) return CommonFatherForTwoNodes(s1->father(), s2->father());

	return NULL;
}

CStrain* CommonFatherForSample(vector<CStrain*> &sample){
	if(sample.size()==0) {return NULL;}
	CStrain *commonfather=sample.at(0);

	for(unsigned int i=1;i<sample.size();i++){
		commonfather=CommonFatherForTwoNodes(sample.at(i),commonfather);
		assert(commonfather!=NULL);
	}

	assert(commonfather!=NULL);
	
	return commonfather;
}

/*
int seg_distance(CStrain *s1, CStrain *s2, CStrain *&f){
	if(s1->visited==visited and s2->visited==visited) {f=CommonFatherForTwoNodes(s1,s2);return 0;}
	if(s1==s2) {f=s1; 
		s1->visited=s2->visited=visited;
		return 0;
		}
	if(s1->gen > s2->gen)swap(s1,s2);

	if(s2->gen > s1->gen){
		if(s2->father()->visited==visited){
			s1->visited=s2->visited=visited; 
			return 1;
		}
		s1->visited=s2->visited=visited; 
		return 1+seg_distance(s1, s2->father(),f);
	}
	s1->visited=s2->visited=visited; 
	if(s2->gen == s1->gen) return 2+seg_distance(s1->father(), s2->father(),f);

	f=NULL; 
	return -1;
}
*/

int seg_distanceShahbanu(CStrain  *s, CStrain *gfather){

	if(s->visitedShahbanu==iteration or s==gfather) return 0;

	s->visitedShahbanu=iteration;
	return 1+seg_distanceShahbanu(s->father(),gfather);
}

int NSegSitesShahbanu(vector<CStrain*> &sample){
	iteration++;

	if(sample.size()==0)return 0;
	CStrain *cfather=CommonFatherForSample(sample);
	//cerr << t << "  " << cfather->gen << "   ";
	assert(cfather!=NULL);
	int d=0;
	for(unsigned int i=0;i<sample.size();i++){
		d+=seg_distanceShahbanu(sample.at(i), cfather);	
	}

	return d;
}

//CStrain* cf=NULL;

/*
int NSegSites(vector<CStrain*> &sample){
	//cerr<< t << "  visited before increment in function " << visited << endl; 
	visited++;
	//cerr<< t << "  visited after increment in function " << visited << endl;
	if(sample.size()==0)return 0;
	CStrain *common_father=sample.at(0);
	int d=0;
	for(unsigned int i=1;i<sample.size();i++){
		
		d+=seg_distance(common_father, sample.at(i), common_father);	
		assert(common_father!=NULL);
	}

	cf=common_father;
	
	return d;
}
*/

class CDiversityOut
	{
	public:
	CDiversityOut(double d=0, int n=0):distance(d),SegSites(n){}
	double distance;
	int SegSites;
 	private:
	};

CDiversityOut GeneticDiversityGlobal(){

	vector<CStrain*> sample;

	std::tr1::uniform_int<int> unif4(0,model.system_state.at(INF).size()-1);
	std::tr1::uniform_int<int> unif5(0,Nd-1);

	for (int i=0; i<nt; i++){
		int chosen=unif4(eng);
		CNode* p_node=model.system_state.at(INF).at(chosen);
		/*if(t==1458)//t==2915 || t==2467){
			cerr << p_node->ID << "  ";
		}*/
		chosen=unif5(eng);
		sample.push_back(p_node->pathogens.at(chosen));
	}
	//cerr<< t << "  visited before function executed " << visited << endl;

//	int n1=NSegSites(sample);
	int n2=NSegSitesShahbanu(sample);

	//if(n1!=n2){
		//cerr<< t <<"   n seg sites= "<< n1-n2 <<"  n strains= "<<sample.size()<< "  allstrains  " << allstrains.size() << " visited " << visited << "   f->gen" << cf->gen <<endl;
	//}

	//cerr << n2 << endl;

	//CStrain *f2=CommonFatherForSample(sample);

	//cerr<< cf->gen <<"\t"<< f2->gen <<endl;
	
	//assert(cf==f2);
	
	/*
	if(t==7){//t==2915 || t==2467)
		for (int i=0; i<nt; i++){	
			cerr << sample.at(i)->visited << "  ";
		}
	}
	*/	

	double dist=0.;
	int k=0;

	for(int i=0; i<nt; i++){
		for(int j=i+1; j<nt; j++){
			k=distance( sample.at(i), sample.at(j) );
			assert(k!=-1);
			dist+=k;
		}
	}

	dist=dist/(nt*(nt-1.)/2.);

	return CDiversityOut(dist, n2);
}

double GeneticDiversityLocal(CNode* p_node){

	map<int,CStrain*> sample;

	std::tr1::uniform_int<int> unif6(0,Nd-1);

	assert(nd<=Nd);

	while (sample.size()<(unsigned int) nd){
		int chosen=unif6(eng);
		sample[chosen]=p_node->pathogens.at(chosen);
	}
	
	double dist=0.;
	int k=0;

	map<int,CStrain*>::iterator it1, it2;

	for(it1=sample.begin(); it1!=sample.end(); it1++){
		for(it2=it1,it2++; it2!=sample.end(); it2++){
			k=distance( it1->second, it2->second );
			assert(k!=-1);
			dist+=k;
		}
	}

	dist=dist/(nd*(nd-1.)/2.);

	return dist;
}

double GeneticDiversityLocalAverage(){

	double dist=0.;

	for (unsigned int i=0; i < model.system_state.at(INF).size(); i++){
		dist+=GeneticDiversityLocal( model.system_state.at(INF).at(i) );
	}

	dist=dist/(double)model.system_state.at(INF).size();

	return dist;
}

void InitialConditions(){

	mig_rate = r0*rec_rate/(double)(Nd*avconnect);

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

		// my version
		//std::tr1::poisson_distribution<double> poisson( Nd*mig_rate*p_node->degree );
		// Diana's version
		std::tr1::poisson_distribution<double> poisson( Nd*mig_rate*p_node->count_neighbours_state(SUS) );
		//fully mixed Diana's version
		//std::tr1::poisson_distribution<double> poisson( Nd*mig_rate*model.system_state.at(SUS).size() );
		//fully mixed my version
		//std::tr1::poisson_distribution<double> poisson( Nd*mig_rate*(pop-1) );

		unsigned int nmigrants = (unsigned int)poisson(eng);
		// check if it makes sense
		if(nmigrants > (unsigned int) Nd){
			//cerr<<nmigrants<<endl;
			nmigrants=(unsigned int) Nd;

		}
		assert( nmigrants <= (unsigned int) Nd);
		assert( p_node->pathogens.size()==(unsigned int) Nd);
		
		unsigned int np=p_node->pathogens.size();
		std::tr1::uniform_int<int> unif1(0,np-1);
		while(p_node->migrants.size()<nmigrants){
			int chosen=unif1(eng);
			p_node->migrants[chosen]=p_node->pathogens.at(chosen);
		}
	}

	
	for (unsigned int i=0; i < model.system_state.at(INF).size(); i++){
		CNode* p_node=model.system_state.at(INF).at(i);

		unsigned int nrecipients;
		if(!fullymixed){
			nrecipients=p_node->count_neighbours_state(SUS);
		}
		else {
			nrecipients=model.system_state.at(SUS).size();
		}

		if(nrecipients==0) {
			p_node->migrants.clear();
			continue;
		}

		std::tr1::uniform_int<int> unif2(1,nrecipients);

		map<int,CStrain*>::iterator itt;

		for(itt=p_node->migrants.begin(); itt != p_node->migrants.end(); itt++){

			int chosen=unif2(eng);

			if(!fullymixed){
				list<CNode*>::iterator it; int count=0;
				for(it=p_node->neighbours.begin(); count!=chosen; it++) {
					if( (*it)->state == SUS ){
						count++;
						if(count==chosen){
							break;
						}
					}
				}
				(*it)->pathogens.push_back(itt->second);
			}
			else {
				model.system_state.at(SUS).at(chosen-1)->pathogens.push_back(itt->second);
			}
		}
		p_node->migrants.clear();

	}
	/*
	int j=0;

	vector<CNode *>::iterator it;
	for (it=model.network->nodes.begin(); it!=model.network->nodes.end(); it++) {
		if( (*it)->pathogens.size()>0 && (*it)->state==SUS ) {
			sus--; inf++;
			model.UpdateSystemState( (*it), SUS, INF);
			j++;
			
		}
		assert( (*it)->migrants.size()==0 );
		if ( (*it)->pathogens.size()>0) { assert ( (*it)->state==INF ); }
	}
	*/

	int j=0;
	vector<CNode*> susceptibles;
	susceptibles=model.system_state.at(SUS);
	for (unsigned int i=0; i < susceptibles.size(); i++){
		CNode* p_node=susceptibles.at(i);
		if( p_node->pathogens.size()>0 ){
			sus--; inf++;
			model.UpdateSystemState( p_node, SUS, INF);
			j++;
		}
	}
	//cerr << "number of migrations  " << j << "  Nd  " << Nd << "  mig_rate  " << mig_rate << "  sus  " << model.system_state.at(SUS).size() << endl;

	vector<CNode *>::iterator it;
	for (it=model.network->nodes.begin(); it!=model.network->nodes.end(); it++) {
		assert( (*it)->migrants.size()==0 );
		if ( (*it)->pathogens.size()>0) { assert ( (*it)->state==INF ); }
	}

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
		if(t%tprint==0){
			//cerr<<"n seg sites= "<<NSegSites(allstrains)<<"  n strains= "<<allstrains.size()<<endl;
			double GDLAverage=GeneticDiversityLocalAverage();
			CDiversityOut GDGlobal=GeneticDiversityGlobal();
			cout<< t <<"\t"<< sus/(double)model.network->get_N() <<"\t"<< inf/(double)model.network->get_N() <<"\t"<< rec/(double)model.network->get_N() <<"\t"<< GDGlobal.distance <<"\t"<< GDLAverage  <<"\t"<< GDGlobal.SegSites<<"\t"<< allstrains.size()<< endl;
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
	if(argc>2)r0=atof(argv[2]);
	eng.seed(seed);
	//cerr<<seed<<endl;

	model.Initial_Conditions();
	InitialConditions();

	Iterate();

return 0;

}
