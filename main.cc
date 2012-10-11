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
int pop=1000;
double avconnect=pop-1;//4.;
double mig_rate;// = r0*rec_rate/(double)(Nd*avconnect);
double mutat_rate = 0.0004;
int nt=50;
int nd=5;
int t=1;
int tmax=1000000;
int tstep=1;
int tprint=100;
//time
bool fullymixed=true;
bool versionD=false;
bool withreplacement=true;

// change avconnect, CNetwork and poisson if fullymixed=true;

//CNetwork *contacts=new CRRGraph(5000, 4);
//CNetwork *contacts=new CLattice(32, 32);
CNetwork *contacts=new CFullymixed(pop);
CModel model(contacts);

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

int seg_distance(CStrain  *s, CStrain *gfather){

	if(s->visited==iteration or s==gfather) return 0;

	s->visited=iteration;
	return 1+seg_distance(s->father(),gfather);
}

int NSegSites(vector<CStrain*> &sample){
	iteration++;

	if(sample.size()==0)return 0;
	CStrain *cfather=CommonFatherForSample(sample);
	assert(cfather!=NULL);
	int d=0;
	for(unsigned int i=0;i<sample.size();i++){
		d+=seg_distance(sample.at(i), cfather);	
	}

	return d;
}

class CDiversityOut{
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
		chosen=unif5(eng);
		sample.push_back(p_node->pathogens.at(chosen));
	}

	int ss=NSegSites(sample);

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

	return CDiversityOut(dist, ss);
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

void UpdateSus(){

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
		if(versionD) {
			assert( (*it)->migrantsV.size()==0 );
		}
		else {
			assert( (*it)->migrants.size()==0 );
		}
		if ( (*it)->pathogens.size()>0) { assert ( (*it)->state==INF ); }
	}

}

void MigrationGF(){

	for (unsigned int i=0; i < model.system_state.at(INF).size(); i++){
		CNode* p_node=model.system_state.at(INF).at(i);

		//fully mixed Diana's version
		std::tr1::poisson_distribution<double> poisson( Nd*mig_rate*model.system_state.at(SUS).size() );
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

		unsigned int nrecipients=model.system_state.at(SUS).size();

		if(nrecipients==0 || p_node->migrants.size()==0) {
			p_node->migrants.clear();
			continue;
		}

		std::tr1::uniform_int<int> unif2(1,nrecipients);

		map<int,CStrain*>::iterator itt;

		for(itt=p_node->migrants.begin(); itt != p_node->migrants.end(); itt++){
			int chosen=unif2(eng);
			model.system_state.at(SUS).at(chosen-1)->pathogens.push_back(itt->second);
		}
		p_node->migrants.clear();
	}
}

void MigrationGFwR(){

	for (unsigned int i=0; i < model.system_state.at(INF).size(); i++){
		CNode* p_node=model.system_state.at(INF).at(i);

		//fully mixed Diana's version
		std::tr1::poisson_distribution<double> poisson( Nd*mig_rate*model.system_state.at(SUS).size() );
		//fully mixed my version
		//std::tr1::poisson_distribution<double> poisson( Nd*mig_rate*(pop-1) );

		unsigned int nmigrants = (unsigned int)poisson(eng);
		
		assert( p_node->pathogens.size()==(unsigned int) Nd);
		
		unsigned int np=p_node->pathogens.size();
		std::tr1::uniform_int<int> unif1(0,np-1);
		while(p_node->migrantsV.size()<nmigrants){
			int chosen=unif1(eng);
			p_node->migrantsV.push_back( p_node->pathogens.at(chosen) );
		}
	}

	
	for (unsigned int i=0; i < model.system_state.at(INF).size(); i++){

		CNode* p_node=model.system_state.at(INF).at(i);

		unsigned int nrecipients=model.system_state.at(SUS).size();

		if(nrecipients==0 || p_node->migrantsV.size()==0) {
			p_node->migrantsV.clear();
			continue;
		}

		std::tr1::uniform_int<int> unif2(1,nrecipients);

		for(unsigned int k=0; k < p_node->migrantsV.size(); k++){
			int chosen=unif2(eng);
			model.system_state.at(SUS).at(chosen-1)->pathogens.push_back(p_node->migrantsV.at(k));
		}
		p_node->migrantsV.clear();
	}
}

void MigrationGnF(){

	for (unsigned int i=0; i < model.system_state.at(INF).size(); i++){
		CNode* p_node=model.system_state.at(INF).at(i);

		// my version
		//std::tr1::poisson_distribution<double> poisson( Nd*mig_rate*p_node->degree );
		// Diana's version
		std::tr1::poisson_distribution<double> poisson( Nd*mig_rate*p_node->count_neighbours_state(SUS) );

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

		unsigned int nrecipients=p_node->count_neighbours_state(SUS);

		if(nrecipients==0 || p_node->migrants.size()==0) {
			p_node->migrants.clear();
			continue;
		}

		std::tr1::uniform_int<int> unif2(1,nrecipients);

		map<int,CStrain*>::iterator itt;

		for(itt=p_node->migrants.begin(); itt != p_node->migrants.end(); itt++){
			int chosen=unif2(eng);
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
		p_node->migrants.clear();
	}
}

void MigrationGnFwR(){

	for (unsigned int i=0; i < model.system_state.at(INF).size(); i++){
		CNode* p_node=model.system_state.at(INF).at(i);

		// my version
		//std::tr1::poisson_distribution<double> poisson( Nd*mig_rate*p_node->degree );
		// Diana's version
		std::tr1::poisson_distribution<double> poisson( Nd*mig_rate*p_node->count_neighbours_state(SUS) );
		unsigned int nmigrants = (unsigned int)poisson(eng);
		assert( p_node->pathogens.size()==(unsigned int) Nd);
		
		unsigned int np=p_node->pathogens.size();
		std::tr1::uniform_int<int> unif1(0,np-1);
		while(p_node->migrantsV.size()<nmigrants){
			int chosen=unif1(eng);
			p_node->migrantsV.push_back( p_node->pathogens.at(chosen) );
		}

	}

	for (unsigned int i=0; i < model.system_state.at(INF).size(); i++){

		CNode* p_node=model.system_state.at(INF).at(i);

		unsigned int nrecipients=p_node->count_neighbours_state(SUS);

		if(nrecipients==0 || p_node->migrantsV.size()==0) {
			p_node->migrantsV.clear();
			continue;
		}

		std::tr1::uniform_int<int> unif2(1,nrecipients);

		for(unsigned int k=0; k < p_node->migrantsV.size(); k++){
			int chosen=unif2(eng);
			list<CNode*>::iterator it; int count=0;
			for(it=p_node->neighbours.begin(); count!=chosen; it++) {
				if( (*it)->state == SUS ){
					count++;
					if(count==chosen){
						break;
					}
				}
			}
			(*it)->pathogens.push_back(p_node->migrantsV.at(k));
		}
		p_node->migrantsV.clear();
	}
}

void MigrationDF(){

	// Diana's version
	for (unsigned int i=0; i < model.system_state.at(INF).size(); i++){
		CNode* p_node=model.system_state.at(INF).at(i);
		assert( p_node->pathogens.size()==(unsigned int) Nd);
		for(unsigned int j=0; j < p_node->pathogens.size(); j++){
			if( unif(eng) < mig_rate ){
				p_node->migrantsV.push_back( p_node->pathogens.at(j) );
			}
		}
	}

	for (unsigned int i=0; i < model.system_state.at(INF).size(); i++){

		CNode* p_node=model.system_state.at(INF).at(i);

		unsigned int nrecipients=model.system_state.at(SUS).size()+model.system_state.at(INF).size()-1;

		if(nrecipients==0 || p_node->migrantsV.size()==0) {
			p_node->migrantsV.clear();
			continue;
		}

		for(int l=SUS; l<REC; l++){			
			for(unsigned int j=0; j < model.system_state.at(l).size(); j++){
				if(l==SUS) {assert ( model.system_state.at(SUS).at(j) != p_node );}
				if (l==INF && model.system_state.at(INF).at(j) == p_node ){continue;}
				for(unsigned int k=0; k < p_node->migrantsV.size(); k++){
					model.system_state.at(l).at(j)->pathogens.push_back(p_node->migrantsV.at(k));
				}
			}
		}	
		p_node->migrantsV.clear();
	}
}

void MigrationDnF(){

	// Diana's version
	for (unsigned int i=0; i < model.system_state.at(INF).size(); i++){
		CNode* p_node=model.system_state.at(INF).at(i);
		assert( p_node->pathogens.size()==(unsigned int) Nd);
		for(unsigned int j=0; j < p_node->pathogens.size(); j++){
			if( unif(eng) < mig_rate ){
				p_node->migrantsV.push_back( p_node->pathogens.at(j) );
			}
		}
	}

	for (unsigned int i=0; i < model.system_state.at(INF).size(); i++){

		CNode* p_node=model.system_state.at(INF).at(i);

		unsigned int nrecipients=p_node->count_neighbours_state(SUS)+p_node->count_neighbours_state(INF);

		if(nrecipients==0 || p_node->migrantsV.size()==0) {
			p_node->migrantsV.clear();
			continue;
		}

		list<CNode*>::iterator it;
		for(it=p_node->neighbours.begin(); it!=p_node->neighbours.end(); it++){
			if( (*it)->state == SUS || (*it)->state == INF ){
				for(unsigned int j=0; j < p_node->migrantsV.size(); j++){
					(*it)->pathogens.push_back( p_node->migrantsV.at(j) );
				}
			}
		}	
		p_node->migrantsV.clear();
	}
}

void Migration(){
	if(versionD){
		if(fullymixed){
			MigrationDF();
		}
		else{
			MigrationDnF();
		}
	}
	else {
		if(fullymixed){
			if(withreplacement){
				MigrationGFwR();
			}
			else{
				MigrationGF();
			}
		}
		else{
			if(withreplacement){
				MigrationGnFwR();
			}
			else{
				MigrationGnF();
			}
		}
	}
	UpdateSus();
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
