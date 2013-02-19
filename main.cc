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
std::tr1::uniform_real<double> unif(0,1);

//To to be used for assigning ID's
int stotal;
//List of all alive strains
vector<CStrain*> allstrains;
vector<CStrain*> strains;
CStrain *top=NULL;

int t=1;
int tmax=1000000;
int tstep=1;
int tprint=5;
int transient=2000;

bool versionD=true;
bool withreplacement=false;

bool fullymixed=true;
bool correct=true;
bool SIR=false;
bool SIRS=true;
bool SI=false;

bool DivLocFullSamp=true;
bool DivGlobFullSamp=true;

ofstream outorigfixtimes("OrigFixTimes");
ofstream outsusinfrec("SusInfRec");
ofstream outdiversity("Diversity");
ofstream outmutant("Mutant");
ofstream outallstrains("AllStrains");
//ofstream outstrains("Strains");

CStrain* mutant=NULL;

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
	CDiversityOut(double d=0, double n=0):distance(d),SegSites(n){}
	double distance;
	double SegSites;
 	private:
};

CDiversityOut GeneticDiversityGlobal(){

	vector<CStrain*> sample;

	if(DivGlobFullSamp){
		sample=strains;
	}
	else{	
		std::tr1::uniform_int<int> unif4(0,model.system_state.at(INF).size()-1);
		std::tr1::uniform_int<int> unif5(0,Nd-1);

		for (int i=0; i<nt; i++){
			int chosen=unif4(eng);
			CNode* p_node=model.system_state.at(INF).at(chosen);
			chosen=unif5(eng);
			sample.push_back(p_node->pathogens.at(chosen));
		}
	}

	double ss=NSegSites(sample);

	double dist=0.;
	int k=0;

	if(DivGlobFullSamp){
		for(unsigned int i=0; i<strains.size(); i++){
			for(unsigned int j=i+1; j<strains.size(); j++){
				k=distance( sample.at(i), sample.at(j) );
				assert(k!=-1);
				dist+=k*(sample.at(i)->NCopies)*(sample.at(j)->NCopies);
			}
		}
	
		dist=dist/(Nd*inf*(Nd*inf-1.)/2.);
	}
	else{
		for(int i=0; i<nt; i++){
			for(int j=i+1; j<nt; j++){
				k=distance( sample.at(i), sample.at(j) );
				assert(k!=-1);
				dist+=k;
			}
		}

		dist=dist/(nt*(nt-1.)/2.);
	}

	return CDiversityOut(dist, ss);
}

double GeneticDiversityLocal(CNode* p_node){

	map<int,CStrain*> sample;
	assert(nd<=Nd);

	if(DivLocFullSamp){
		assert(nd==Nd);
		for(int i=0; i<Nd; i++){
			sample[i]=p_node->pathogens.at(i);
		}
	}
	else{
		std::tr1::uniform_int<int> unif6(0,Nd-1);
		while (sample.size()<(unsigned int) nd){
			int chosen=unif6(eng);
			sample[chosen]=p_node->pathogens.at(chosen);
		}
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

CDiversityOut GeneticDiversityLocalAverage(){

	double dist=0.;
	double ss=0.;

	for (unsigned int i=0; i < model.system_state.at(INF).size(); i++){
		dist+=GeneticDiversityLocal( model.system_state.at(INF).at(i) );
		ss+=NSegSites( model.system_state.at(INF).at(i)->pathogens );
	}

	dist=dist/(double)model.system_state.at(INF).size();
	ss=(double)ss/(double)model.system_state.at(INF).size();

	return CDiversityOut(dist, ss);
}

void InitialConditions(){

	stotal=0;
	top = new CStrain(stotal,NULL);
	allstrains.push_back(top);
	strains.push_back(top);
	stotal++;

	top->t_origination=t;

	for (unsigned int i=0; i < model.system_state.at(INF).size(); i++){
		for(int j=0; j < Nd; j++){
			model.system_state.at(INF).at(i)->add_pathogen(top);
		}
	}
}

void RecoveryBirthFromI(){

	int j=0;

	vector<CNode*> infectives;

	infectives=model.system_state.at(INF);

	for (unsigned int i=0; i < infectives.size(); i++){
		double rand_num = unif(eng);
		if( rand_num < rec_rate ){
			//cerr << rand_num << endl;
			inf--; rec++;
			j++;
			CNode* p_node=infectives.at(i);
			model.UpdateSystemState( p_node, INF, REC);
			p_node->pathogens.clear();
		}
		else {
			if( rand_num < (rec_rate+birth_rate) ){
				inf--; sus++;
				j++;
				CNode* p_node=infectives.at(i);
				model.UpdateSystemState( p_node, INF, SUS);
				p_node->pathogens.clear();
			}
		}
	}	

	//cerr << "number of RecoveryBirthFromI " << j << endl;
	//cerr << "number of infectives " << inf << endl;
	//cerr << "rec_rate+birth_rate " << rec_rate+birth_rate << endl; 
}

void BirthFromR(){

	vector<CNode*> recovereds;

	recovereds=model.system_state.at(REC);

	for (unsigned int i=0; i < recovereds.size(); i++){
		if( unif(eng) < birth_rate ){
			rec--; sus++;
			CNode* p_node=recovereds.at(i);
			model.UpdateSystemState( p_node, REC, SUS);
			assert( p_node->pathogens.size()==0);
			//p_node->pathogens.clear();
		}
	}	
}

void BirthDiseaseDeath(){
	
	vector<CNode*> infectives;

	infectives=model.system_state.at(INF);

	for (unsigned int i=0; i < infectives.size(); i++){
		if( unif(eng) < disease_death_rate ){
			inf--; sus++;
			CNode* p_node=infectives.at(i);
			model.UpdateSystemState( p_node, INF, SUS);
			p_node->pathogens.clear();
		}
	}

}

void ImmunityLossBirthFromR(){

	vector<CNode*> recovereds;

	recovereds=model.system_state.at(REC);

	for (unsigned int i=0; i < recovereds.size(); i++){

		if( unif(eng) < (wan_rate+birth_rate) ){
			rec--; sus++;
			CNode* p_node=recovereds.at(i);
			model.UpdateSystemState( p_node, REC, SUS);
			assert( p_node->pathogens.size()==0);
			//p_node->pathogens.clear();
		}
	}

}

void UpdateSus(){

	int j=0;
	vector<CNode*> susceptibles;
	susceptibles=model.system_state.at(SUS);
	for (unsigned int i=0; i < susceptibles.size(); i++){
		CNode* p_node=susceptibles.at(i);
		if( p_node->pathogens.size()>0 ){
	
			//cout<< t <<"  "<<p_node->pathogens.size()<<endl;

			for(unsigned int k=0; k < p_node->pathogens.size(); k++){
				//assert(p_node->pathogens.at(k)!=NULL);
				if(p_node->laststrain!=NULL && p_node->pathogens.at(k)!=NULL){
				//cout<< t <<"   "<<distance(p_node->laststrain,p_node->pathogens.at(k))<<endl;	
				}
			}

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
			model.system_state.at(SUS).at(chosen-1)->add_pathogen(itt->second);
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
			model.system_state.at(SUS).at(chosen-1)->add_pathogen(p_node->migrantsV.at(k));
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
			(*it)->add_pathogen(itt->second);
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
			(*it)->add_pathogen(p_node->migrantsV.at(k));
		}
		p_node->migrantsV.clear();
	}
}

void MigrationCorrectF(){
	// correct version

	unsigned int nrecipients=model.system_state.at(SUS).size();

	if(nrecipients==0) { return;}

	for (unsigned int i=0; i < model.system_state.at(INF).size(); i++){
		CNode* p_node=model.system_state.at(INF).at(i);
		assert( p_node->pathogens.size()==(unsigned int) Nd);
		for(unsigned int k=0; k < p_node->pathogens.size(); k++){
			for(unsigned int j=0; j < nrecipients; j++){
				assert ( model.system_state.at(SUS).at(j) != p_node );

				if(p_node->pathogens.at(k)!=mutant){
					if( unif(eng) < mig_rate ){
						model.system_state.at(SUS).at(j)->add_pathogen(p_node->pathogens.at(k));
					}
				}
				else{
					if( unif(eng) < mig_rate_mutant ){
						model.system_state.at(SUS).at(j)->add_pathogen(p_node->pathogens.at(k));
					}	
				}
			}
		}
	}
}

void MigrationCorrectnF(){
	// correct version

	for (unsigned int i=0; i < model.system_state.at(INF).size(); i++){
		CNode* p_node=model.system_state.at(INF).at(i);
		assert( p_node->pathogens.size()==(unsigned int) Nd);

		unsigned int nrecipients=p_node->count_neighbours_state(SUS);
		if(nrecipients==0) { continue;}

		for(unsigned int k=0; k < p_node->pathogens.size(); k++){

			list<CNode*>::iterator it;
			for(it=p_node->neighbours.begin(); it!=p_node->neighbours.end(); it++){
				if( (*it)->state == SUS){
					if(p_node->pathogens.at(k)!=mutant){
						if( unif(eng) < mig_rate ){
							(*it)->add_pathogen(p_node->pathogens.at(k));
						}
					}
					else{
						if( unif(eng) < mig_rate_mutant ){
							(*it)->add_pathogen(p_node->pathogens.at(k));
						}	
					}
				}
			}
		}
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
					model.system_state.at(l).at(j)->add_pathogen(p_node->migrantsV.at(k));
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
					(*it)->add_pathogen( p_node->migrantsV.at(j) );
				}
			}
		}	
		p_node->migrantsV.clear();
	}
}

void Migration(){

	if(correct){
		if(fullymixed){
			MigrationCorrectF();
		}
		else{
			MigrationCorrectnF();
		}
	}
	else{
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
				s->t_origination=t;
			}
			newgeneration.push_back(s);
		}
		p_node->pathogens.clear();
		p_node->pathogens=newgeneration;

		//assert(p_node->pathogens.at(0)!=NULL);
		assert(p_node->pathogens.size()==1); // works only for Nd=1
		p_node->laststrain=p_node->pathogens.at(0); // the first one is copied

		assert(p_node->pathogens.size()==(unsigned int)Nd);
	}

	//cerr << "number of mutations " << k << endl;

}

void IntroduceMutant(){

	std::tr1::uniform_int<int> unif4(0,model.system_state.at(INF).size()-1);
	std::tr1::uniform_int<int> unif5(0,Nd-1);

	int chosen=unif4(eng);
	CNode* p_node=model.system_state.at(INF).at(chosen);
	chosen=unif5(eng);
	//mutant=p_node->pathogens.at(chosen);

	CStrain *s = new CStrain(stotal, p_node->pathogens.at(chosen) );
	stotal++;
	allstrains.push_back(s);
	mutant=s;

	mutant->t_origination=t;

	p_node->pathogens.clear();

	for(int j=0; j < Nd; j++){
		p_node->add_pathogen(mutant);
	}

}

void UpdateDynamics(){

	if (SIR){ // childhood diseases
		RecoveryBirthFromI(); // I->R, I->S
		BirthFromR(); // R->S
		Migration(); // S->I
		Reproduction();
	}
	else{
		if (SIRS){ // influenza
			RecoveryBirthFromI(); // I->R, I->S
			ImmunityLossBirthFromR(); // R->S
			Migration(); // S->I
			Reproduction();
		}
		else{ // HIV
			if (SI) {
				BirthDiseaseDeath(); // I->S
				Migration(); // S->I
				Reproduction();
			}
		}
	}

}

void PrintSIR(){

	outsusinfrec << t <<"\t"<< sus/(double)model.network->get_N() <<"\t"<< inf/(double)model.network->get_N() <<"\t"<< rec/(double)model.network->get_N() << endl;

}

void PrintAllStrains(){

	outallstrains << t <<"\t"<< strains.size() <<"\t"<< allstrains.size() <<"\t"<< stotal << endl;
}

void PrintDiversity(){

	CDiversityOut GDLAverage=GeneticDiversityLocalAverage();
        CDiversityOut GDGlobal=GeneticDiversityGlobal();

	outdiversity << t <<"\t"<< GDGlobal.distance <<"\t"<< GDGlobal.SegSites <<"\t"<< GDLAverage.distance <<"\t"<< GDLAverage.SegSites << endl;

}

void PrintMutant(){

	if(mutant!=NULL && mutant->NCopies == 0) {
		outmutant << t <<"\t"<< mutant->NCopies/(double)(Nd*inf) << endl;
		cerr << "lost" << endl; 
		exit(0);
	}
	if(mutant!=NULL && mutant->NCopies == (Nd*inf) ) {
		outmutant << t <<"\t"<< mutant->NCopies/(double)(Nd*inf) << endl;
		cerr<< "fixation"<<endl; 
		exit(0);
	}

	if(mutant!=NULL) outmutant << t <<"\t"<< mutant->NCopies/(double)(Nd*inf) << endl;
}

void UpdateStrainsNCopies(){

	vector<CStrain *>::iterator it=allstrains.begin();
	for (; it != allstrains.end(); it++){
		(*it)->NCopies=0;
		(*it)->MemoryCopies=0;
	}

	for (unsigned int i=0; i < model.system_state.at(INF).size(); i++){
		it=model.system_state.at(INF).at(i)->pathogens.begin();
		for (; it != model.system_state.at(INF).at(i)->pathogens.end(); it++){
			(*it)->NCopies++;
		}
	}

	//PrintMutant();
	
	vector<CNode *>::iterator itt; //going through nodes of the network
	for (itt=model.network->nodes.begin(); itt!=model.network->nodes.end(); itt++) {
		CStrain* ss=(*itt)->laststrain;
		
		//cerr<<"before"<<endl;
		//ss->MemoryCopies++;
		//cerr<<"after"<<endl;

	}
}

void RemoveDeadStrains(){
	vector<CStrain *> newall;
	vector<CStrain *> newlive;
	vector<CStrain *>::iterator it=allstrains.begin();
	for (; it != allstrains.end(); it++){

		if( (*it)->NCopies>0 ) {
			newlive.push_back(*it);
		}
		if(!(*it)->dead){
			newall.push_back(*it);
			}
		else{
			delete *it;
			*it=NULL;
			}
	}
	allstrains=newall;
	strains=newlive;
}

void PrintOriginationsFixations(){
	
	vector<CStrain *>::iterator it=allstrains.begin();

	for (; it != allstrains.end(); it++){
		if((*it)->notprinted && !((*it)->notfixed)){
			(*it)->notprinted=false;
			outorigfixtimes<< (*it)->t_origination <<"\t"<< (*it)->t_fixation <<"\t"<< (*it)->ID <<endl;
		}
	}
	/*
	int tot=0;

	for (it=allstrains.begin(); it != allstrains.end(); it++){
		tot+=(*it)->NCopies;
	}
	*/
	//cerr<<tot-(inf*Nd)<< "\t"<< top->subtotal_ncopies() - (inf*Nd)<<endl;
}

void UpdateMemory(){

	UpdateStrainsNCopies();
	top->delete_dead_branches();
	RemoveDeadStrains();
	top->subtotal_ncopies( inf*Nd, t );
}

void Iterate(){

	while(t<=tmax and inf>0 and sus>=0 and rec>=0){

		UpdateMemory();

		PrintSIR();
		PrintAllStrains();
		PrintDiversity();
		PrintOriginationsFixations();
		PrintMutant();

		if(t%(tmax+1)==0) IntroduceMutant();

		UpdateDynamics(); //Dynamics running
		PrintSIR();
		//cerr << network->average_degree() <<"\t"<< model.state_degree(SUS) <<"\t"<< model.state_degree(INF) <<"\t" << model.state_degree(REC) <<endl;
		t+=tstep;
		
		assert( (unsigned int)(sus+inf+rec)==model.network->get_N());
		assert((unsigned int)sus==model.system_state.at(SUS).size() && (unsigned int)inf==model.system_state.at(INF).size() && (unsigned int)rec==model.system_state.at(REC).size() );
	}

}

int main(int argc, char **argv){

	if(argc>1)seed=atoi(argv[1]);
	//if(argc>2)r0=atof(argv[2]);
	eng.seed(seed);
	//cerr<<seed<<endl;

	model.Initial_Conditions();
	InitialConditions();

	Iterate();

return 0;

}

//ssh -v theoserv.phy.umist.ac.uk -l rozhnova
//1ynxg4v
//ssh -v thpc55
//git push origin master