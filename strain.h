#ifndef STRAIN_H
#define STRAIN_H
#include <vector>
#include<assert.h>
#include"parameters.h"
#include"tools.h"

using namespace std;


class CStrain{
	public:
	explicit CStrain(int i, CStrain *f, double im_d=1);
	~CStrain();
	void SetAlive();
	double sumM();
	//double WeightedSumM(double chi(double) );
	int count_neigh();
	void die();
	void trim();
	void trim_links();
	CStrain *father();
	void add_neighbour(CStrain *ps, int d, double im_d=1){neighbours.push_back(CLink<CStrain>(ps, d, im_d));is_leaf=false;};
	void add_link(CStrain *ps, int d, double im_d=1){links.push_back(CLink<CStrain>(ps, d, im_d));};
	void get_diversity(double &div, size_t distance=0, CStrain *exclude=NULL);
	double WeightedSumS(double kapa(double), double distance=0, CStrain *exclude=NULL);
	double WeightedSumSup(double kapa(double), double distance=0, bool upward_only=false,  CStrain *original=NULL,CStrain *exclude=NULL);
	
	double WeightedSumI(double kapa(double), double distance=0, CStrain *original=NULL, CStrain *exclude=NULL);
	void calSubMeanFitness(double &sumfitness, double &totaln);
	bool delete_dead_branches();
	void remove_dead_children();
	double calSubN();
	unsigned int NNeigh();
	double sumNeighI();
	void print(ostream &out);
	void print_node(ostream &out)const;
	void print(ostream &out, double x, double y);
	double cal_print_widths();
	COffset &cal_offsets();
	
	//static unsigned int max_dist;
	std::vector<CLink<CStrain> > neighbours;
	std::vector<CLink<CStrain> > links;
	COffset offset;
	double print_width;
	bool is_leaf, dead;
	int ID;
	unsigned int visited;
	int NCopies;
	float x,y;
	float color;
	//double I, dI,dII, dSS, S, dS;
	int gen;
	private:
};

const double base_print_width=0.005;

unsigned int CStrain::NNeigh(){
	if(neighbours.at(0).head==NULL){
		return neighbours.size()-1;
		}
	return neighbours.size();
	}

//this is just used in function get_infected()
//to find out what was the maximum distance
//ever reached between the alive nodes

//Constructor take an int for ID and the point of the
//father node
CStrain::CStrain(int i, CStrain *f, double im_d){
	NCopies=0;
	gen=0;
	visited=-1;
	print_width = base_print_width;
	ID=i;
	//dI=dII=dS=dSS=I=S=0;
	if(f!=NULL){
		f->add_neighbour(this, 1, im_d);
		f->add_link(this,1, im_d);
		gen=f->gen+im_d;
	}
	add_neighbour(f,1, im_d);
	add_link(f,1, im_d);
	

	dead=true;
	is_leaf=true;
	if(ID>=0) {
		SetAlive();
	}

	x=y=0;
	color=0;
}

void CStrain::SetAlive(){
	dead=false;
}

//Cleans the allocated memory if not yet cleaned 
CStrain::~CStrain(){
}

//Returns the pointer to the father if not root node
CStrain* CStrain::father(){
	return neighbours.at(0).head;
}


//Returns the sum of N of all strains at distances up to rmax
//sumN[0] is the N of the strain itself
//sumN[1] is the sum of N's of all strains at distance 1
//....
//sumN[rmax] is the sum of N's of all strains at distance rmax

/*
double CStrain::WeightedSumI(double kapa(double), double distance, CStrain *original, CStrain *exclude){
	//if(distance>=rmax) return 0;
	//if(distance>max_dist)max_dist=distance;
	double weightedsum=kapa(distance)*I;
	if(original!=NULL)dSS+=original->S*kapa(distance);
	
	
	//comment not to go through the father -- infection forward only
	if(links.at(0).head!=NULL and links.at(0).head!=exclude) weightedsum+=links.at(0).head->WeightedSumI(kapa, distance+links.at(0).immune_distance, original,this);

	//if(is_leaf) return weightedsum;

	for(size_t i=1; i<links.size(); i++){
		if(links.at(i).head==exclude) continue;
		weightedsum+=links.at(i).head->WeightedSumI(kapa, distance+links.at(i).immune_distance,original, this);
	}
	return weightedsum;
}


double CStrain::WeightedSumS(double kapa(double), double distance, CStrain *exclude){
	//if(distance>=rmax) return 0;
	//if(distance>max_dist)max_dist=distance;
	double weightedsum=kapa(distance)*S;

	
	if(links.at(0).head!=NULL and links.at(0).head!=exclude)
		weightedsum+=links.at(0).head->WeightedSumS(kapa, distance+links.at(0).immune_distance, this);

	//if(is_leaf) return weightedsum;

	for(size_t i=1; i<links.size(); i++){
		if(links.at(i).head==exclude) continue;
		weightedsum+=links.at(i).head->WeightedSumS(kapa, distance+links.at(i).immune_distance, this);
	}
	return weightedsum;
}


double CStrain::WeightedSumSup(double kapa(double), double distance, bool upward_only, CStrain * original, CStrain *exclude){
	//if(distance>=rmax) return 0;
    	//if(distance>max_dist)max_dist=distance;
	distance=fabs(gen-original->gen);
    	double weightedsum=kapa(distance)*S;
	if(original!=NULL)dII+=original->I*kapa(distance);

   
    	if(links.at(0).head!=NULL and links.at(0).head!=exclude)
        	weightedsum+=links.at(0).head->WeightedSumSup(kapa, distance+links.at(0).immune_distance, false, original, this);

    	//if(upward_only) return weightedsum;

	if(gen>=original->gen)return weightedsum;

	// check if this function can be substituted by a simple loops over the strains and a check if(gen>=original->gen)

	    //if(is_leaf) return weightedsum;

    	for(size_t i=1; i<links.size(); i++){
        	if(links.at(i).head==exclude) continue;
        	weightedsum+=links.at(i).head->WeightedSumSup(kapa, distance+links.at(i).immune_distance, false, original, this);
    	}
    	return weightedsum;
}
*/

void CStrain::get_diversity(double &diversity, size_t distance, CStrain *exclude){
	//diversity+=N*distance;

	for(size_t i=0; i<links.size(); i++){
		if(links.at(i).head==NULL) continue;
		if(links.at(i).head==exclude) continue;
		links.at(i).head->get_diversity(diversity, distance+links.at(i).length, this);
	}
	return;
}

//Returns the number of alive neighbours of the strain
int CStrain::count_neigh(){
	int count=0;
	for (size_t i=0; i<neighbours.size();i++){
		if(neighbours.at(i).head==NULL) continue;
		if(! neighbours.at(i).head->dead) count++;
	}
	return count;
}

// To save some memory we clean M of 
// dead strains
void CStrain::die(){
	assert(!dead);
	dead=true;
}


void CStrain::trim_links(){

       int alive_branches=0;
       for(size_t i=1; i<links.size(); i++){
               if(!(links.at(i).head->is_leaf) or !(links.at(i).head->dead)){
                       links.at(i).head->trim_links();
                       alive_branches++;
               }
       }
       if(alive_branches==0 and dead){
                is_leaf=true;
                 if(links.size()>1) color=2;//color for dead branch
       }
}


void CStrain::remove_dead_children(){
	if( neighbours.size()==1) return;
	vector<CLink<CStrain> > newneighbours;
	newneighbours.push_back(neighbours.at(0));

	//delete children
	std::vector<CLink<CStrain> >::iterator it=neighbours.begin();
	for(it++; it!=neighbours.end(); it++){
		if(!(*it).head->dead) 
			newneighbours.push_back(*it);
	}

	neighbours=newneighbours;
	if(neighbours.size()==1) is_leaf=true;

	return;
}

bool CStrain::delete_dead_branches(){
	std::vector<CLink<CStrain> >::iterator it=neighbours.begin();
	for(it++; it!=neighbours.end(); it++){
		(*it).head->dead=(*it).head->delete_dead_branches();
	}
	remove_dead_children();
	if(neighbours.size()==1 and NCopies==0){
		return true;
		}
	return false;

}

void CStrain::trim(){
	if(is_leaf)return;
	std::vector<CLink<CStrain> >::iterator it, it0;
	int alive_branches=0;
	assert(neighbours.size()>1);
	for(it=neighbours.begin(), it++; it!=neighbours.end(); it++){
		if(!(*it).head->is_leaf or !(*it).head->dead){
			(*it).head->trim();
			alive_branches++;
		}
	}
	if(alive_branches==0 and dead){
		is_leaf=true;
		if(neighbours.size()>1) color=2;//color for dead branch
	}

}
/*
double CStrain::sumNeighI(){
	std::vector<CLink<CStrain> >::iterator it;
	double totI=0;
	int nn=0;
	for(it=neighbours.begin(); it!=neighbours.end(); it++){
		if((*it).head!=NULL){totI+=(*it).head->I; nn++;}
		}
	//totI-=(double)nn*I;
	return totI;
}
*/
/*
//prints the whole tree starting from this node
void CStrain::print(ostream &out){
	out<<ID<< "   "<<N<<"   "<<neighbours.size()-1<<"  ";
	std::vector<CLink<CStrain> >::iterator it;
	for(it=neighbours.begin(); it!=neighbours.end(); it++){
		if((*it).head==NULL)continue;
		out<<(*it).head->ID<< "   ";
	}
	out<<endl;
	for(it=neighbours.begin(), it++; it!=neighbours.end(); it++){
		(*it).head->print(out);
	//	out<<endl;
	}
}
*/

//prints N, number of neighbours and their IDs
/*
void CStrain::print_node(ostream &out)const{
	
	out<<ID<< "   "<<N<<"   ";
	out<<fitness<<"  ";
	if(!dead) {
		assert(N>0);
		for(size_t i=0; i<rmax+1; i++){
	 		out<<M[i]<<"  ";
		}
	}

	if(is_leaf) {
		out<<0<<endl;	
		return;
	}

	std::vector<CStrain*>::const_iterator it=neighbours.begin();
	int neigh=neighbours.size();
	if((*it)==NULL)neigh--;
	out<<neigh<<"  ";
	for(it=neighbours.begin(); it!=neighbours.end(); it++){
		if((*it)==NULL)continue;
		out<<(*it)->ID<< "   ";
	}
	out<<endl;
}
*/
void CStrain::print(ostream &out, double x, double y){
	static double dL=0.02;
	out<<"c "<<x<<"  "<<y<<"  "<<base_print_width/2<<"  "<<color<<endl;
	this->x=x;this->y=y;

	if(neighbours.size()<2) return;
	size_t nn=neighbours.size()-1.;
	size_t ind[neighbours.size()];
	mysort(neighbours,ind, neighbours.size());
	if(nn==1){
		out<<"l "<<x<<"  "<<y<<"  "<<x+dL<<"   "<<y<<endl;
		neighbours.at(ind[1]).head->print(out,x+dL, y);
		return;
	}
	double L=offset.width();
	L-=neighbours.at(ind[1]).head->offset.top;
	L-=neighbours.at(ind[nn]).head->offset.bottom;
	y+=L/2.;
	out<<"l "<<x<<"  "<<y<<"  ";
	out<<x<<"   "<<y-L<<endl;
	
	for(size_t i=1; i<=nn; i++){
		out<<"l "<<x<<"  "<<y<<"  "<<x+dL<<"   "<<y<<endl;
		neighbours.at(ind[i]).head->print(out,x+dL, y);
		y-=neighbours.at(ind[i]).head->offset.bottom;
		if(i<nn)y-=neighbours.at(ind[i+1]).head->offset.top;
	}
}
COffset &CStrain::cal_offsets(){
	offset.top=0.0;
	offset.bottom=0.0;
	size_t n=neighbours.size()-1;//no. of children
	cal_print_widths();
	if(n==0){
		offset.top=print_width/2;
		offset.bottom=print_width/2;
		return offset;
		}
	
	size_t ind[neighbours.size()];
	mysort(neighbours,ind, neighbours.size());

	COffset &off1=neighbours.at(ind[1]).head->cal_offsets();
	if(n==1){
		offset=off1;
		return offset;
		}
	
	double size=0;
	size+=off1.bottom;
	for(size_t i=2; i<n; i++){
		COffset &off=neighbours.at(ind[i]).head->cal_offsets();
		size+=off.width();
		}
	COffset &off2=neighbours.at(ind[n]).head->cal_offsets();
	size+=off2.top;

	offset.top=size/2+off1.top;
	offset.bottom=size/2+off2.bottom;

	return offset;
}

double CStrain::cal_print_widths(){
	if(neighbours.size()==1)return base_print_width;
	
	print_width=0;
	for(size_t i=1; i<neighbours.size(); i++)
		print_width+=neighbours.at(i).head->cal_print_widths();

	return print_width;
}

#endif
