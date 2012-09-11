#include<iostream>
#include<math.h>
#include<stdlib.h>
using namespace std;

int main(int argn, char **args){

if(argn!=4){
	cerr<<"ERROR: Usage: "<<args[0]<<" begin end number_of_points"<<endl;
	exit(1);
	}

double rstart=atof(args[1]), rend=atof(args[2]);
double N=atof(args[3]);

double factor=exp((log(rend)-log(rstart))/N);

double x=rstart;
for(int i=0; i<=N;i++){
	cout<<x<<" ";
	x*=factor;
	}
	cout<<endl;
return 0;
}
