// if Cf is small (ex. 0.1) then all mutations get fixed and the dynamics goes to extinction very slowly
// if Cf is bigger then dynamics goes to extinction quite fast
// if Cf is very big then there are no fixations

// if Cf is dist and small (ex. 0.1) then all mutations get fixed and the dynamics goes to extinction very slowly
// 1334441118 Cf=1.5, cp=0.015

#ifndef PARAMETERS_H
#define PARAMETERS_H

const int N0=10;
const double nPop=1e+8;
const double R0=3.;
const unsigned int inf_period=4;
const double nu=365.0/(double)inf_period; // units 1/year
const int Lep=120;
const int Lnep=160;
const int Lne=300;
const int L=Lep+Lnep+Lne;//580;//120+160+300=580; 60 epitope codons, 80 non-epitope codons and 300 neutral sites
const double mu=5.8*1e-3; // units 1/year
double mut_rate=mu*L/365.;
const double beta=R0*nu;
const double f0=nu*(R0-1.); // beta-nu 1.7
const double beta0=beta/(nPop*(beta-nu)); // 0.00001/f0
const double dt=1./(365.);
const double tMax=77.999;//38.999;
const int nbins=100;
const int timestart=1;
const int timeend=41;
const int maxmut=1000;

unsigned int iTimeMax=10000;
unsigned int iTime=0;
double a=5;
double Ub=0.001;


// ./plot.sh [0:10] 1:2 single00*

// git stash -- to delete current
// git pull reza-github ganna:ganna -- to download
// git log -- to check the version

// git commit -a -m "your comment"
// git push origin ganna:ganna

// git config --global color.log
// git config --global color.diff always

// git add signal.h

// git checkout branch

// ./strain &
// kill -30 27241 (the latter is the process number)
// ps
// fg

//gawk '{ print $0 > "line"NR}' fitness

//:wq

//mplayer N_distribution.avi -loop 5
//mplayer N_distribution.avi -speed 0.5

// ./strain 4 or ./strain 4 0 1st run
// ./strain 4 1 2nd run

#endif /* PARAMETERS_H */
