#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

const double pi = M_PI;

double error(double* ave, double* av2, int n){
	if(n==0)
		return 0;
	else
		return sqrt((av2[n]-pow(ave[n],2))/double(n));
}
 
int main (int argc, char *argv[]){
	//Inizialization of the Random Number Generator
	Random rnd;
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();
	ifstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
			input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
			rnd.SetRandom(seed,p1,p2);
			}
		}
	input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
	rnd.SaveSeed();
   
	//Inizialize Variables, Arrays and Outputs
	unsigned int M = 1E4;
	unsigned int step = 100;
	
	int *x_disc = new int [3];
	double *x_cont = new double [3];
	double p = 0.;
	double theta = 0.;
	double fi = 0.;
	double x = 0.;
	double y = 0.;
	double *r_disc = new double [step];
	double *r_disc_2 = new double [step];
	double *r_cont = new double [step];
	double *r_cont_2 = new double [step];
	double *sigma_disc = new double [step];
	double *sigma_cont = new double [step];
	
	ofstream disc;
	ofstream cont;
	ofstream test;
	disc.open("disc.dat");
	cont.open("cont.dat");
	test.open("test.dat");
	int t = 0;
	
	for(unsigned int i=0;i<M;i++){
		x_disc[0]=0.;
		x_disc[1]=0.;
		x_disc[2]=0.;
		x_cont[0]=0.;
		x_cont[1]=0.;
		x_cont[2]=0.;
		for(unsigned int j=0;j<step;j++){
			//Generate Random Steps
			p=rnd.Rannyu(0.,6.); //<-Discreete
			//Continuum  
			theta=rnd.Rannyu(0., pi);
			x=rnd.Rannyu(-1.,1.);
			y=rnd.Rannyu(-1.,1.);
			while(pow(x,2)+pow(y,2)>1){//keeping the points inside the unit circle
				x=rnd.Rannyu(-1., 1.);
				y=rnd.Rannyu(-1.,1.);
			}
			
			if(y>=0)
				fi=acos(x/sqrt(pow(x,2)+pow(y,2)));
			else 
				fi=2*pi-acos(x/sqrt(pow(x,2)+pow(y,2)));
			
			//Fi should be correctly uniformly distributed in (0,2*pi) after the above code
			
			//Choosing the discreete random path
			if(p>=0 and p<1)
				x_disc[0]++;
			if(p>=1 and p<2)
				x_disc[0]--;
			if(p>=2 and p<3)
				x_disc[1]++;
			if(p>=3 and p<4)
				x_disc[1]--;
			if(p>=4 and p<5)
				x_disc[2]++;
			if(p>=5 and p<6)
				x_disc[2]--;
			
			//Cumulate the discreete mean	
			r_disc[j]+=sqrt(pow(x_disc[0],2)+pow(x_disc[1],2)+pow(x_disc[2],2));
			if(j==99){
				test << t << ";" << sqrt(pow(x_disc[0],2)+pow(x_disc[1],2)+pow(x_disc[2],2)) << endl;
				t++;
			}
			r_disc_2[j]+=pow(x_disc[0],2)+pow(x_disc[1],2)+pow(x_disc[2],2);
			
			//Choosing the continuum random path
			x_cont[0]+=cos(fi)*sin(theta);
			x_cont[1]+=sin(fi)*sin(theta);
			x_cont[2]+=cos(theta);
			
			//Cumulate the continuum mean
			r_cont[j]+=sqrt(pow(x_cont[0],2)+pow(x_cont[1],2)+pow(x_cont[2],2));
			r_cont_2[j]+=pow(x_cont[0],2)+pow(x_cont[1],2)+pow(x_cont[2],2);
			
			//Finalize calculations
			if(i==M-1){
				r_disc[j]/=M;
				r_disc_2[j]/=M;
				sigma_disc[j]=sqrt((r_disc_2[j]-pow(r_disc[j],2)));
				//cout << pow(r_disc[j],2) << " " << r_disc_2[j] << " " << sigma_disc[j] <<  endl;
				r_cont[j]/=M;
				r_cont_2[j]/=M;
				sigma_cont[j]=sqrt((r_cont_2[j]-pow(r_cont[j],2)));
				
				disc << j << ";" << r_disc[j] << ";" << sigma_disc[j] << endl;
				cont << j << ";" << r_cont[j] << ";" << sigma_cont[j] << endl;
			}
		
		}
		
		
	}
	
	disc.close();
	cont.close();
	test.close();
	
	return 0;
}
