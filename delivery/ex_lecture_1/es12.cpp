/*************************************************************************
**************************************************************************
    _/_/_/_/       _/_/_/_/_/       Numerical Simulation Laboratory
   _/     _/      _/               Physics Department
  _/     _/      _/  _/_/_/       Università degli Studi di Milano
 _/     _/ _    _/      _/  _    Daniele Gilio
_/_/_/_/  /_/  _/_/_/_/_/  /_/  email: daniele.gilio@studenti.unimi.it
**************************************************************************
**************************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std; 
 
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
	
	//Initialize variables and outputs
	
	unsigned int realizations = 1E4;
	unsigned int N = 1;
	double sum_st = 0., sum_exp = 0., sum_l = 0.;
	double lambda =1., mu = 0., gamma = 1.; 	
	
	ofstream standard;
	ofstream exponential;
	ofstream lorentz;
	
	//Start of the actual program
	//N=1
	standard.open("standard_N_1.dat");
	exponential.open("exp_N_1.dat");
	lorentz.open("lorentz_N_1.dat");
	
	
	for(unsigned int i=0;i<realizations;i++){
		standard << rnd.Rannyu() << endl;
		exponential << rnd.Exp(lambda) << endl;
		lorentz << rnd.Cauchy_Lorentz(gamma,mu) << endl;
	}
	
	standard.close();
	exponential.close();
	lorentz.close();
	
	//N=2
	N=2;
	standard.open("standard_N_2.dat");
	exponential.open("exp_N_2.dat");
	lorentz.open("lorentz_N_2.dat");
	
	
	for(unsigned int i=0;i<realizations;i++){
		sum_st = 0.;
		sum_exp = 0.;
		sum_l = 0.;
		for(unsigned int j=0;j<N;j++){
			sum_st+=rnd.Rannyu();
			sum_exp+=rnd.Exp(lambda);
			sum_l+=rnd.Cauchy_Lorentz(gamma,mu);
		}
		standard << sum_st/N << endl;
		exponential << sum_exp/N << endl;
		lorentz << sum_l/N << endl;
		
	}
	
	standard.close();
	exponential.close();
	lorentz.close();
	
	//N=10
	N=10;
	standard.open("standard_N_10.dat");
	exponential.open("exp_N_10.dat");
	lorentz.open("lorentz_N_10.dat");
	
	
	for(unsigned int i=0;i<realizations;i++){
		sum_st = 0.;
		sum_exp = 0.;
		sum_l = 0.;
		for(unsigned int j=0;j<N;j++){
			sum_st+=rnd.Rannyu();
			sum_exp+=rnd.Exp(lambda);
			sum_l+=rnd.Cauchy_Lorentz(gamma,mu);
		}
		standard << sum_st/N << endl;
		exponential << sum_exp/N << endl;
		lorentz << sum_l/N << endl;
	}
	
	standard.close();
	exponential.close();
	lorentz.close();
	
	//N=100
	N=100;
	standard.open("standard_N_100.dat");
	exponential.open("exp_N_100.dat");
	lorentz.open("lorentz_N_100.dat");
	
	
	for(unsigned int i=0;i<realizations;i++){
		sum_st = 0.;
		sum_exp = 0.;
		sum_l = 0.;
		for(unsigned int j=0;j<N;j++){
			sum_st+=rnd.Rannyu();
			sum_exp+=rnd.Exp(lambda);
			sum_l+=rnd.Cauchy_Lorentz(gamma,mu);
		}
		standard << sum_st/N << endl;
		exponential << sum_exp/N << endl;
		lorentz << sum_l/N << endl;
	}
	
	standard.close();
	exponential.close();
	lorentz.close();
	
	return 0;
}	
