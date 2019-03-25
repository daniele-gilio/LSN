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

double f_uni(double x){//Uniform Distribution Integrand
	return pi*0.5*cos(pi*(x*0.5));
}

double f_imp(double x){//Importance sampling Integrand
	return f_uni(x)/(-2*x+2);
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
	unsigned int M = 1E6;
	unsigned int block_number = 100;
	unsigned int block_size = M/block_number;
	double x = 0.;
	double y = 0.;
	double *I_u = new double [block_number];//Uniform Integral
	double *I_u_2 = new double[block_number];//Squared Uniform Integral for statistical analysis
	double *I_i = new double [block_number];//Importance sampling Integral
	double *I_i_2 = new double[block_number];//Squared Importance sampling Integral for statistical analysis
	double *sigma_u = new double [block_number];//Uniform error
	double *sigma_i = new double [block_number];//Importance sampling error
	double *prog_sum_u = new double [block_number];
	double *prog_sum_i = new double [block_number];
	double *prog_sum2_u = new double [block_number];
	double *prog_sum2_i = new double [block_number];
	ofstream st_sample;
	ofstream imp_sample;
	
	for(unsigned int i=0;i<block_number;i++){
		I_u[i]=0;
		I_i[i]=0;
		I_u_2[i]=0;
		I_i_2[i]=0;
		sigma_u[i]=0;
		sigma_i[i]=0;
		prog_sum_u[i]=0;
		prog_sum_i[i]=0;
		prog_sum2_u[i]=0;
		prog_sum2_i[i]=0;	
	}
		
	//Start of the actual program 
	
	st_sample.open("st_sample.dat");
	imp_sample.open("imp_sample.dat");
	
	//Uniform Distribution and Importance Sampling
	for(unsigned int i=0;i<block_number;i++){
		//
		for(unsigned int j=0;j<block_size;j++){
			x=rnd.Rannyu();//<-Uniform Distribution
			y=1.-sqrt(1.-x);//<-Importance sampling with p(x)=-2x+2
			I_u[i]+=f_uni(x);
			I_i[i]+=f_imp(y);
		}
		I_u[i]/=block_size;
		I_i[i]/=block_size;
		I_u_2[i]=pow(I_u[i],2);
		I_i_2[i]=pow(I_i[i],2);
		//Compute cumulative mean and error
		for(unsigned int j=0;j<i+1;j++){
			prog_sum_u[i]+=I_u[j];
			prog_sum_i[i]+=I_i[j];
			prog_sum2_u[i]+=I_u_2[j];
			prog_sum2_i[i]+=I_i_2[j];
		}
		
		prog_sum_u[i]/=(i+1);
		prog_sum_i[i]/=(i+1);
		prog_sum2_u[i]/=(i+1);
		prog_sum2_i[i]/=(i+1);
		
		if(i==0){
			sigma_u[i]=0.;
			sigma_i[i]=0.;
		}
		else{
			sigma_u[i]=error(prog_sum_u,prog_sum2_u,i);
			sigma_i[i]=error(prog_sum_i,prog_sum2_i,i);
		}
		st_sample << (i+1)*block_size << ";" << prog_sum_u[i] << ";" << sigma_u[i] << endl;
		imp_sample << (i+1)*block_size << ";" << prog_sum_i[i] << ";" << sigma_i[i] << endl;
	}
	
	
	st_sample.close();
	imp_sample.close();
	
	return 0;
}

