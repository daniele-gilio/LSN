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

const double pi = M_PI;

double psi(double x, double mu, double sigma){
	return exp(-pow(x-mu,2.)/(2.*pow(sigma,2.)))+exp(-pow(x+mu,2.)/(2.*pow(sigma,2.)));
}

double v(double x){
	return pow(x,4.)-2.5*pow(x,2.);
}

double d2_psi(double x, double mu, double sigma){
	return exp(-pow(x-mu,2.)/(2.*pow(sigma,2.)))*(pow(x-mu,2.)/pow(sigma,4.)-pow(sigma,-2.))+exp(-pow(x+mu,2.)/(2.*pow(sigma,2.)))*(pow(x+mu,2.)/pow(sigma,4.)-pow(sigma,-2.));
}

double psi_2(double x, double mu, double sigma){
	return pow(psi(x,mu,sigma),2.);
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

  //Block Variables
  unsigned int M = 1E6;
  unsigned int block_number = 100;
  unsigned int block_size = M/block_number;

  //Simulation Variables
	double x = 1., xnew;
	double mu, sigma;
  double delta = 2.7;

	//Data Blocking
  double *pot = new double [block_number];
  double *pot_2 = new double [block_number];
  double *prog0 = new double [block_number];
  double *prog0_2 = new double [block_number];
  double *s_pot = new double [block_number];
  double px = 0.;
  double py = 0.;
  double judge = 0.;
  double a = 0.;

	//Set starting points
 	x=0.;
	mu=0.75;
	sigma=0.7;

 	ofstream out;

	//Optimize Variables
	double eMin=1E5, mu_opt=0., sigma_opt=0.;
	double accepted;

	//Optimize Algorithm
	cout << "Optimizing..." << endl;
	for(mu=0.75;mu<0.8;mu+=0.01){
		for(sigma=0.55;sigma<0.65;sigma+=0.01){
			accepted = 0.;
			for(unsigned int i=0;i<block_size*block_number;i++){
				xnew=rnd.Rannyu(x-delta,x+delta);
				//Metropolis
				px=pow(psi(x,mu,sigma),2.);
	   		py=pow(psi(xnew,mu,sigma),2.);
		   	a=min(1.,(py/px));
		   	judge = rnd.Rannyu();
		   	if(judge<=a){
		   			x=xnew;
						accepted++;
				}

	      pot[0]+=(-0.5*d2_psi(x,mu,sigma)+v(x)*psi(x,mu,sigma))/(psi(x,mu,sigma));
			}
			pot[0]/=block_number*block_size;
			if(pot[0]<eMin){
				eMin=pot[0];
				mu_opt=mu;
				sigma_opt=sigma;
			}

		}

	}

	cout << "Acceptance Rate for Optimized Values: " << accepted*100./double(block_size*block_number) << "%" << endl;
	cout << "Optimized Mu: " << mu_opt << endl;
	cout << "Optimized Sigma: " << sigma_opt << endl;

	mu=mu_opt;
	sigma=sigma_opt;
	x=0.;
	int nbin = 400;
	double *hist = new double [nbin];
	double x_min=-5., x_max=5.;
	double bin_size = (x_max-x_min)/double(nbin);
	out.open("potential.dat");
	ofstream h;
	h.open("hist.dat");
	accepted=0.;
	cout << "Performing Black Magic..." << endl;
	for(unsigned int i=0;i<block_number;i++){
	 	for(unsigned int j=0;j<block_size;j++){
	       xnew=rnd.Rannyu(x-delta,x+delta);
				//Uniform Metropolis
				px=psi_2(x,mu,sigma);
	   		py=psi_2(xnew,mu,sigma);
		   	a=min(1.,(py/px));
		   	judge = rnd.Rannyu();
		   	if(judge<=a)
		   			x=xnew;
	       pot[i]+=(-0.5*d2_psi(x,mu,sigma)+v(x)*psi(x,mu,sigma))/(psi(x,mu,sigma));

				 //h << x << endl;
				for(int l=0;l<nbin;l++)
					if(x>x_min+l*bin_size and x<x_min+(l+1)*bin_size){
						hist[l]+=1.;
						accepted++;
					}
   	}

		pot[i]/=double(block_size);
		pot_2[i]=pow(pot[i],2.);

		for(unsigned int k=0;k<i+1;k++){
				prog0[i]+=pot[k];
				prog0_2[i]+=pot_2[k];
			}

			prog0[i]/=double(i+1);
			prog0_2[i]/=double(i+1);


			s_pot[i]=sqrt((prog0_2[i]-pow(prog0[i],2))/double(i+1)); //i+1 is cheating, I know...

			out << i+1 << ";" << prog0[i] << ";" << s_pot[i] << endl;
	 }
	 	//h.close();
   	out.close();
		cout << accepted << endl;
		double z=0;
		out.open("hist.dat");
		for(int i=0;i<nbin;i++){
			out << hist[i]*double(nbin)/(10.*accepted) << endl;
			z+=hist[i]*double(nbin)/(10.*accepted);
		}

		cout << z << endl;
		out.close();

		cout << "Done" << endl;

	return 0;
}
