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
   	
   	unsigned int M = 1E6;
   	unsigned int block_number = 100;
   	unsigned int block_size = M/block_number;
   	unsigned int test = 10000;
   	
   	double *x_0 = new double [3];
   	double *y_0 = new double [3];
   	double *x_1 = new double [3];
   	double *y_1 = new double [3];
   	double delta_0 = 1.2;
   	double delta_1 = 3.;
   	double a0 = 0.;
   	double ave_a0 = 0.;
   	double a1 = 0.;
   	double ave_a1 = 0.;
   	double p0x = 0.;
   	double p1x = 0.;
   	double p0y = 0.;
   	double p1y = 0.;
   	double judge_0 = 0.;
   	double judge_1 = 0.;
   	double *ave_r0 = new double [block_number];
   	double *ave_r0_2 = new double [block_number];
   	double *ave_r1 = new double [block_number];
   	double *ave_r1_2 = new double [block_number];
   	double *prog0 = new double [block_number];
   	double *prog1 = new double [block_number];
   	double *prog0_2 = new double [block_number];
   	double *prog1_2 = new double [block_number];
   	double *s_r0 = new double [block_number];
   	double *s_r1 = new double [block_number];
   	
   	//Set starting point
   	x_0[0]=1.;
   	x_0[1]=1.;
   	x_0[2]=1.;
   	
   	x_1[0]=1.;
   	x_1[1]=1.;
   	x_1[2]=1.;
   	
   	ofstream r_0, r_1;
   	r_0.open("test_run_r0.dat");
   	r_1.open("test_run_r1.dat");
   	
   	//Test run to check for the 50% acceptance rate
	for(unsigned int i=0;i<test;i++){
		for(unsigned int j=0;j<3;j++){
			y_0[j]=rnd.Rannyu(x_0[j]-delta_0,x_0[j]+delta_0);
			y_1[j]=rnd.Rannyu(x_1[j]-delta_1,x_1[j]+delta_1);
		}	
   		p0x=exp(-2.*sqrt(x_0[0]*x_0[0]+x_0[1]*x_0[1]+x_0[2]*x_0[2]))/pi;	   			
   		p0y=exp(-2.*sqrt(y_0[0]*y_0[0]+y_0[1]*y_0[1]+y_0[2]*y_0[2]))/pi;
	   	a0=min(1.,(p0y/p0x));
	   	judge_0 = rnd.Rannyu();
	   	if(judge_0<=a0)
	   		for(unsigned int j=0;j<3;j++)
	   			x_0[j]=y_0[j];
	   	ave_a0+=a0;
	   	//ave_r0+=sqrt(x_0[0]*x_0[0]+x_0[1]*x_0[1]+x_0[2]*x_0[2]);
	   	r_0 << x_0[0] << ";" << x_0[1] << ";" << x_0[2] << endl;
	   	
	   	p1x=exp(-sqrt(x_1[0]*x_1[0]+x_1[1]*x_1[1]+x_1[2]*x_1[2]))*x_1[2]*x_1[2]/(32.*pi);	   			
   		p1y=exp(-sqrt(y_1[0]*y_1[0]+y_1[1]*y_1[1]+y_1[2]*y_1[2]))*y_1[2]*y_1[2]/(32.*pi);
	   	a1=min(1.,(p1y/p1x));
	   	judge_1 = rnd.Rannyu();
	   	if(judge_1<=a1)
	   		for(unsigned int j=0;j<3;j++)
	   			x_1[j]=y_1[j];
	   	ave_a1+=a1;
	   	//ave_r1+=sqrt(x_1[0]*x_1[0]+x_1[1]*x_1[1]+x_1[2]*x_1[2]);
	   	r_1 << x_1[0] << ";" << x_1[1] << ";" << x_1[2] << endl;
	}
	
	r_0.close();
	r_1.close();   	
	ave_a0/=test/100.;
	ave_a1/=test/100.;
	//ave_r0/=test;
	//ave_r1/=test;
	   	
	cout << "The acceptance rate is: " << ave_a0 << "%" << endl;
	//cout << ave_r0 << endl;
	
	cout << ave_a1 << endl;
	//cout << ave_r1 << endl;
	
	r_0.open("r0.dat");
	r_1.open("r1.dat");
	
	x_0[0]=1.;
   	x_0[1]=1.;
   	x_0[2]=1.;
   	
   	x_1[0]=1.;
   	x_1[1]=1.;
   	x_1[2]=1.;   	
   	
   	for(unsigned int i=0;i<block_number;i++){
   		for(unsigned int j=0;j<block_size;j++){
   			for(unsigned int j=0;j<3;j++){
				y_0[j]=rnd.Rannyu(x_0[j]-delta_0,x_0[j]+delta_0);
				y_1[j]=rnd.Rannyu(x_1[j]-delta_1,x_1[j]+delta_1);
			}	
	   		p0x=exp(-2.*sqrt(x_0[0]*x_0[0]+x_0[1]*x_0[1]+x_0[2]*x_0[2]))/pi;	   			
	   		p0y=exp(-2.*sqrt(y_0[0]*y_0[0]+y_0[1]*y_0[1]+y_0[2]*y_0[2]))/pi;
		   	a0=min(1.,(p0y/p0x));
		   	judge_0 = rnd.Rannyu();
		   	if(judge_0<=a0)
		   		for(unsigned int j=0;j<3;j++)
		   			x_0[j]=y_0[j];
		  	
		   	ave_r0[i]+=sqrt(x_0[0]*x_0[0]+x_0[1]*x_0[1]+x_0[2]*x_0[2]);
		   	
		   	p1x=exp(-sqrt(x_1[0]*x_1[0]+x_1[1]*x_1[1]+x_1[2]*x_1[2]))*x_1[2]*x_1[2]/(32.*pi);	   			
	   		p1y=exp(-sqrt(y_1[0]*y_1[0]+y_1[1]*y_1[1]+y_1[2]*y_1[2]))*y_1[2]*y_1[2]/(32.*pi);
		   	a1=min(1.,(p1y/p1x));
		   	judge_1 = rnd.Rannyu();
		   	if(judge_1<=a1)
		   		for(unsigned int j=0;j<3;j++)
		   			x_1[j]=y_1[j];
		   	ave_a1+=a1;
		   	ave_r1[i]+=sqrt(x_1[0]*x_1[0]+x_1[1]*x_1[1]+x_1[2]*x_1[2]);
		   		
	   	}
	   	
	   	ave_r0[i]/=double(block_size);
	   	ave_r0_2[i]=pow(ave_r0[i],2.);
	   	ave_r1[i]/=double(block_size);
	   	ave_r1_2[i]=pow(ave_r1[i],2.);
	   	
	   	for(unsigned int k=0;k<i+1;k++){
			prog0[i]+=ave_r0[k];
			prog1[i]+=ave_r1[k];
			prog0_2[i]+=ave_r0_2[k];
			prog1_2[i]+=ave_r1_2[k];
		}
		
		prog0[i]/=double(i+1);
		prog1[i]/=double(i+1);
		prog0_2[i]/=double(i+1);
		prog1_2[i]/=double(i+1);
		
		
		if(i==0){
			s_r0[i]=0.;
			s_r1[i]=0.;
		}
		
		else{
			s_r0[i]=sqrt((prog0_2[i]-pow(prog0[i],2))/double(i));
			s_r1[i]=sqrt((prog1_2[i]-pow(prog1[i],2))/double(i));
		}
		
		r_0 << i+1 << ";" << prog0[i] << ";" << s_r0[i] << endl;
		r_1 << i+1 << ";" << prog1[i] << ";" << s_r1[i] << endl;
   	}
   	
   	r_0.close();
   	r_1.close();
   	
   	//Gauss, can be done in the upper cycles, we need to double the variables though. This way is just a little slower computationally, but faster to code. Anyway, don't forget to correct it.
   	
   	r_0.open("r0_gauss.dat");
   	r_1.open("r1_gauss.dat");
   	
   	for(unsigned int i=0;i<block_number;i++){
   		ave_r0[i]=0.;
   		ave_r0_2[i]=0.;
   		ave_r1[i]=0.;
   		ave_r1_2[i]=0.;
   		prog0[i]=0.;
   		prog0_2[i]=0.;
   		prog1[i]=0.;
   		prog1_2[i]=0.;
   		s_r0[i]=0.;
   		s_r1[i]=0.;
   		for(unsigned int j=0;j<block_size;j++){
   			for(unsigned int j=0;j<3;j++){
				y_0[j]=rnd.Gauss(x_0[j], delta_0);
				y_1[j]=rnd.Gauss(x_1[j], delta_1);
			}	
	   		p0x=exp(-2.*sqrt(x_0[0]*x_0[0]+x_0[1]*x_0[1]+x_0[2]*x_0[2]))/pi;	   			
	   		p0y=exp(-2.*sqrt(y_0[0]*y_0[0]+y_0[1]*y_0[1]+y_0[2]*y_0[2]))/pi;
		   	a0=min(1.,(p0y/p0x));
		   	judge_0 = rnd.Rannyu();
		   	if(judge_0<=a0)
		   		for(unsigned int j=0;j<3;j++)
		   			x_0[j]=y_0[j];
		  	
		   	ave_r0[i]+=sqrt(x_0[0]*x_0[0]+x_0[1]*x_0[1]+x_0[2]*x_0[2]);
		   	
		   	p1x=exp(-sqrt(x_1[0]*x_1[0]+x_1[1]*x_1[1]+x_1[2]*x_1[2]))*x_1[2]*x_1[2]/(32.*pi);	   			
	   		p1y=exp(-sqrt(y_1[0]*y_1[0]+y_1[1]*y_1[1]+y_1[2]*y_1[2]))*y_1[2]*y_1[2]/(32.*pi);
		   	a1=min(1.,(p1y/p1x));
		   	judge_1 = rnd.Rannyu();
		   	if(judge_1<=a1)
		   		for(unsigned int j=0;j<3;j++)
		   			x_1[j]=y_1[j];
		   	ave_a1+=a1;
		   	ave_r1[i]+=sqrt(x_1[0]*x_1[0]+x_1[1]*x_1[1]+x_1[2]*x_1[2]);
		   		
	   	}
	   	
	   	ave_r0[i]/=double(block_size);
	   	ave_r0_2[i]=pow(ave_r0[i],2.);
	   	ave_r1[i]/=double(block_size);
	   	ave_r1_2[i]=pow(ave_r1[i],2.);
	   	
	   	for(unsigned int k=0;k<i+1;k++){
			prog0[i]+=ave_r0[k];
			prog1[i]+=ave_r1[k];
			prog0_2[i]+=ave_r0_2[k];
			prog1_2[i]+=ave_r1_2[k];
		}
		
		prog0[i]/=double(i+1);
		prog1[i]/=double(i+1);
		prog0_2[i]/=double(i+1);
		prog1_2[i]/=double(i+1);
		
		
		if(i==0){
			s_r0[i]=0.;
			s_r1[i]=0.;
		}
		
		else{
			s_r0[i]=sqrt((prog0_2[i]-pow(prog0[i],2))/double(i));
			s_r1[i]=sqrt((prog1_2[i]-pow(prog1[i],2))/double(i));
		}
		
		r_0 << i+1 << ";" << prog0[i] << ";" << s_r0[i] << endl;
		r_1 << i+1 << ";" << prog1[i] << ";" << s_r1[i] << endl;
   	}
   	
   	r_0.close();
   	r_1.close();
   	
	return 0;
}
