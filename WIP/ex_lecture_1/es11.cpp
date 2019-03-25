#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

double error(double* ave, double* av2, int n){
	if(n==0)
		return 0;
	else
		return sqrt((av2[n]-pow(ave[n],2))/n);
}; 
 
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
   
	//Initialize Variables
	int M = 100000, block_number = 100;
	int block_size = M/block_number; //Take M as a multiple of block_number for simplicity and for this code to work
	int k=0; //Nested index
	double sum = 0;
	   
	//Initialize and fill arrays
	double *r = new double [M]; //Random array
	   
	for(int i=0;i<M;i++)
		r[i]=rnd.Rannyu();

	int *x = new int [block_number]; // Range array
	double *ave = new double [block_number]; // Average array
	double *av2 = new double [block_number]; // Average squared array
	double *prog_sum = new double [block_number]; // Progressive sum array
	double *prog_sum2 = new double [block_number]; // Progressive squared sum array
	double *prog_err = new double [block_number]; // Progressive error array
	   
	for(int i=0;i<block_number;i++){
		//Fill them fo safety
		ave[i]=0;
	   	av2[i]=0;
	   	prog_sum[i]=0;
	   	prog_sum2[i]=0;
	   	prog_err[i]=0;
	   	x[i]=i;
	}
	//End of inizialization
	   
	//Start of the actual program
	//First section   
	
	for(int i=0;i<block_number;i++){
	   	sum=0;
	   	for(int j=0;j<block_size;j++){
	   		k=j+(i*block_size);
	   		sum+=r[k];
	   	}
	   	ave[i]=sum/block_size;
	   	av2[i]=pow(ave[i],2);
	}
	   
	for(int i=0;i<block_number;i++){
	   	for(int j=0;j<i+1;j++){
	   		prog_sum[i] += ave[j];
	   		prog_sum2[i] += av2[j];
	   	}
	   	prog_sum[i]/=(i+1);
	   	prog_sum2[i]/=(i+1);
	   	prog_err[i] = error(prog_sum, prog_sum2, i);
	}
	   
	//First output to file   
	ofstream out;
	out.open("out1.dat");
	   
	for(int i=0;i<block_number;i++){
	   	out << x[i]*block_size << ";" << prog_sum[i] << ";" << prog_err[i] << endl;
	}
	   
	out.close();
	   
	//Second section
	//Let's start by cleaning up some stuff
	for(int i=0;i<block_number;i++){
	   	ave[i]=0;
	   	av2[i]=0;
	   	prog_sum[i]=0;
	   	prog_sum2[i]=0;
	   	prog_err[i]=0;
	   	x[i]=i;
	}
	   
	//Basically we do the same thing as before but we analyze the sigma
	
	for(int i=0;i<block_number;i++){
	   	sum=0;
	   	for(int j=0;j<block_size;j++){
	   		k=j+(i*block_size);
	   		sum+=pow(r[k]-0.5,2); //<-That's the main (and only) change
	   	}
	   	ave[i]=sum/block_size;
	   	av2[i]=pow(ave[i],2);
	}
	   
	for(int i=0;i<block_number;i++){
	   	for(int j=0;j<i+1;j++){
	   		prog_sum[i] += ave[j];
	   		prog_sum2[i] += av2[j];
	   	}
	   	prog_sum[i]/=(i+1);
	   	prog_sum2[i]/=(i+1);
	   	prog_err[i] = error(prog_sum, prog_sum2, i);
	}
	
	//Second output to file   
	   
	out.open("out2.dat");
	   
	for(int i=0;i<block_number;i++){
	   	out << x[i]*block_size << ";" << prog_sum[i] << ";" << prog_err[i] << endl;
	}
	   
	out.close();
	   
	//Third section
	
	//Initialize variables and arrays
	int chi_block_number = 100; //Divide [0,1] into chi_block_number intervals
	block_size=10000; //Reuse block_size 
	double expected = double(block_size)/chi_block_number; //It's worth defining it here since its constant
	double *r1 = new double [block_size];
	double *intervals = new double [chi_block_number+1];
	double *hits = new double [chi_block_number];
	double chi = 0;
	   
	for (int i=0;i<chi_block_number+1;i++)
	   	intervals[i]+=i*(1/double(chi_block_number));
	
	//Chi squared implementation
	out.open("out3.dat");
	//This section below is kinda bulky, further optimization may be required
	for(int l=0;l<chi_block_number;l++){
		//Reset everything before every chi computation
		chi=0;
		
		for(int i=0;i<block_size;i++)
		   	r1[i]=rnd.Rannyu();
		
		for(int i=0;i<chi_block_number;i++)
			hits[i]=0;
			
		//Calculate chi squared
		for(int i=0;i<chi_block_number;i++){
		   	for(int j=0;j<block_size;j++)
		   		if(r1[j] > intervals[i] and r1[j]<intervals[i+1]){
		   			hits[i]+=1;
		   		}
			chi+=pow(hits[i]-expected,2)/double(expected);
		}
		out << l+1 << ";" << chi << endl;		
	}
	
	out.close();
	
	return 0;
}

