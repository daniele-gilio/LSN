#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

bool locate (double point_1, double point_2, double lenght, double d){
	bool k = false;
	//Calculate the angle of the stick
	double angle = atan(point_1/point_2);
	//Verify if it intersects a line
	if(d<=lenght*sin(angle)/2.)
		k=true;
	
	return k; 	
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
		while (!input.eof()){
			input >> property;
			if(property == "RANDOMSEED"){
			input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
			rnd.SetRandom(seed,p1,p2);
			}
		}
	input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
	rnd.SaveSeed();
	
	//Inizialize variables and arrays

	double point_1 = 0.;
	double point_2 = 0.;
	unsigned int block_number = 100;
	unsigned int M = 5*1E6;
	unsigned int block_size = M/block_number; //As always take M as a multiple of block_number
	unsigned int ruler = 100; //Number of (0,1) divisions
	double distance = 1./ruler; //Actual distance between grid lines
	double lenght = distance; //Needle lenght
	double d = 0.;  //Distance from center stick point to the nearest line
	double *pi = new double [block_number]; //Pi result's array
	int hits = 0;
	
	for(unsigned int i=0;i<block_number;i++){
		//Fill it for safety
		pi[i]=0;
	} 
	
	//Start of the actual program
	
	for(unsigned int i=0;i<block_number;i++){
		//Generate the needle position and check if it intersects the grid
		hits=0;
		for(unsigned int j=0;j<block_size;j++){
			point_1=rnd.Rannyu();
			point_2=rnd.Rannyu();
			
			while(pow(point_1,2)+pow(point_2,2)>1){ //<-Keeping the points inside the unit circle
				point_1=rnd.Rannyu();
				point_2=rnd.Rannyu();
			}
			
			d = rnd.Rannyu(0, distance/2.); //Take d (distance of the center of the needle from the nearest grid line) randomly from 0 to distance/2
			
			if(locate(point_1,point_2,lenght,d)==true)
				hits++;
		}
		
		pi[i]=2.*lenght*block_size/(hits*distance);
	}
	
	//Initialize variables, arrays and output for statistical analysis
	
	double prog_mean = 0;
	double prog_mean_2 = 0;
	double *mean = new double [block_number];
	double *sigma = new double [block_number];
	
	ofstream out;
	out.open("pi.dat");
	
	//Statistical analysis

	for(unsigned int i=0;i<block_number;i++){
		prog_mean+=pi[i];
		prog_mean_2+=pow(pi[i],2);
		if(i==0)
			sigma[i]=0; //<-Can't compute the first sigma due to the n-1 at the denominator, we put it at 0 for convenience
		else
			sigma[i]=sqrt((prog_mean_2/(i+1)-pow(prog_mean/(i+1),2))/i);
		mean[i]=prog_mean/(i+1);
		out << i << ";" << mean[i] << ";" << sigma[i] << endl;
	}
	
	out.close();	  
	
	
	return 0;
}
