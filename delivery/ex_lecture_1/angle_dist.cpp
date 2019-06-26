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
	
	//Start of the actual program

	double point_1 = 0, point_2 = 0, point_3 = 0, point_4 = 0;
	unsigned int block_number = 100;
	unsigned int M = 1E6;
	unsigned int block_size = M/block_number; //as always take M as a multiple of block_number
	
	ofstream standard;
	ofstream circle;
	standard.open("angle_dist_st.dat");
	circle.open("angle_dist_circle.dat");
	
	for(unsigned int j=0;j<block_size;j++){
		//Circle Points
		point_1=rnd.Rannyu();
		point_2=rnd.Rannyu();
		//Standard Points
		point_3=rnd.Rannyu();
		point_4=rnd.Rannyu();
			
		while(pow(point_1,2)+pow(point_2,2)>1){ //keeping the points inside the unit circle (make a distribution graph to explain!!)
			point_1=rnd.Rannyu();
			point_2=rnd.Rannyu();
		}
		
		circle << atan(point_1/point_2) << endl;	
		standard << atan(point_3/point_4) << endl;	
	}	
	
	standard.close();
	circle.close();
	return 0;
}
