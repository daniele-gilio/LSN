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
  //Block Variables
  unsigned int M = 1E6;
  unsigned int block_number = 100;
  unsigned int block_size = M/block_number;
  unsigned int test = 10000;
  //Uniform Variables
  double *x_0 = new double [3];
  double *y_0 = new double [3];
  double *x_1 = new double [3];
  double *y_1 = new double [3];
  double ave_a0 = 0.;
  double ave_a1 = 0.;
  double delta_0 = 1.2;
  double delta_1 = 2.95;
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
  //Gaussian Variables
  double ave_a0_g = 0.;
  double ave_a1_g = 0.;
  double *x_0_g = new double [3];
  double *y_0_g = new double [3];
  double *x_1_g = new double [3];
  double *y_1_g = new double [3];
  double *ave_r0_g = new double [block_number];
  double *ave_r0_2_g = new double [block_number];
  double *ave_r1_g = new double [block_number];
  double *ave_r1_2_g = new double [block_number];
  double *prog0_g = new double [block_number];
  double *prog1_g = new double [block_number];
  double *prog0_2_g = new double [block_number];
  double *prog1_2_g = new double [block_number];
  double *s_r0_g = new double [block_number];
  double *s_r1_g = new double [block_number];
  //Common Variables
  double delta_0_g = 0.75;
  double delta_1_g = 1.85;
  double px = 0.;
  double py = 0.;
  double judge = 0.;
  double a = 0.;
	//Set starting points
 	x_0[0]=1.;  x_0_g[0]=1.;
 	x_0[1]=1.;  x_0_g[1]=1.;
 	x_0[2]=1.;  x_0_g[2]=1.;

 	x_1[0]=1.;  x_1_g[0]=1.;
 	x_1[1]=1.;  x_1_g[1]=1.;
 	x_1[2]=1.;  x_1_g[2]=1.;

 	ofstream r_0, r_1, r_0_g, r_1_g;
 	r_0.open("test_run_r0.dat");
 	r_1.open("test_run_r1.dat");
  r_0_g.open("test_run_r0_g.dat");
  r_1_g.open("test_run_r1_g.dat");

  //Test run to check for the 50% acceptance rate
	for(unsigned int i=0;i<test;i++){
		for(unsigned int j=0;j<3;j++){
			y_0[j]=rnd.Rannyu(x_0[j]-delta_0,x_0[j]+delta_0);
			y_1[j]=rnd.Rannyu(x_1[j]-delta_1,x_1[j]+delta_1);
      y_0_g[j]=rnd.Gauss(x_0_g[j],delta_0_g);
			y_1_g[j]=rnd.Gauss(x_1_g[j],delta_1_g);
		}
      //Uniform Metropolis
      px=exp(-2.*sqrt(x_0[0]*x_0[0]+x_0[1]*x_0[1]+x_0[2]*x_0[2]))/pi;
   		py=exp(-2.*sqrt(y_0[0]*y_0[0]+y_0[1]*y_0[1]+y_0[2]*y_0[2]))/pi;
	   	a=min(1.,(py/px));
	   	judge = rnd.Rannyu();
	   	if(judge<=a)
	   		for(unsigned int j=0;j<3;j++)
	   			x_0[j]=y_0[j];
	   	ave_a0+=a;

      //Gauss Metropolis
      px=exp(-2.*sqrt(x_0_g[0]*x_0_g[0]+x_0_g[1]*x_0_g[1]+x_0_g[2]*x_0_g[2]))/pi;
   		py=exp(-2.*sqrt(y_0_g[0]*y_0_g[0]+y_0_g[1]*y_0_g[1]+y_0_g[2]*y_0_g[2]))/pi;
	   	a=min(1.,(py/px));
	   	judge = rnd.Rannyu();
	   	if(judge<=a)
	   		for(unsigned int j=0;j<3;j++)
	   			x_0_g[j]=y_0_g[j];
	   	ave_a0_g+=a;
      //Output to file
	   	r_0 << x_0[0] << ";" << x_0[1] << ";" << x_0[2] << endl;
      r_0_g << x_0_g[0] << ";" << x_0_g[1] << ";" << x_0_g[2] << endl;
      //Uniform Metropolis
	   	px=exp(-sqrt(x_1[0]*x_1[0]+x_1[1]*x_1[1]+x_1[2]*x_1[2]))*x_1[2]*x_1[2]/(32.*pi);
   		py=exp(-sqrt(y_1[0]*y_1[0]+y_1[1]*y_1[1]+y_1[2]*y_1[2]))*y_1[2]*y_1[2]/(32.*pi);
	   	a=min(1.,(py/px));
	   	judge = rnd.Rannyu();
	   	if(judge<=a)
	   		for(unsigned int j=0;j<3;j++)
	   			x_1[j]=y_1[j];
	   	ave_a1+=a;
      //Gauss Metropolis
      px=exp(-sqrt(x_1_g[0]*x_1_g[0]+x_1_g[1]*x_1_g[1]+x_1_g[2]*x_1_g[2]))*x_1_g[2]*x_1_g[2]/(32.*pi);
   		py=exp(-sqrt(y_1_g[0]*y_1_g[0]+y_1_g[1]*y_1_g[1]+y_1_g[2]*y_1_g[2]))*y_1_g[2]*y_1_g[2]/(32.*pi);
	   	a=min(1.,(py/px));
	   	judge = rnd.Rannyu();
	   	if(judge<=a)
	   		for(unsigned int j=0;j<3;j++)
	   			x_1_g[j]=y_1_g[j];
	   	ave_a1_g+=a;
	   	r_1 << x_1[0] << ";" << x_1[1] << ";" << x_1[2] << endl;
      r_1_g << x_1_g[0] << ";" << x_1_g[1] << ";" << x_1_g[2] << endl;
	}

	r_0.close();
	r_1.close();
  r_0_g.close();
  r_1_g.close();
	ave_a0/=test/100.;
  ave_a0_g/=test/100.;
	ave_a1/=test/100.;
  ave_a1_g/=test/100.;

	cout << "The uniform acceptance rate for the s state is: " << ave_a0 << "%" << endl;
	cout << "The uniform acceptance rate for the p state is: " << ave_a1 << "%" << endl;
  cout << "The gaussian acceptance rate for the s state is: " << ave_a0_g << "%" << endl;
	cout << "The gaussian acceptance rate for the p state is: " << ave_a1_g << "%" << endl;

  r_0.open("r0.dat");
	r_1.open("r1.dat");
  r_0_g.open("r0_gauss.dat");
  r_1_g.open("r1_gauss.dat");

  //Reset starting points
 	x_0[0]=1.;  x_0_g[0]=1.;
 	x_0[1]=1.;  x_0_g[1]=1.;
 	x_0[2]=1.;  x_0_g[2]=1.;

 	x_1[0]=1.;  x_1_g[0]=1.;
 	x_1[1]=1.;  x_1_g[1]=1.;
 	x_1[2]=1.;  x_1_g[2]=1.;

  for(unsigned int i=0;i<block_number;i++){
  	for(unsigned int j=0;j<block_size;j++){
  		for(unsigned int j=0;j<3;j++){
        y_0[j]=rnd.Rannyu(x_0[j]-delta_0,x_0[j]+delta_0);
				y_1[j]=rnd.Rannyu(x_1[j]-delta_1,x_1[j]+delta_1);
        y_0_g[j]=rnd.Gauss(x_0_g[j],delta_0_g);
  			y_1_g[j]=rnd.Gauss(x_1_g[j],delta_1_g);
			}
        //Uniform Metropolis
	   		px=exp(-2.*sqrt(x_0[0]*x_0[0]+x_0[1]*x_0[1]+x_0[2]*x_0[2]))/pi;
	   		py=exp(-2.*sqrt(y_0[0]*y_0[0]+y_0[1]*y_0[1]+y_0[2]*y_0[2]))/pi;
		   	a=min(1.,(py/px));
		   	judge = rnd.Rannyu();
		   	if(judge<=a)
		   		for(unsigned int j=0;j<3;j++)
		   			x_0[j]=y_0[j];
        ave_r0[i]+=sqrt(x_0[0]*x_0[0]+x_0[1]*x_0[1]+x_0[2]*x_0[2]);

        //Gauss Metropolis
        px=exp(-2.*sqrt(x_0_g[0]*x_0_g[0]+x_0_g[1]*x_0_g[1]+x_0_g[2]*x_0_g[2]))/pi;
     		py=exp(-2.*sqrt(y_0_g[0]*y_0_g[0]+y_0_g[1]*y_0_g[1]+y_0_g[2]*y_0_g[2]))/pi;
  	   	a=min(1.,(py/px));
  	   	judge = rnd.Rannyu();
  	   	if(judge<=a)
  	   		for(unsigned int j=0;j<3;j++)
  	   			x_0_g[j]=y_0_g[j];
        ave_r0_g[i]+=sqrt(x_0[0]*x_0[0]+x_0[1]*x_0[1]+x_0[2]*x_0[2]);

        //Uniform Metropolis
		   	px=exp(-sqrt(x_1[0]*x_1[0]+x_1[1]*x_1[1]+x_1[2]*x_1[2]))*x_1[2]*x_1[2]/(32.*pi);
	   		py=exp(-sqrt(y_1[0]*y_1[0]+y_1[1]*y_1[1]+y_1[2]*y_1[2]))*y_1[2]*y_1[2]/(32.*pi);
		   	a=min(1.,(py/px));
		   	judge = rnd.Rannyu();
		   	if(judge<=a)
		   		for(unsigned int j=0;j<3;j++)
		   			x_1[j]=y_1[j];
		   	ave_r1[i]+=sqrt(x_1[0]*x_1[0]+x_1[1]*x_1[1]+x_1[2]*x_1[2]);

        //Gauss Metropolis
        px=exp(-sqrt(x_1_g[0]*x_1_g[0]+x_1_g[1]*x_1_g[1]+x_1_g[2]*x_1_g[2]))*x_1_g[2]*x_1_g[2]/(32.*pi);
     		py=exp(-sqrt(y_1_g[0]*y_1_g[0]+y_1_g[1]*y_1_g[1]+y_1_g[2]*y_1_g[2]))*y_1_g[2]*y_1_g[2]/(32.*pi);
  	   	a=min(1.,(py/px));
  	   	judge = rnd.Rannyu();
  	   	if(judge<=a)
  	   		for(unsigned int j=0;j<3;j++)
  	   			x_1_g[j]=y_1_g[j];
        ave_r1_g[i]+=sqrt(x_1[0]*x_1[0]+x_1[1]*x_1[1]+x_1[2]*x_1[2]);

	   	}

	   	ave_r0[i]/=double(block_size);
	   	ave_r0_2[i]=pow(ave_r0[i],2.);
	   	ave_r1[i]/=double(block_size);
	   	ave_r1_2[i]=pow(ave_r1[i],2.);
      ave_r0_g[i]/=double(block_size);
	   	ave_r0_2_g[i]=pow(ave_r0_g[i],2.);
	   	ave_r1_g[i]/=double(block_size);
	   	ave_r1_2_g[i]=pow(ave_r1_g[i],2.);

	   	for(unsigned int k=0;k<i+1;k++){
			prog0[i]+=ave_r0[k];
			prog1[i]+=ave_r1[k];
			prog0_2[i]+=ave_r0_2[k];
			prog1_2[i]+=ave_r1_2[k];
      prog0_g[i]+=ave_r0_g[k];
			prog1_g[i]+=ave_r1_g[k];
			prog0_2_g[i]+=ave_r0_2_g[k];
			prog1_2_g[i]+=ave_r1_2_g[k];
		}

		prog0[i]/=double(i+1);
		prog1[i]/=double(i+1);
		prog0_2[i]/=double(i+1);
		prog1_2[i]/=double(i+1);
    prog0_g[i]/=double(i+1);
    prog1_g[i]/=double(i+1);
    prog0_2_g[i]/=double(i+1);
    prog1_2_g[i]/=double(i+1);

		if(i==0){
			s_r0[i]=0.;
			s_r1[i]=0.;
      s_r0_g[i]=0.;
			s_r1_g[i]=0.;
		}

		else{
			s_r0[i]=sqrt((prog0_2[i]-pow(prog0[i],2))/double(i));
			s_r1[i]=sqrt((prog1_2[i]-pow(prog1[i],2))/double(i));
      s_r0_g[i]=sqrt((prog0_2_g[i]-pow(prog0_g[i],2))/double(i));
			s_r1_g[i]=sqrt((prog1_2_g[i]-pow(prog1_g[i],2))/double(i));
		}

		r_0 << i+1 << ";" << prog0[i] << ";" << s_r0[i] << endl;
		r_1 << i+1 << ";" << prog1[i] << ";" << s_r1[i] << endl;
    r_0_g << i+1 << ";" << prog0_g[i] << ";" << s_r0_g[i] << endl;
		r_1_g << i+1 << ";" << prog1_g[i] << ";" << s_r1_g[i] << endl;
    }

   	r_0.close();
   	r_1.close();
    r_0_g.close();
   	r_1_g.close();

	return 0;
}
