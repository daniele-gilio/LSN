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
   	
   	unsigned int M = 1E4;
   	unsigned int block_number = 100;
   	unsigned int block_size = M/block_number;
   	double S_0 = 100.;
   	double T = 1.;
	double K = 100.;
	double r = 0.1;
	double volatility = 0.25;
	double *call_d = new double [block_number]; //call array for direct calculation
	double *call_s = new double [block_number]; //call array for discretized calculation
	double *call_2_d = new double [block_number]; //call squared array (direct)
	double *call_2_s = new double [block_number]; //call squared array (sampled)
	double *put_d = new double [block_number]; //put array for direct calculation
	double *put_s = new double [block_number]; //put array for sampled calculation
	double *put_2_d = new double [block_number]; //put squared array (direct)
	double *put_2_s = new double [block_number]; //put squared array (sampled)
	double *p_sum_d = new double [block_number]; //put progressive sum for direct calculation 
	double *p_sum_s = new double [block_number]; //put progressive sum for sampled calculation
	double *p_sum_d_2 = new double [block_number]; //put progressive squared sum (direct)
	double *p_sum_s_2 = new double [block_number]; //put progressive squared sum (sampled)
	double *c_sum_d = new double [block_number]; //call progressive sum for direct calculation
	double *c_sum_s = new double [block_number]; //call progressive sum for sampled calculation
	double *c_sum_d_2 = new double [block_number]; //call progressive squared sum (direct)
	double *c_sum_s_2 = new double [block_number]; //call progressive squared sum (sampled)
	double *s_p_d = new double [block_number]; //direct put sigma
	double *s_p_s = new double [block_number]; //sampled put sigma
	double *s_c_d = new double [block_number]; //direct call sigma
	double *s_c_s = new double [block_number]; //sampled call sigma
	double S_t_d = 0.; //direct sampled value
	double *S_t_i = new double [100];  //sampled price array
	S_t_i[0]=S_0; //<- Put the 0th value at S_0
	unsigned int time_size = 100;
	double t = T/double(time_size); //Time intervals lenght
	
	
	ofstream direct;
	ofstream sample;
	direct.open("direct.dat");
	sample.open("sample.dat");
	
	for(unsigned int i=0;i<block_number;i++){
		call_d[i]=0.;
		put_d[i]=0.;
		call_2_d[i]=0.;
		put_2_d[i]=0.;
		p_sum_d[i]=0.;
		p_sum_d_2[i]=0;
		c_sum_d[i]=0.;
		c_sum_d_2[i]=0.;
		s_p_d[i]=0.;
		s_c_d[i]=0.;
		call_s[i]=0.;
		put_s[i]=0.;
		call_2_s[i]=0.;
		put_2_s[i]=0.;
		p_sum_s[i]=0.;
		p_sum_s_2[i]=0;
		c_sum_s[i]=0.;
		c_sum_s_2[i]=0.;
		s_p_s[i]=0.;
		s_c_s[i]=0.;
		for(unsigned int j=0;j<block_size;j++){
			//Sample prices every t
			for(unsigned int l=0;l<time_size;l++){
				S_t_i[l+1]=0.;
				S_t_i[l+1]=S_t_i[l]*exp((r-0.5*pow(volatility,2.))*t+volatility*rnd.Gauss(0.,T)*sqrt(t));
			}
			
			if(S_t_i[time_size]-K<0.)
				put_s[i]+=exp(-r*T)*(K-S_t_i[time_size]);
			else	
				call_s[i]+=exp(-r*T)*(S_t_i[time_size]-K);
				
			//Directly calculate the final price
			S_t_d=0.;
			S_t_d=S_0*exp((r-0.5*pow(volatility,2.))*T+volatility*rnd.Gauss(0.,T));
			
			if(S_t_d-K<0.)
				put_d[i]+=exp(-r*T)*(K-S_t_d);
			else 
				call_d[i]+=exp(-r*T)*(S_t_d-K);
			
		}
		
		//Statistical Analysis 
		put_d[i]/=double(block_size);
		put_2_d[i]=pow(put_d[i],2.);
		call_d[i]/=double(block_size);
		call_2_d[i]=pow(call_d[i],2.);
		
		put_s[i]/=double(block_size);
		put_2_s[i]=pow(put_s[i],2.);
		call_s[i]/=double(block_size);
		call_2_s[i]=pow(call_s[i],2.);
		
		for(unsigned int k=0;k<i+1;k++){
			p_sum_d[i]+=put_d[k];
			p_sum_d_2[i]+=put_2_d[k];
			c_sum_d[i]+=call_d[k];
			c_sum_d_2[i]+=call_2_d[k];
			p_sum_s[i]+=put_s[k];
			p_sum_s_2[i]+=put_2_s[k];
			c_sum_s[i]+=call_s[k];
			c_sum_s_2[i]+=call_2_s[k];
		}
		
		p_sum_d[i]/=double(i+1);
		p_sum_d_2[i]/=double(i+1);
		c_sum_d[i]/=double(i+1);
		c_sum_d_2[i]/=double(i+1);
		p_sum_s[i]/=double(i+1);
		p_sum_s_2[i]/=double(i+1);
		c_sum_s[i]/=double(i+1);
		c_sum_s_2[i]/=double(i+1);
		
		if(i==0){
			s_p_d[i]=0.;
			s_c_d[i]=0.;
			s_p_s[i]=0.;
			s_c_s[i]=0.;
		}
		
		else{
			s_p_d[i]=sqrt((p_sum_d_2[i]-pow(p_sum_d[i],2))/double(i));
			s_c_d[i]=sqrt((c_sum_d_2[i]-pow(c_sum_d[i],2))/double(i));
			s_p_s[i]=sqrt((p_sum_s_2[i]-pow(p_sum_s[i],2))/double(i));
			s_c_s[i]=sqrt((c_sum_s_2[i]-pow(c_sum_s[i],2))/double(i));
		}
		direct << i+1 << ";" << p_sum_d[i] << ";" << s_p_d[i] << ";" << c_sum_d[i] << ";" << s_c_d[i] << endl; //column 0 is shared than we have columns 1 and 2 for puts and columns 3 and 4 for calls
		sample << i+1 << ";" << c_sum_s[i] << ";" << s_c_s[i] << ";" << p_sum_s[i] << ";" << s_p_s[i] << endl; //column 0 is shared, columns 1 and 2 for calls and 3 and 4 for puts 
		
	}
	
	direct.close();
	sample.close();
	return 0;
}
