/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"
#include "random.h"

using namespace std;

int main(){
  	bool restart = true;
  	unsigned int counter = 0; 
  	while(restart==true){
  		Input(restart,counter);  //Inizialization
  		cout << "Do you want to rescale velocities? (1=yes,0=no): ";
		cin >> restart;
		counter++;
	}
 	int nconf = 1;
  	for(int istep=1; istep <= nstep; ++istep){
     		Move();           //Move particles with Verlet algorithm
     		if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     		if(istep%10 == 0){
        		Measure();     //Properties measurement
        		//ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        		nconf += 1;
     		}
  	}
  	ConfFinal();         //Write final configuration to restart

  	return 0;
}


void Input(bool restart, unsigned int counter){ //Prepare all stuff for the simulation
  	ifstream ReadInput,ReadConf;
  	double ep, ek, pr, et, vir;
  
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
  
	
	
	//Conglomerate couts to make only one control (output on terminal is a little prettier)
	if(counter==0){
		ReadInput.open("input.dat"); //Read input

		ReadInput >> temp;

		ReadInput >> npart;
	

		ReadInput >> rho;
	
		vol = (double)npart/rho;
	
		box = pow(vol,1.0/3.0);
	
		ReadInput >> rcut;
		ReadInput >> delta;
		ReadInput >> nstep;
		ReadInput >> iprint;
		cout << "Classic Lennard-Jones fluid        " << endl;
	  	cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
	 	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
	  	cout << "The program uses Lennard-Jones units " << endl;
	
		cout << "Number of particles = " << npart << endl;
		cout << "Density of particles = " << rho << endl;
		cout << "Volume of the simulation box = " << vol << endl;
		cout << "Edge of the simulation box = " << box << endl;
		cout << "The program integrates Newton equations with the Verlet method " << endl;
		cout << "Time step = " << delta << endl;
		cout << "Number of steps = " << nstep << endl << endl;
		ReadInput.close();
	}
	

	//Prepare array for measurements
	iv = 0; //Potential energy
	ik = 1; //Kinetic energy
	ie = 2; //Total energy
	it = 3; //Temperature
	n_props = 4; //Number of observables

	//Read initial configuration (r(t))
  	cout << "Read initial configuration from file config.0 " << endl << endl;
  	ReadConf.open("config.0");
  	for (int i=0; i<npart; ++i){
    		ReadConf >> x[i] >> y[i] >> z[i];
    		x[i] = x[i] * box;
    		y[i] = y[i] * box;
    		z[i] = z[i] * box;
  	}
  	ReadConf.close();
	
	double sumv[3] = {0.0, 0.0, 0.0};
	if(counter==0){
		//Prepare initial velocities (first time)
	   	cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
	   	
	   	for (int i=0; i<npart; ++i){
	     		vx[i] = rnd.Rannyu() - 0.5;
	     		vy[i] = rnd.Rannyu() - 0.5;
	     		vz[i] = rnd.Rannyu() - 0.5;

	     		sumv[0] += vx[i];
	     		sumv[1] += vy[i];
	     		sumv[2] += vz[i];
	   	}
	
	   	for (int idim=0; idim<3; ++idim) 
	   		sumv[idim] /= (double)npart;
	   	double sumv2 = 0.0, fs;
	   	for (int i=0; i<npart; ++i){
	     		vx[i] = vx[i] - sumv[0];
	     		vy[i] = vy[i] - sumv[1];
	     		vz[i] = vz[i] - sumv[2];

	     		sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	   	}
	   	sumv2 /= (double)npart;
		//Rescaling due to poor random sampling (can be avoided if we use maxwell-boltzmann distribution for velocities)
	   	fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
	   	for (int i=0; i<npart; ++i){
	     		vx[i] *= fs;
	     		vy[i] *= fs;
	     		vz[i] *= fs;
			//Used to start verlet algorithm (can be used to reach equilibrium , not very good per se)
	     		xold[i] = x[i] - vx[i] * delta;
	     		yold[i] = y[i] - vy[i] * delta;
	     		zold[i] = z[i] - vz[i] * delta;
	   	}
	}
   
   	//Restart implementation
   	if(restart==true and counter!=0){
   	   ifstream test;
   	   //ofstream new_old; 
	   unsigned int test_step = 500;
	   for(unsigned int i=0;i<test_step;i++){
	   	Move();
	   	if(i==test_step-2)
	   		ConfOld();//save r(t-dt)
	   }
	   Move();//arrive at r(t+dt)
	   test.open("config.old");
	   double x_t=0., y_t=0., z_t=0., k=0., fs=0.;
	   for(unsigned int i=0;i<npart;i++){
	   	test >> x_t >> y_t >> z_t;
	   	if(int(x_t)%1==0){
	   	vx[i]=(x[i]-x_t*box)/(2.*delta);
	   	vy[i]=(y[i]-y_t*box)/(2.*delta);
	   	vz[i]=(z[i]-z_t*box)/(2.*delta);
	   	k += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
	   	cout << k << endl;
	   	}
	   }
	   stima_temp = (2.0 / 3.0) * k/(double)npart; //Temperature
	   fs=temp/stima_temp;
	   test.close();
	   //new_old.open("config.old");
	   
	   /*for(unsigned int i=0;i<npart;i++){  
	   	//new_old << (x[i]-vx[i]*s*2*delta)/box << "	" << (y[i]-vy[i]*s*2*delta)/box << "	" << (z[i]-vz[i]*s*2*delta)/box << endl;
	   	
	   }
	   
	   //new_old.close();*/
	   for (int idim=0; idim<3; ++idim) 
	   		sumv[idim] /= (double)npart;
	   	double sumv2 = 0.0;
	   	for (int i=0; i<npart; ++i){
	     		vx[i] = vx[i] - sumv[0];
	     		vy[i] = vy[i] - sumv[1];
	     		vz[i] = vz[i] - sumv[2];

	     		sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	   	}
	   	sumv2 /= (double)npart;
	   	for (int i=0; i<npart; ++i){
	     		vx[i] *= fs;
	     		vy[i] *= fs;
	     		vz[i] *= fs;
			//Used to start verlet algorithm (can be used to reach equilibrium , not very good per se)
	     		xold[i] = (x[i]-vx[i]*2*delta)/box;
	     		yold[i] = (y[i]-vy[i]*2*delta)/box;
	     		zold[i] = (z[i]-vz[i]*2*delta)/box;
	   	}
	   
	   cout << "-----------------------------------------------------" << endl;
   	}
	   
   
   	return;
}


void Move(void){ 
	//Move particles with Verlet algorithm
  	double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  	for(int i=0; i<npart; ++i){ 
  		//Force acting on particle i
   		fx[i] = Force(i,0);
    		fy[i] = Force(i,1);
    		fz[i] = Force(i,2);
 	 }

  	for(int i=0; i<npart; ++i){ 
  		//Verlet integration scheme

    		xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    		ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    		znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    		vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    		vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    		vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

   		xold[i] = x[i];
    		yold[i] = y[i];
    		zold[i] = z[i];

    		x[i] = xnew;
    		y[i] = ynew;
    		z[i] = znew;
  	}
  	return;
}

double Force(int ip, int idir){ 
	//Compute forces as -Grad_ip V(r)
  	double f=0.0;
 	double dvec[3], dr;

  	for (int i=0; i<npart; ++i){
    		if(i != ip){
      			dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      			dvec[1] = Pbc( y[ip] - y[i] );
      			dvec[2] = Pbc( z[ip] - z[i] );

     			dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      			dr = sqrt(dr);
	
      			if(dr < rcut){
        			f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      			}	
      		}
      }
  
  	return f;
}

void Measure(){ 
	//Properties measurement
  	int bin;
  	double v, t, vij;
  	double dx, dy, dz, dr;
  	ofstream Epot, Ekin, Etot, Temp;

  	Epot.open("output_epot.dat",ios::app);
  	Ekin.open("output_ekin.dat",ios::app);
  	Temp.open("output_temp.dat",ios::app);
  	Etot.open("output_etot.dat",ios::app);

 	v = 0.0; //Reset observables
 	t = 0.0;

	//Cycle over pairs of particles
  	for (int i=0; i<npart-1; ++i){
    		for (int j=i+1; j<npart; ++j){

     			dx = Pbc( x[i] - x[j] );
     			dy = Pbc( y[i] - y[j] );
     			dz = Pbc( z[i] - z[j] );

     			dr = dx*dx + dy*dy + dz*dz;
     			dr = sqrt(dr);

     			if(dr < rcut){
       				vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

				//Potential energy
      				v += vij;
     			}
    		}          
    	}

	//Kinetic energy
  	for (int i=0; i<npart; ++i) 
  		t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
   	stima_pot = v/(double)npart; //Potential energy
    	stima_kin = t/(double)npart; //Kinetic energy
    	stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    	stima_etot = (t+v)/(double)npart; //Total energy

    	Epot << stima_pot  << endl;
    	Ekin << stima_kin  << endl;
    	Temp << stima_temp << endl;
    	Etot << stima_etot << endl;

   	Epot.close();
    	Ekin.close();
    	Temp.close();
    	Etot.close();

    	return;
}


void ConfFinal(void){ 
	//Write final configuration
  	ofstream WriteConf;

  	cout << "Print final configuration to file config.final " << endl << endl;
  	WriteConf.open("config.final");

  	for (int i=0; i<npart; ++i){
    		WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  	}	
  	WriteConf.close();
  	return;
}

void ConfOld(void){ 
	//Write final-1 configuration for rescaling
  	ofstream WriteConf;

  	cout << "Print second to last configuration to file config.old for rescaling purposes " << endl << endl;
  	WriteConf.open("config.old");

  	for (int i=0; i<npart; ++i){
    		WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  	}
  	WriteConf.close();
  	return;
}

void ConfXYZ(int nconf){ 
	//Write configuration in .xyz format
  	ofstream WriteXYZ;

  	WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  	WriteXYZ << npart << endl;
  	WriteXYZ << "This is only a comment!" << endl;
  	for (int i=0; i<npart; ++i){
    		WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  	}
 	 WriteXYZ.close();
}

double Pbc(double r){  
	//Algorithm for periodic boundary conditions with side L=box
    	return r - box * rint(r/box);
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
