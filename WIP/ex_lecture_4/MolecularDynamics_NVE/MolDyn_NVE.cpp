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
	unsigned int phase = 0.; //0 for generic, 1 for solid, 2 for liquid, 3 for gas
  	char s[] = "sk"; //argon or krypton
  	Execute(s);

  	return 0;
}

void Execute(char []){
	bool restart = true;
  	unsigned int counter = 0;
  	ofstream pot_ave, kin_ave, temp_ave, etot_ave,p_ave; 
  	while(restart==true){
  		Input(restart,counter);  //Inizialization
  		cout << "Do you want to rescale velocities? (1=yes,0=no): ";
		cin >> restart;
		counter++;
	}
	unsigned int block_size = nstep/block_number;
	restart=false;
	counter = 1984;
	Input(restart, counter);
 	int nconf = 1;
 	int h=1;
 	pot_ave.open("ave_pot.dat");
 	kin_ave.open("ave_kin.dat");
 	temp_ave.open("ave_temp.dat");
 	etot_ave.open("ave_etot.dat");
 	p_ave.open("ave_p.dat");
  	for(unsigned int i=1;i<block_number+1;i++){
  		for(unsigned int j=1;j<block_size+1;j++){
	     		Move();           //Move particles with Verlet algorithm
	     		if(h%iprint == 0) cout << "Number of time-steps: " << i*block_size << endl;
	     		Measure(i);     //Properties measurement
	     		if(h%10 == 0 or h==1){
				ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
				nconf += 1;
	     		}
	     		h++;
	     	}
	     	//Save averages for statistical analysis
	     	ave_pot[i]/=double(block_size);
	     	ave_pot_2[i]=pow(ave_pot[i],2.);
	     	ave_kin[i]/=double(block_size);
	     	ave_kin_2[i]=pow(ave_kin[i],2.);
	     	ave_temp[i]/=double(block_size);
	     	ave_temp_2[i]=pow(ave_temp[i],2.);
	     	ave_etot[i]/=double(block_size);
	     	ave_etot_2[i]=pow(ave_etot[i],2.);
	     	ave_p[i]/=double(block_size);
	     	ave_p_2[i]=pow(ave_p[i],2.);
	     	for(unsigned int j=0;j<i;j++){
	     		prog_pot[i]+=ave_pot[j];
	     		prog_kin[i]+=ave_kin[j];
	     		prog_temp[i]+=ave_temp[j];
	     		prog_etot[i]+=ave_etot[j];
	     		prog_p[i]+=ave_p[j];
	     		prog_pot_2[i]+=ave_pot_2[j];
	     		prog_kin_2[i]+=ave_kin_2[j];
	     		prog_temp_2[i]+=ave_temp_2[j];
	     		prog_etot_2[i]+=ave_etot_2[j];
	     		prog_p_2[i]+=ave_p_2[j];
	     	}
	     	prog_pot[i]/=double(i);
	     	prog_kin[i]/=double(i);
	     	prog_temp[i]/=double(i);
	     	prog_etot[i]/=double(i);
	     	prog_p[i]/=double(i);
	     	prog_pot_2[i]/=double(i);
	     	prog_kin_2[i]/=double(i);
	     	prog_temp_2[i]/=double(i);
	     	prog_etot_2[i]/=double(i);
	     	prog_p_2[i]/=double(i);
	     	if(i==1){
	     		sigma_pot[i]=0.;
	     		sigma_kin[i]=0.;
	     		sigma_temp[i]=0.;
	     		sigma_etot[i]=0.;
	     		sigma_p[i]=0.;
	     	}
	     	else{
	     		sigma_pot[i]=error(prog_pot,prog_pot_2,i-1);
	     		sigma_kin[i]=error(prog_kin,prog_kin_2,i-1);
	     		sigma_temp[i]=error(prog_temp,prog_temp_2,i-1);
	     		sigma_etot[i]=error(prog_etot,prog_etot_2,i-1);
	     		sigma_p[i]=error(prog_p,prog_p_2,i-1);
	     	}
	     	
	     	pot_ave << i << ";" << prog_pot[i] << ";" << sigma_pot[i] << endl;
	     	kin_ave << i << ";" << prog_kin[i] << ";" << sigma_kin[i] << endl;
	     	temp_ave << i << ";" << prog_temp[i] << ";" << sigma_temp[i] << endl;
	     	etot_ave << i << ";" << prog_etot[i] << ";" << sigma_etot[i] << endl;
	     	p_ave << i << ";" << prog_p[i] << ";" << sigma_p[i] << endl;
	     		
  	}
  	pot_ave.close();
  	kin_ave.close();
  	temp_ave.close();
  	etot_ave.close();
  	p_ave.close();
  	ConfFinal();         //Write final configuration to restart
}

void Input(bool restart, unsigned int counter){ //Prepare all stuff for the simulation
  	ifstream ReadInput,ReadConf;
  	double ep, ek, pr, et, vir,fs;
  
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
  
	
	
	//Conglomerate to make only one control (output on terminal is a little prettier)
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
  		//Read r(t)
    		ReadConf >> x[i] >> y[i] >> z[i];
    		x[i] = x[i] * box;
    		x_t[i] = x[i];
    		y[i] = y[i] * box;
    		y_t[i] = y[i];
    		z[i] = z[i] * box;
    		z_t[i] = z[i];
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
	   	double sumv2 = 0.0;
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
	     		xold[i] = Pbc(x[i] - vx[i] * delta);
	     		xold_t[i]=xold[i];
	     		yold[i] = Pbc(y[i] - vy[i] * delta);
	     		yold_t[i] = yold[i];
	     		zold[i] = Pbc(z[i] - vz[i] * delta);
	     		zold_t[i] = zold[i];
	   	}
	   	
	   	ConfOld(); //Save first old configuration
	}
	
	if(restart==true and counter!=0){
		ifstream old;
		old.open("config.old"); 
		for(unsigned int i=0;i<npart;i++){
			//Read r(t-dt)
			old >> xold_t[i] >> yold_t[i] >> zold_t[i];
			xold_t[i]*=box;
			yold_t[i]*=box;
			zold_t[i]*=box;
			xold[i]=xold_t[i];
			yold[i]=yold_t[i];
			zold[i]=zold_t[i];
			if(i==0)
				Move();//arrive at r(t+dt) and compute first velocities
			//Save r(t+dt) for later use
			x_t[i]=x[i];
			y_t[i]=y[i];
			z_t[i]=z[i];
			//Save first velocities for later use
			vx_t[i] = vx[i];
    			vy_t[i] = vy[i];
    			vz_t[i] = vz[i];
    			
    			cout << vx_t[i] << " " << vy_t[i] << " " << vz_t[i] << endl;
		}
		old.close();
		unsigned int test_steps = 500;
		double k = 0.;
		
		//Make a test simulation to get the "real" temperature
		for(unsigned int i=0;i<test_steps;i++)
			Move();
		
		//Compute "real" temperature	
		for (unsigned int i=0; i<npart; ++i) 
  			{k += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);cout << vx[i] << endl;}
  		
  		stima_temp = (2.0 / 3.0) * k/(double)npart; //Temperature
  		//cout << stima_temp << endl;
		
		fs=temp/stima_temp; //Rescale factor for velocities
		//cout << fs << endl;
		
		for(unsigned int i=0;i<npart;i++){
			//Rescale velocities
			vx_t[i]*=fs;
			vy_t[i]*=fs;
			vz_t[i]*=fs;
			//Calculate new r(t-dt)
			xold_t[i]=Pbc(x_t[i]-2.*delta*vx_t[i]);
			yold_t[i]=Pbc(y_t[i]-2.*delta*vy_t[i]);
			zold_t[i]=Pbc(z_t[i]-2.*delta*vz_t[i]);
		}
		
		ConfOld(); //Save new old configuration
	}
	
	if(restart==false and counter==1984){
		ifstream old;
		old.open("config.old");
		Move(); //arrive at r(t+dt)
		for(unsigned int i=0;i<npart;i++){
			//Read r(t-dt)
			old >> xold_t[i] >> yold_t[i] >> zold_t[i];
			xold_t[i]*=box;
			yold_t[i]*=box;
			zold_t[i]*=box;
			//Compute first velocities
			vx[i] = Pbc(x[i] - xold_t[i])/(2.0 * delta);
    			vy[i] = Pbc(y[i] - yold_t[i])/(2.0 * delta);
    			vz[i] = Pbc(z[i] - zold_t[i])/(2.0 * delta);
		}
		old.close();
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

void Measure(unsigned int i){ 
	//Properties measurement
  	int bin;
  	double v, t, vij, q, modulus;
  	double dx, dy, dz, dr;
  	ofstream Epot, Ekin, Etot, Temp, P;

  	Epot.open("output_epot.dat",ios::app);
  	Ekin.open("output_ekin.dat",ios::app);
  	Temp.open("output_temp.dat",ios::app);
  	Etot.open("output_etot.dat",ios::app);
  	P.open("output_p.dat",ios::app);

 	v = 0.0; //Reset observables
 	t = 0.0;
 	q=0.;
 	

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
     			
     			q+=pow(dr,-12.)-0.5*pow(dr,-6.);
    		}          
    	}

	//Kinetic energy
  	for (int i=0; i<npart; ++i) 
  		t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  		
   
   	stima_pot = v/(double)npart; //Potential energy
   	ave_pot[i]+=stima_pot;
    	stima_kin = t/(double)npart; //Kinetic energy
    	ave_kin[i]+=stima_kin;
    	stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    	ave_temp[i]+=stima_temp; 
    	stima_etot = (t+v)/(double)npart; //Total energy
    	ave_etot[i]+=stima_etot;
    	
    	stima_p = rho*stima_temp+q*48./(3.*double(npart)*vol);
    	ave_p[i]+=stima_p;
    	//cout << stima_p << endl;

    	Epot << stima_pot  << endl;
    	Ekin << stima_kin  << endl;
    	Temp << stima_temp << endl;
    	Etot << stima_etot << endl;
    	P << stima_p << endl;

   	Epot.close();
    	Ekin.close();
    	Temp.close();
    	Etot.close();
    	P.close();

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
    		WriteConf << xold_t[i]/box << "   " <<  yold_t[i]/box << "   " << zold_t[i]/box << endl;
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

double error(double* ave, double* av2, int n){
	if(n==0)
		return 0.;
	else
		return sqrt((av2[n]-pow(ave[n],2))/double(n));
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
