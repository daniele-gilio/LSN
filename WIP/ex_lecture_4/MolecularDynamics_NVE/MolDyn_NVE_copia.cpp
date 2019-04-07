/*************************************************************************
**************************************************************************
    _/_/_/_/       _/_/_/_/_/       Numerical Simulation Laboratory
   _/     _/      _/               Physics Department
  _/     _/      _/  _/_/_/       Università degli Studi di Milano
 _/     _/ _    _/      _/  _    Daniele Gilio
_/_/_/_/  /_/  _/_/_/_/_/  /_/  email: daniele.gilio@studenti.unimi.it
**************************************************************************
**************************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"
#include "random.h"

using namespace std;

int main(){
  	unsigned int phase = 0; //0 for generic, 1 for solid, 2 for liquid, 3 for gas
  	//#pragma omp parallel for
  	for(phase=1;phase<4;phase++){
  		cout << endl << endl << "Phase = " << phase << endl << endl << endl;
	  	string generic = "generic";
	  	Execute(generic,phase);
	}
	
	//Inputs and outputs for real values conversion
	ofstream argon_k, argon_pot, argon_p, argon_temp, argon_etot;
	ofstream krypton_k, krypton_pot, krypton_p, krypton_temp, krypton_etot;
	ifstream kinetic, potential, pressure , temperature, total;
	
	double kin, pot, pres, tempe, tot;
	double s_kin, s_pot, s_pres, s_temp, s_tot;
	//Set real values for Argon and Krypton
	double s_a = 0.34, e_a = 120., m_a = 39.948, s_k = 0.364, e_k =164., m_k = 83.798, k_b = 1.38*1E-23;
	unsigned int n;
	
	for(unsigned int j=1;j<4;j++){
		//Open all files according to their phases
		kinetic.open("ave_kin_generic_"+to_string(j)+".dat");
		potential.open("ave_pot_generic_"+to_string(j)+".dat");
		pressure.open("ave_p_generic_"+to_string(j)+".dat");
		temperature.open("ave_temp_generic_"+to_string(j)+".dat");
		total.open("ave_etot_generic_"+to_string(j)+".dat");
		argon_k.open("ave_kin_argon_"+to_string(j)+".dat");
		argon_pot.open("ave_pot_argon_"+to_string(j)+".dat");
		argon_p.open("ave_p_argon_"+to_string(j)+".dat");
		argon_temp.open("ave_temp_argon_"+to_string(j)+".dat");
		argon_etot.open("ave_etot_argon_"+to_string(j)+".dat");
		krypton_k.open("ave_kin_krypton_"+to_string(j)+".dat");
		krypton_pot.open("ave_pot_krypton_"+to_string(j)+".dat");
		krypton_p.open("ave_p_krypton_"+to_string(j)+".dat");
		krypton_temp.open("ave_temp_krypton_"+to_string(j)+".dat");
		krypton_etot.open("ave_etot_krypton_"+to_string(j)+".dat");
		for(unsigned int i=0;i<block_number;i++){
			//Read Generic Average Values
			kinetic >> n >> kin >> s_kin;
			potential >> n >> pot >> s_pot;
			pressure >> n >> pres >> s_pres;
			temperature >> n >> tempe >> s_temp;
			total >> n >> tot >> s_tot;
			//Convert them to Argon (SI units)
			argon_k << n << ";" << kin*e_a*k_b << ";" << s_kin*e_a*k_b << endl;
			argon_pot << n  << ";" << pot*e_a*k_b << ";" << s_pot*e_a*k_b << endl;
			argon_p << n << ";" << pres*e_a*k_b*pow(s_a,-3.)*pow(1E9,3.)*1E-5 << ";" << s_pres*e_a*k_b*pow(s_a,-3.)*pow(1E9,3.)*1E-5 << endl;
			argon_temp << n << ";" << tempe*e_a << ";" <<  s_temp*e_a << endl;
			argon_etot << n << ";" << tot*e_a*k_b << ";" << s_tot*e_a*k_b << endl;
			//Convert them to Krypton (SI units)
			krypton_k << n << ";" << kin*e_k*k_b << ";" << s_kin*e_k*k_b << endl;
			krypton_pot << n << ";" << pot*e_k*k_b << ";" << s_pot*e_k*k_b << endl;
			krypton_p << n << ";" << pres*e_k*k_b*pow(s_k,-3.)*pow(1E9,3.)*1E-5 << ";" << s_pres*e_k*k_b*pow(s_k,-3.)*pow(1E9,3.)*1E-5 << endl;
			krypton_temp << n << ";" << tempe*e_k << ";" << s_temp*e_k << endl;
			krypton_etot << n << ";" << tot*e_k*k_b << ";" << s_tot*e_k*k_b << endl;
			
		}
		kinetic.close();
		potential.close();
		pressure.close();
		temperature.close();
		total.close();
		argon_k.close();
		argon_pot.close();
		argon_p.close();
		argon_temp.close();
		argon_etot.close();
		krypton_k.close();
		krypton_pot.close();
		krypton_p.close();
		krypton_temp.close();
		krypton_etot.close();
		
	}

  	return 0;
}

void Execute(string s, unsigned int phase){
	bool restart = true;
  	unsigned int counter = 0;
  	ofstream pot_ave, kin_ave, temp_ave, etot_ave,p_ave; 
  	while(restart==true){
  		Input(restart,counter,phase);  //Inizialization
  		if(ex==false){
  			cout << "Do you want to rescale velocities? (1=yes,0=no): ";
			cin >> restart;
			counter++;
		}
		else
			return;
	}
	
	unsigned int block_size = nstep/block_number;
	/*		  //-------------------------------------------------------------------------------------------------
	restart=false;   //Create a flag that makes the program start over as the first time with correct old configuration
	counter = 1984; //The Big Brother is watching you!
		       //----------------------------------------------------------------------------------------------------
	Input(restart, counter, phase);*/
 	int nconf = 1;
 	int h=1;
 	//Open generic average outputs
 	pot_ave.open("ave_pot_"+s+"_"+to_string(phase)+".dat");
 	kin_ave.open("ave_kin_"+s+"_"+to_string(phase)+".dat");
 	temp_ave.open("ave_temp_"+s+"_"+to_string(phase)+".dat");
 	etot_ave.open("ave_etot_"+s+"_"+to_string(phase)+".dat");
 	p_ave.open("ave_p_"+s+"_"+to_string(phase)+".dat");
 	
 	//Actual excecution of the simulation via data blocking
  	for(unsigned int i=1;i<block_number+1;i++){
  		for(unsigned int j=1;j<block_size+1;j++){
	     		Move();           //Move particles with Verlet algorithm
	     		if(h%iprint == 0) cout << "Number of time-steps: " << i*block_size << endl;
	     		Measure(i);     //Properties measurement
	     		if(h%10 == 0 or h==1){
				//ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
				//The line above can be improved to output different phases in different files 
				nconf += 1;
	     		}
	     		h++; //Global counter 
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
	     	//Compute averages
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
	     	
	     	//Compute errors
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
	     	
	     	//Output averages
	     	pot_ave << i << " " << prog_pot[i] << " " << sigma_pot[i] << endl;
	     	kin_ave << i << " " << prog_kin[i] << " " << sigma_kin[i] << endl;
	     	temp_ave << i << " " << prog_temp[i] << " " << sigma_temp[i] << endl;
	     	etot_ave << i << " " << prog_etot[i] << " " << sigma_etot[i] << endl;
	     	p_ave << i << " " << prog_p[i] << " " << sigma_p[i] << endl;
	     		
  	}
  	pot_ave.close();
  	kin_ave.close();
  	temp_ave.close();
  	etot_ave.close();
  	p_ave.close();
  	
  	ConfFinal(); //Save final configuration
}

void Input(bool restart, unsigned int counter, unsigned int phase){ //Prepare all stuff for the simulation
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
		if(phase==0){
			ReadInput.open("input.dat"); //Read generic input
		}
		
		if(phase==1){
			ReadInput.open("input.solid"); //Read Solid input
		}
		
		if(phase==2){
			ReadInput.open("input.liquid"); //Read Liquid input
		}
		
		if(phase==3){
			ReadInput.open("input.gas"); //Read Gas input
		}
		
		else if(phase>3){
			cerr << "Wrong Phase input!!" << endl; //Throw out an error if phase number is wrong
			ex = true;
			return;
		}

		//Read input phase
		ReadInput >> temp;
		ReadInput >> npart;
		ReadInput >> rho;
		ReadInput >> rcut;
		ReadInput >> delta;
		ReadInput >> nstep;
		ReadInput >> iprint;
		
		//Compute necessary dimensions
		vol = (double)npart/rho;
		box = pow(vol,1.0/3.0);
		
		//Spit out informations on the program itself and on phase configuration
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
	
	//First Start Setup
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
   	   ofstream new_old; 
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
	   new_old.open("config.old");
	   
	   for(unsigned int i=0;i<npart;i++){  
	   	new_old << (x[i]-vx[i]*s*2*delta)/box << "	" << (y[i]-vy[i]*s*2*delta)/box << "	" << (z[i]-vz[i]*s*2*delta)/box << endl;
	   	
	   }
	   
	   new_old.close();
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

void Measure(unsigned int i){ 
	//Properties measurement
  	int bin;
  	double v, t, vij, q;
  	double dx, dy, dz, dr;
  	ofstream Epot, Ekin, Etot, Temp;

  	Epot.open("output_epot.dat",ios::app);
  	Ekin.open("output_ekin.dat",ios::app);
  	Temp.open("output_temp.dat",ios::app);
  	Etot.open("output_etot.dat",ios::app);

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
     			//q is used to compute the virial
     			q+=pow(dr,-12.)-0.5*pow(dr,-6.);
    		}          
    	}

	//Kinetic energy
  	for (int i=0; i<npart; ++i) 
  		t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
   	//Potential energy and average
   	stima_pot = v/(double)npart; 
   	ave_pot[i]+=stima_pot;
   	//Kinetic energy and average
    	stima_kin = t/(double)npart; 
    	ave_kin[i]+=stima_kin;
    	//Temperature and average
    	stima_temp = (2.0 / 3.0) * t/(double)npart; 
    	ave_temp[i]+=stima_temp;
    	//Total energy 
    	stima_etot = (t+v)/(double)npart; 
    	ave_etot[i]+=stima_etot;
    	//Pressure and average
    	stima_p = rho*stima_temp+q*48./(3.*double(npart)*vol);
    	ave_p[i]+=stima_p;

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

double error(double* ave, double* av2, int n){
	if(n==0)
		return 0.;
	else
		return sqrt((av2[n]-pow(ave[n],2))/double(n));
}

//Original Code by:
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
