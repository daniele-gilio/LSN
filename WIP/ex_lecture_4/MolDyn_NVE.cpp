/*************************************************************************
**************************************************************************
    _/_/_/_/       _/_/_/_/_/       Numerical Simulation Laboratory
   _/     _/      _/               Physics Department
  _/     _/      _/  _/_/_/       Università degli Studi di Milano
 _/     _/ _    _/      _/  _    Daniele Gilio
_/_/_/_/  /_/  _/_/_/_/_/  /_/  email: daniele.gilio@studenti.unimi.it
**************************************************************************
**************************************************************************/
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"
#include "random.h"

using namespace std;

int main(){
  Input();             //Inizialization
  int nconf = 1;
  h=1;
  block_size=nstep/block_number;
  while(c<block_number){
    for(int j=0; j<block_size; j++){
       Move();           //Move particles with Verlet algorithm
       if(h%iprint == 0) cout << "Number of time-steps: " << h << endl;
       if(h%10 == 0){
          Measure();     //Properties measurement
  //        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
          nconf += 1;
       }
       if(h==nstep-1)
        ConfOld(); //Write second to last configuration r(t-dt)
       h++;
      }
      if(real_measure==1){
        DataBlock();
      }
      c++;

    }
    if(phase!="generic"){
      Argon();
      Krypton();
    }
    ConfFinal();         //Write final configuration to restart r(t)

  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

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

  ReadInput.open("input.dat"); //Read input

  ReadInput >> restart;
  ReadInput >> real_measure;
  ReadInput >> phase;
  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl;
  cout << "Phase: " << phase << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

  if(restart==false){
  //Read initial configuration from config.0
    cout << "Read initial configuration from file config.0 " << endl << endl;
    ReadConf.open("config.0");
    for (int i=0; i<npart; ++i){
      ReadConf >> x[i] >> y[i] >> z[i];
      x[i] = x[i] * box;
      y[i] = y[i] * box;
      z[i] = z[i] * box;
    }
    ReadConf.close();

  //Prepare initial velocities
     cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
     double sumv[3] = {0.0, 0.0, 0.0};
     for (int i=0; i<npart; ++i){
       vx[i] = rnd.Rannyu() - 0.5;
       vy[i] = rnd.Rannyu() - 0.5;
       vz[i] = rnd.Rannyu() - 0.5;

       sumv[0] += vx[i];
       sumv[1] += vy[i];
       sumv[2] += vz[i];
     }
     for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
     double sumv2 = 0.0, fs;
     for (int i=0; i<npart; ++i){
       vx[i] = vx[i] - sumv[0];
       vy[i] = vy[i] - sumv[1];
       vz[i] = vz[i] - sumv[2];

       sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
     }
     sumv2 /= (double)npart;

     fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
     for (int i=0; i<npart; ++i){
       vx[i] *= fs;
       vy[i] *= fs;
       vz[i] *= fs;

       xold[i] = Pbc(x[i] - vx[i] * delta);
       yold[i] = Pbc(y[i] - vy[i] * delta);
       zold[i] = Pbc(z[i] - vz[i] * delta);
     }
   }

   else{
      //Read initial configuration from old.final r(t)
       cout << "Read initial configuration from file old.final " << endl << endl;
       ReadConf.open("old_"+phase+".final");
       for (int i=0; i<npart; ++i){
         ReadConf >> x[i] >> y[i] >> z[i];
         x[i] = x[i] * box;
         y[i] = y[i] * box;
         z[i] = z[i] * box;
       }
       ReadConf.close();
       //Read second to last configuration from old.0 r(t-dt)
       cout << "Read second to last configuration from file old.0 " << endl << endl;
       ReadConf.open("old_"+phase+".0");
       for (int i=0; i<npart; ++i){
         ReadConf >> xold[i] >> yold[i] >> zold[i];
         xold[i] = xold[i] * box;
         yold[i] = yold[i] * box;
         zold[i] = zold[i] * box;
       }
       ReadConf.close();

       Move(); //Arrive at r(t+dt) and compute velocities

       double sumv[3] = {0.0, 0.0, 0.0};
       for (int i=0; i<npart; ++i){
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
        fs = sqrt(3. * temp / sumv2);   // fs = velocity scale factor
        cout << "Velocity Scale Factor = " << fs << endl << endl;
        for (int i=0; i<npart; ++i){
          vx[i] *= fs;
          vy[i] *= fs;
          vz[i] *= fs;

          xold[i] = Pbc(x[i] - vx[i] * delta);
          yold[i] = Pbc(y[i] - vy[i] * delta);
          zold[i] = Pbc(z[i] - vz[i] * delta);
        }
   }
   return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

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

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
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

void Measure(){ //Properties measurement
  double v, t, vij, q;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Press;

  Epot.open("output_epot_"+phase+".dat",ios::app);
  Ekin.open("output_ekin_"+phase+".dat",ios::app);
  Temp.open("output_temp_"+phase+".dat",ios::app);
  Etot.open("output_etot_"+phase+".dat",ios::app);
  Press.open("output_press_"+phase+".dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;
  q=0.;

//cycle over pairs of particles
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
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    q/=double(npart);
    stima_pot = v/(double)npart; //Potential energy
    stima_kin = t/(double)npart; //Kinetic energy
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy
    stima_press = rho*stima_temp+(48.*q/(3.*vol));
    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Press << stima_press << endl;

    if(real_measure==true){
      ave_pot[c]+=stima_pot;
      ave_kin[c]+=stima_kin;
      ave_t[c]+=stima_temp;
      ave_tot[c]+=stima_etot;
      ave_press[c]+=stima_press;
    }

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Press.close();

    return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final and old.final" << endl << endl;
  WriteConf.open("config_"+phase+".final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

  WriteConf.open("old_"+phase+".final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfOld(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print second to last configuration to file old.0 " << endl << endl;
  WriteConf.open("old_"+phase+".0");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

void DataBlock(){
  ofstream Epot, Ekin, Etot, Temp, Press;

  Epot.open("ave_epot_"+phase+".dat",ios::app);
  Ekin.open("ave_ekin_"+phase+".dat",ios::app);
  Temp.open("ave_temp_"+phase+".dat",ios::app);
  Etot.open("ave_etot_"+phase+".dat",ios::app);
  Press.open("ave_press_"+phase+".dat",ios::app);

  ave_pot[c]/=double(block_size)/10.; //the 10 is there because we're measuring once every ten times
  ave_pot_2[c]=pow(ave_pot[c],2.);
  ave_kin[c]/=double(block_size)/10.;
  ave_kin_2[c]=pow(ave_kin[c],2.);
  ave_tot[c]/=double(block_size)/10.;
  ave_tot_2[c]=pow(ave_tot[c],2.);
  ave_t[c]/=double(block_size)/10.;
  ave_t_2[c]=pow(ave_t[c],2.);
  ave_press[c]/=double(block_size)/10.;
  ave_press_2[c]=pow(ave_press[c],2.);
  for(int k=0;k<c+1;k++){
    prog_pot[c]+=ave_pot[k];
    prog_pot_2[c]+=ave_pot_2[k];
    prog_kin[c]+=ave_kin[k];
    prog_kin_2[c]+=ave_kin_2[k];
    prog_tot[c]+=ave_tot[k];
    prog_tot_2[c]+=ave_tot_2[k];
    prog_t[c]+=ave_t[k];
    prog_t_2[c]+=ave_t_2[k];
    prog_press[c]+=ave_press[k];
    prog_press_2[c]+=ave_press_2[k];
  }
  prog_pot[c]/=double(c+1);
  prog_pot_2[c]/=double(c+1);
  prog_kin[c]/=double(c+1);
  prog_kin_2[c]/=double(c+1);
  prog_tot[c]/=double(c+1);
  prog_tot_2[c]/=double(c+1);
  prog_t[c]/=double(c+1);
  prog_t_2[c]/=double(c+1);
  prog_press[c]/=double(c+1);
  prog_press_2[c]/=double(c+1);

  if(c==0){
    s_pot[c]=0.;
    s_kin[c]=0.;
    s_tot[c]=0.;
    s_t[c]=0.;
    s_press[c]=0.;
  }

  else{
    s_pot[c]=sqrt((prog_pot_2[c]-pow(prog_pot[c],2))/double(c));
    s_kin[c]=sqrt((prog_kin_2[c]-pow(prog_kin[c],2))/double(c));
    s_tot[c]=sqrt((prog_tot_2[c]-pow(prog_tot[c],2))/double(c));
    s_t[c]=sqrt((prog_t_2[c]-pow(prog_t[c],2))/double(c));
    s_press[c]=sqrt((prog_press_2[c]-pow(prog_press[c],2))/double(c));
  }

  Epot << c+1 << " " << prog_pot[c] << " " << s_pot[c] << endl;
  Ekin << c+1 << " " << prog_kin[c] << " " << s_kin[c] << endl;
  Etot << c+1 << " " << prog_tot[c] << " " << s_tot[c] << endl;
  Temp << c+1 << " " << prog_t[c] << " " << s_t[c] << endl;
  Press << c+1 << " " << prog_press[c] << " " << s_press[c] << endl;

  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
  Press.close();

  return;
}

void Argon(){
  double sigma=0.34*1E-9;
  double epsilon=120.*k_be;
  double t = 120.;
  double p = t*k_bj*1E-5/pow(sigma,3.);
  double in;
  double s;
  int m;
  ifstream pot,kin,tot,pres,temp;
  pot.open("ave_epot_"+phase+".dat");
  kin.open("ave_ekin_"+phase+".dat");
  temp.open("ave_temp_"+phase+".dat");
  tot.open("ave_etot_"+phase+".dat");
  pres.open("ave_press_"+phase+".dat");

  ofstream Epot, Ekin, Etot, Temp, Press;

  Epot.open("argon_ave_epot_"+phase+".dat");
  Ekin.open("argon_ave_ekin_"+phase+".dat");
  Temp.open("argon_ave_temp_"+phase+".dat");
  Etot.open("argon_ave_etot_"+phase+".dat");
  Press.open("argon_ave_press_"+phase+".dat");

  for(int i=0;i<block_number;i++){
    pot >> m >> in >> s;
    Epot << m << ";" << in*epsilon << ";" << s*epsilon << endl;
    kin >> m >> in >> s;
    Ekin << m << ";" << in*epsilon << ";" << s*epsilon << endl;
    tot >> m >> in >> s;
    Etot << m << ";" << in*epsilon << ";" << s*epsilon << endl;
    temp >> m >> in >> s;
    Temp << m << ";" << in*t << ";" << s*t << endl;
    pres >> m >> in >> s;
    Press << m << ";" << in*p << ";" << s*p << endl;
  }
  pot.close();
  kin.close();
  temp.close();
  tot.close();
  pres.close();
  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
  Press.close();

  return;
}

void Krypton(){
  double sigma=0.364*1E-9;
  double epsilon=164.*k_be;
  double t = 164.;
  double p = t*k_bj*1E-5/pow(sigma,3.);
  double in;
  double s;
  int m;
  ifstream pot,kin,tot,pres,temp;
  pot.open("ave_epot_"+phase+".dat");
  kin.open("ave_ekin_"+phase+".dat");
  temp.open("ave_temp_"+phase+".dat");
  tot.open("ave_etot_"+phase+".dat");
  pres.open("ave_press_"+phase+".dat");

  ofstream Epot, Ekin, Etot, Temp, Press;

  Epot.open("krypton_ave_epot_"+phase+".dat");
  Ekin.open("krypton_ave_ekin_"+phase+".dat");
  Temp.open("krypton_ave_temp_"+phase+".dat");
  Etot.open("krypton_ave_etot_"+phase+".dat");
  Press.open("krypton_ave_press_"+phase+".dat");

  for(int i=0;i<block_number;i++){
    pot >> m >> in >> s;
    Epot << m << ";" << in*epsilon << ";" << s*epsilon << endl;
    kin >> m >> in >> s;
    Ekin << m << ";" << in*epsilon << ";" << s*epsilon << endl;
    tot >> m >> in >> s;
    Etot << m << ";" << in*epsilon << ";" << s*epsilon << endl;
    temp >> m >> in >> s;
    Temp << m << ";" << in*t << ";" << s*t << endl;
    pres >> m >> in >> s;
    Press << m << ";" << in*p << ";" << s*p << endl;
  }
  pot.close();
  kin.close();
  temp.close();
  tot.close();
  pres.close();
  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
  Press.close();

  return;
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
