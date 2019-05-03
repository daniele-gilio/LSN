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
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{
  Input(); //Inizialization
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();

//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;

  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> restart;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility

  n_props = 4; //Number of observables

//initial configuration
  if(restart==false){
    for (int i=0; i<nspin; ++i)
    {
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;
    }

  //Evaluate energy etc. of the initial configuration
    Measure();

  //Print initial values for the potential energy and virial
    cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
  }

  else{
    ReadInput.open("config.final");
    for (int i=0; i<nspin; ++i)
    {
      ReadInput >> s[i];
    }
    ReadInput.close();
  //Evaluate energy etc. of the initial configuration
    Measure();

  //Print initial values for the potential energy and virial
    cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
  }
}


void Move(int metro)
{
  int o;
  double /*p,*/ energy_old/*, energy_new, sm*/;
  //double energy_up, energy_down;
  double delta, r, A;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);
    energy_old=Boltzmann(s[o],o);
    delta=Boltzmann(-s[o],o)-energy_old;

    if(metro==1){ //Metropolis
      if(delta<0){
        s[o]=-s[o];
        accepted++;
      }
      else{
        r=rnd.Rannyu();
        A=min(1., exp(-beta*delta));
        if(r<A){
          s[o]=-s[o];
          accepted++;
        }
      }
      attempted++;
    }
    else{ //Gibbs sampling
        A = pow(1.+exp(-beta*delta),-1.);
        r = rnd.Rannyu();
        if(r<A){
          accepted++;
        }
        else{
          s[o]=-s[o];
          accepted++;
        }
      attempted++;
    }
  }
}

double Boltzmann(int sm, int ip)
{
    double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
    return ene;
}

void Measure()
{
  //int bin;
  double u = 0.0/*, m = 0.0*/, u_2=0., s_i=0.;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     u_2 += pow(-J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]), 2.);
     s_i += s[i];
  }
  walker[iu] = u;
  walker[ic] = u_2;
  walker[ix] = pow(s_i, 2.);
  walker[im] = s_i;
}


void Reset(int iblk) //Reset block averages
{

   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{

   ofstream Ene, Heat, Mag, Chi;
   const int wd=12;

    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    if(h==0){
      Ene.open("output.ene.0",ios::app);
      stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
      glob_av[iu]  += stima_u;
      glob_av2[iu] += stima_u*stima_u;
      err_u=Error(glob_av[iu],glob_av2[iu],iblk);
      Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
      Ene.close();

      Heat.open("output.heat.0",ios::app);
      stima_c = ((blk_av[ic]/blk_norm/(double)nspin)-pow(stima_u, 2.))*pow(beta, 2.); //Heat Capacity
      glob_av[ic]  += stima_c;
      glob_av2[ic] += stima_c*stima_c;
      err_c=Error(glob_av[ic],glob_av2[ic],iblk);
      Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
      Heat.close();

      Chi.open("output.x.0",ios::app);
      stima_x = beta*(blk_av[ix]/blk_norm/(double)nspin); //Chi
      glob_av[ix]  += stima_x;
      glob_av2[ix] += stima_x*stima_x;
      err_x=Error(glob_av[ix],glob_av2[ix],iblk);
      Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
      Chi.close();

      if(iblk==nblk){
        if(metro==1){
          Ene.open("output.ene_m.t",ios::app);
          Ene << temp << ";" << stima_u << ";" << err_u << endl;
          Ene.close();

          Heat.open("output.heat_m.t",ios::app);
          Heat << temp << ";" << stima_c << ";" << err_c << endl;
          Heat.close();

          Chi.open("output.x_m.t",ios::app);
          Chi << temp << ";" << stima_x << ";" << err_x << endl;
          Chi.close();
        }

        else{
          Ene.open("output.ene_g.t",ios::app);
          Ene << temp << ";" << stima_u << ";" << err_u << endl;
          Ene.close();

          Heat.open("output.heat_g.t",ios::app);
          Heat << temp << ";" << stima_c << ";" << err_c << endl;
          Heat.close();

          Chi.open("output.x_g.t",ios::app);
          Chi << temp << ";" << stima_x << ";" << err_x << endl;
          Chi.close();
        }

      }
    }

    else{
      Mag.open("output.mag.0",ios::app);
      stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetization
      glob_av[im]  += stima_m;
      glob_av2[im] += stima_m*stima_m;
      err_m=Error(glob_av[im],glob_av2[im],iblk);
      Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd+2) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
      Mag.close();

      if(iblk==nblk){
        if(metro==1){
          Mag.open("output.mag_m.t",ios::app);
          Mag << temp << ";" << stima_m << ";" << err_m << endl;
          Mag.close();
        }

        else{
          Mag.open("output.mag_g.t",ios::app);
          Mag << temp << ";" << stima_m << ";" << err_m << endl;
          Mag.close();
        }
      }
    }



    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
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
