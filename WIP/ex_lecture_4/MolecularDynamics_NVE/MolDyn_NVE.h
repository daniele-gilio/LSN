/*************************************************************************
**************************************************************************
    _/_/_/_/       _/_/_/_/_/       Numerical Simulation Laboratory
   _/     _/      _/               Physics Department
  _/     _/      _/  _/_/_/       Università degli Studi di Milano
 _/     _/ _    _/      _/  _    Daniele Gilio
_/_/_/_/  /_/  _/_/_/_/_/  /_/  email: daniele.gilio@studenti.unimi.it
**************************************************************************
**************************************************************************/

#include <string> //---------------------------------
//Needed here to use string in Execute function
using namespace std;//-------------------------------

//parameters, observables
const int m_props=4;
unsigned int block_number = 100;
int n_props;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_p;
double *ave_pot = new double [block_number];
double *prog_pot = new double [block_number];
double *ave_kin = new double [block_number];
double *prog_kin = new double [block_number];
double *ave_etot = new double [block_number];
double *prog_etot = new double [block_number];
double *ave_temp = new double [block_number];
double *prog_temp = new double [block_number];
double *ave_p = new double [block_number];
double *prog_p = new double [block_number];
double *ave_pot_2 = new double [block_number];
double *prog_pot_2 = new double [block_number];
double *ave_kin_2 = new double [block_number];
double *prog_kin_2 = new double [block_number];
double *ave_etot_2 = new double [block_number];
double *prog_etot_2 = new double [block_number];
double *ave_temp_2 = new double [block_number];
double *prog_temp_2 = new double [block_number];
double *ave_p_2 = new double [block_number];
double *prog_p_2 = new double [block_number];
double *sigma_pot = new double [block_number];
double *sigma_kin = new double [block_number];
double *sigma_etot = new double [block_number];
double *sigma_temp = new double [block_number];
double *sigma_p = new double [block_number];

// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part],x_t[m_part],y_t[m_part],z_t[m_part],xold_t[m_part],yold_t[m_part],zold_t[m_part];
double vx[m_part],vy[m_part],vz[m_part],vx_t[m_part],vy_t[m_part],vz_t[m_part];

// thermodynamical state
unsigned int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed;
double delta;

//functions
void Input(bool restart, unsigned int counter, unsigned int phase);
void Move(void);
void ConfFinal(void);
void ConfOld(void);
void ConfXYZ(int);
void Measure(unsigned int i);
double Force(int, int);
double Pbc(double);
double error(double* ave, double* av2, int n);
void Execute(string s, unsigned int phase);

//exceptions
bool ex = false;

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
