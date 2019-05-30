/*************************************************************************
**************************************************************************
    _/_/_/_/       _/_/_/_/_/       Numerical Simulation Laboratory
   _/     _/      _/               Physics Department
  _/     _/      _/  _/_/_/       Università degli Studi di Milano
 _/     _/ _    _/      _/  _    Daniele Gilio
_/_/_/_/  /_/  _/_/_/_/_/  /_/  email: daniele.gilio@studenti.unimi.it
**************************************************************************
**************************************************************************/
#include "tsp.h"

int main(){
  Input();
  if(city_dist==false)
    cout << "Circumference City Distribution" << endl;
  else
    cout << "Square City Distribution" << endl;
  cout << "Acting as a God... (creating people, killing people, being merciful sometimes, you know... Go(o)d stuff)" << endl;
  for(int i=0;i<max_iterations;i++){
    temp=1./double(i+1.);
    generate(temp);

  }

  save();

  cout << "Done." << endl;

  return 0;
}
