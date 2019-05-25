#include "tsp.h"

int main(){
  //ind *population = new ind [900];
  Input();
  if(city_dist==false)
    cout << "Circumference City Distribution" << endl;
  else
    cout << "Square City Distribution" << endl;
  cout << "Acting as a God... (creating people, killing people, being merciful sometimes, you know... Go(o)d stuff)" << endl;
  for(int i=0;i<max_iterations;i++){
    generate(i);
  }

  save();

  cout << "Done." << endl;

  return 0;
}
