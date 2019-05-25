#include "tsp.h"

int main(){
  //ind *population = new ind [900];
  Input();
  if(city_dist==false)
    cout << "Circumference City Distribution" << endl;
  else
    cout << "Square City Distribution" << endl;
  cout << "Acting as a God ..." << endl;
  for(int i=0;i<max_iterations;i++){
    generate(i);
    /*if(i%10==0)
      cout << double(i)*100./double(max_iterations) << " %" << endl;*/
  }

  save();

  cout << "Done." << endl;

  return 0;
}
