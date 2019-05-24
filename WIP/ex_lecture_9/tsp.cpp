#include "tsp.h"

int main(){
  ind *population = new ind [900];
  Input(population);
  cout << "Acting as a god ..." << endl;
  for(int i=0;i<max_iterations;i++){
    generate(i, population);
    if(i%10==0)
      cout << double(i)*100./double(max_iterations) << " %" << endl;
  }

  return 0;
}
