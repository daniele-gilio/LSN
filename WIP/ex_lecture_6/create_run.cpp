#include <iostream>
#include <fstream>

using namespace std;

int main() {
  ofstream out;
  out.open("run.sh");
  double min = 0.5, max = 2.;
  double pass = (max-min)/10.;
  double count=min;
  //Set Gibbs Sampling and h=0
  out << "sed -i '4s/.*/0.0/' input.dat" << endl;
  out << "sed -i '5s/.*/0/' input.dat" << endl;

  for(int i=0;i<11;i++){
    out << "sed -i '1s/.*/"+to_string(count)+"/' input.dat" << endl << "./Monte_Carlo_ISING_1D.exe" << endl;
    count+=pass;
  }

  //Set h=0.02
  out << "sed -i '4s/.*/0.02/' input.dat" << endl;
  count = min;
  for(int i=0;i<11;i++){
    out << "sed -i '1s/.*/"+to_string(count)+"/' input.dat" << endl << "./Monte_Carlo_ISING_1D.exe" << endl;
    count+=pass;
  }

  //Set Metropolis Sampling and h=0
  out << "sed -i '5s/.*/1/' input.dat" << endl;
  out << "sed -i '4s/.*/0.0/' input.dat" << endl;
  count = min;
  
  for(int i=0;i<11;i++){
    out << "sed -i '1s/.*/"+to_string(count)+"/' input.dat" << endl << "./Monte_Carlo_ISING_1D.exe" << endl;
    count+=pass;
  }

  //Set h=0.02
  out << "sed -i '4s/.*/0.02/' input.dat" << endl;
  count = min;
  for(int i=0;i<11;i++){
    out << "sed -i '1s/.*/"+to_string(count)+"/' input.dat" << endl << "./Monte_Carlo_ISING_1D.exe" << endl;
    count+=pass;
  }

  out.close();
  return 0;
}
