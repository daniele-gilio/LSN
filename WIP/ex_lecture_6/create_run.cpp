#include <iostream>
#include <fstream>

using namespace std;

int main() {
  ofstream out;
  out.open("run.sh");
  double min = 0.5, max = 2.;
  double pass = (max-min)/10.;
  double count=min;
  //Set Gibbs Sampling, h=0 and no restart
  out << "sed -i '4s/.*/0.0/' input.dat" << endl;
  out << "sed -i '5s/.*/0/' input.dat" << endl;
  out << "sed -i '8s/.*/0/' input.dat" << endl;
  for(int i=0;i<11;i++){
    //Set Input and do a run for equilibration
    out << "sed -i '1s/.*/"+to_string(count)+"/' input.dat" << endl << "./Monte_Carlo_ISING_1D.exe" << endl;
    //Set restart to 1
    //out << "sed -i '8s/.*/1/' input.dat" << endl;
    //Rerun to take measures
    //out << "./Monte_Carlo_ISING_1D.exe" << endl;
    count+=pass;
  }

  //Set h=0.02 and no restart
  out << "sed -i '4s/.*/0.02/' input.dat" << endl;
  //out << "sed -i '8s/.*/0/' input.dat" << endl;
  count = min;
  for(int i=0;i<11;i++){
    out << "sed -i '1s/.*/"+to_string(count)+"/' input.dat" << endl << "./Monte_Carlo_ISING_1D.exe" << endl;
    //Set restart to 1
    //out << "sed -i '8s/.*/1/' input.dat" << endl;
    //Rerun to take measures
    //out << "./Monte_Carlo_ISING_1D.exe" << endl;
    count+=pass;
  }

  //Set Metropolis Sampling, h=0 and no restart
  out << "sed -i '5s/.*/1/' input.dat" << endl;
  out << "sed -i '4s/.*/0.0/' input.dat" << endl;
  //out << "sed -i '8s/.*/0/' input.dat" << endl;
  count = min;

  for(int i=0;i<11;i++){
    out << "sed -i '1s/.*/"+to_string(count)+"/' input.dat" << endl << "./Monte_Carlo_ISING_1D.exe" << endl;
    //Set restart to 1
    //out << "sed -i '8s/.*/1/' input.dat" << endl;
    //Rerun to take measures
    //out << "./Monte_Carlo_ISING_1D.exe" << endl;
    count+=pass;
  }

  //Set h=0.02 and no restart
  out << "sed -i '4s/.*/0.02/' input.dat" << endl;
  //out << "sed -i '8s/.*/0/' input.dat" << endl;
  count = min;
  for(int i=0;i<11;i++){
    out << "sed -i '1s/.*/"+to_string(count)+"/' input.dat" << endl << "./Monte_Carlo_ISING_1D.exe" << endl;
    //Set restart to 1
    //out << "sed -i '8s/.*/1/' input.dat" << endl;
    //Rerun to take measures
    //out << "./Monte_Carlo_ISING_1D.exe" << endl;
    count+=pass;
  }

  out.close();
  return 0;
}
