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
#include "mpi.h"

int main(int argc, char* argv[]){
  MPI::Init(argc, argv);
  int size = MPI::COMM_WORLD.Get_size();
  int rank = MPI::COMM_WORLD.Get_rank();

  //Initialize
  Input(rank);

  double *best = new double [size];
  if(rank==0){
    if(city_dist==false)
      cout << "Circumference City Distribution" << endl;
    else
      cout << "Square City Distribution" << endl;
  }

  for(int i=0;i<max_iterations;i++){
    temp=1./double(i+1.);

    if(rank>0)
      generate(temp);

    double cost = population[0].get_cost();
    MPI_Gather(&cost,1,MPI_DOUBLE_PRECISION,best,1,MPI_DOUBLE_PRECISION,0,MPI::COMM_WORLD);

    if(rank==0)
      master_save_cost(best,size,temp);

  }

  MPI_Bcast(&best_rank,1,MPI_INTEGER,0,MPI::COMM_WORLD);

  if(rank==best_rank)
    save();

  if(rank==0)
    cout << "Done." << endl;

  MPI::Finalize();

  return 0;
}
