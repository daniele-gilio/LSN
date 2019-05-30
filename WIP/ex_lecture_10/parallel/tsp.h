/*************************************************************************
**************************************************************************
    _/_/_/_/       _/_/_/_/_/       Numerical Simulation Laboratory
   _/     _/      _/               Physics Department
  _/     _/      _/  _/_/_/       Università degli Studi di Milano
 _/     _/ _    _/      _/  _    Daniele Gilio
_/_/_/_/  /_/  _/_/_/_/_/  /_/  email: daniele.gilio@studenti.unimi.it
**************************************************************************
**************************************************************************/
#ifndef _TSP_
#define _TSP_
#include <cmath>
#include <iostream>
#include <fstream>
#include "random.h"

using namespace std;

//Global Variables Needed for Classes to Work

int n_cities=0;

//------Classes---------

class city{
public:
  //Constructor
  city(){
    x=0.;
    y=0.;
  }
  //Functions
  void setcoord(double a, double b){
    x=a;
    y=b;
  }

  double getx(){
    return x;
  }

  double gety(){
    return y;
  }

private:
  double x,y;
};

class ind{
public:
  ind(){
    c = new int [n];
    for(int i=0;i<n;i++){
      c[i]=0;
    }
  }

  void setc(int s, int index){
    c[index]=s;
  }

  void setc(int *s){
    for(int i=0;i<n;i++){
      c[i]=s[i];
    }
  }

  int* getc(){
    return c;
  }

  int getc(int index){
    return c[index];
  }

  void calc_cost(city *city){
    //cost=cost(cities, c, n);
    cost=0.;
    for(int i=0;i<n-1;i++){
      cost+=distance(city[c[i]], city[c[i+1]]);
    }
    cost+=distance(city[c[n-1]], city[c[0]]);
  }

  double get_cost(){
    return cost;
  }

private:
  int n = n_cities;
  int *c = NULL;
  double cost=0.; //embed cost into individual class
  double distance(city a, city b){
    return sqrt(pow(a.getx()-b.getx(), 2.)+pow(a.gety()-b.gety(), 2.));
  }

};

//---------------------------Global Variables--------------------------------------------------------------

int pop=0, max_iterations=0, n_moves=0., best_rank = 0.;
bool city_dist = false; //false for circumference, true for square
Random rnd;
ind *population = NULL;
city *cities = NULL;
ofstream out, mean;
double temp = 1.;

//-------------------------Functions------------------------------------------------------------------------

int Pbc(int index){
  if(index>n_cities-1){
    return index-n_cities;
  }

  else return index;
}

void check(ind s){
  //Check if the sum of all the cities is correct
  double k=0.;
  for(int i=0;i<n_cities;i++)
    k+=s.getc(i);
  if(k==(n_cities-1)*n_cities/2){
    return;
  }

  else{
    std::cerr << "I fucked up :(" << '\n';
  }
}

void swap(ind* population, int i, int j){
  ind t=population[i];
  population[i]=population[j];
  population[j]=t;
}

bool find(int a, int *b, int lenght){
  for(int j=0;j<lenght;j++){
      if(a==b[j]){
        return true;
        break;
      }
    }
  return false;
}

double Boltzmann(double e_new, double e_old, double T){
  return exp(-(e_new-e_old)/T);
}

void creation(){

  int *c = new int [n_cities];
  for(int i=0;i<pop;i++){
    for(int j=0;j<n_cities;j++){
      c[j]=j;
    }

    for(int j=0;j<100;j++){
      int a = rnd.Rannyu(0.,n_cities-1.)+0.5;
      int b = rnd.Rannyu(0.,n_cities-1.)+0.5;
      swap(c[a],c[b]);
    }

    population[i].setc(c);
    check(population[i]);
    population[i].calc_cost(cities);
  }

  delete []c;

}

void Input(int rank){
  //Read Parameters from file
  ifstream in;
  in.open("input.dat");
  in >> pop;
  in >> max_iterations;
  in >> n_cities;
  in >> city_dist;
  in >> n_moves;
  in.close();

  //Initialize Population and Cities
  population = new ind [pop];
  cities = new city [n_cities];

  //Inizialization of the Random Number Generator (every node has the same generator -> everyone has the same set of cities)
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


  //Initialize cities coordinates

  if(city_dist==false){
    for(int i=0;i<n_cities;i++){
      double x=0., y=0., t=0.;
      x=rnd.Rannyu(-1.,1.);
      t=rnd.Rannyu();

      if(t<0.5){
        y=-sqrt(1.-x*x);
      }

      else{
        y=sqrt(1.-x*x);
      }

      cities[i].setcoord(x,y);
    }
  }

  else{
    for(int i=0;i<n_cities;i++){
      double x=0., y=0.;
      x=rnd.Rannyu();
      y=rnd.Rannyu();
      cities[i].setcoord(x,y);
    }

  }

  out.open("cities.dat");
  for(int i=0;i<n_cities;i++){
    out << cities[i].getx() << " " << cities[i].gety() << endl;
  }
  out.close();

  //Reinizialization of the Random Number Generator (everyone has a different generator)
  Primes.open("Primes");
  if (Primes.is_open()){
    //Read different primes for each node
    for(int i=0;i<rank+1;i++){
          Primes >> p1 >> p2 ;
    }
  } else cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();
  input.open("seed.in");
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

  creation();
}



void generate(double T){

  ind test;

  //Mutations
  double p_mut = 0.;
  for(int i=0;i<pop;i++){

    for(int k=0;k<n_moves;k++){
      test.setc(population[i].getc());
      p_mut = rnd.Rannyu();
      //Pair Swap
      if(p_mut<0.4){
        int *c = new int [n_cities];
        for(int j=0;j<n_cities;j++)
          c[j]=test.getc(j);
        int h = rnd.Rannyu(0.,n_cities-2.)+0.5;
        swap(c[h],c[h+1]);
        test.setc(c);
        check(test);
        delete []c;
      }

      //Contiguous Cities Translation
      if(p_mut>0.4 and p_mut<0.8){
        int *c = new int [n_cities];
        for(int j=0;j<n_cities;j++)
          c[j]=test.getc(j);

        int h = rnd.Rannyu(0.,n_cities/2.-0.5)+0.5; //Start
        int l = rnd.Rannyu(n_cities/2.+0.5, n_cities-2.)+0.5; //Finish
        int g = rnd.Rannyu(0.,n_cities/2.)+0.5; //Lenght of Translation

        for(int j=h;j<l;j++)
          swap(c[Pbc(j+g)],c[j]);
        test.setc(c);
        check(test);
        delete []c;
      }

      //Contiguous city inversion
      if(p_mut>0.8 and p_mut<1.){
        int *c = new int [n_cities];
        for(int j=0;j<n_cities;j++)
          c[j]=test.getc(j);

        int h = rnd.Rannyu(0.,n_cities/2.-0.5)+0.5; //Start
        int l = rnd.Rannyu(n_cities/2.+0.5, n_cities-2.)+0.5; //Finish
        int g = l;

        for(int j=h;j<int(l/2.+0.5);j++){
          swap(c[Pbc(j)],c[Pbc(g)]);
          g--;
        }
        test.setc(c);
        check(test);
        delete []c;
      }

    }

    test.calc_cost(cities);

    if(test.get_cost()<population[i].get_cost()){
      population[i]=test;
      check(population[i]);
    }

    else if(rnd.Rannyu()<=Boltzmann(test.get_cost(),population[i].get_cost(),T)){
      population[i]=test;
      check(population[i]);
    }
  }

}

void save(){
  //Save best path
  if(city_dist==false)
    out.open("best_path_circ.dat");
  else
    out.open("best_path_sq.dat");

  for(int i=0;i<n_cities;i++){
    out << cities[population[0].getc(i)].getx() << " " << cities[population[0].getc(i)].gety() << endl;
  }

  out << cities[population[0].getc(0)].getx() << " " << cities[population[0].getc(0)].gety() << endl;


  out.close();
}

void master_save_cost(double *best, int size, double t){
  if(city_dist==0)
    out.open("best_cost_circ.dat", ios::app);
  else
    out.open("best_cost_sq.dat", ios::app);
  double min = best[0];
  best_rank = 0;
  for(int i=1;i<size-1;i++)
    if(best[i]<min){
      min=best[i];
      best_rank=i;
    }

  out << t << " " << min << endl;
  out.close();
}

#endif
