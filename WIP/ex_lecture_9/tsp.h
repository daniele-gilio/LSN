#ifndef _TSP_
#define _TSP_
#include <cmath>
#include <iostream>
#include <fstream>
#include "random.h"

using namespace std;

int pop=900, n_matings=0, n_children=0, max_iterations=0, mutation_number=0, mutation_swaps=0, n_cities=0;
bool city_dist = false; //false for circumference, true for square
Random rnd;

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

city *cities = new city [30]; //also hard coded

class ind{
public:
  ind(){
    for(int i=0;i<n;i++){
      c[i]=0;
    }
  }

  //Copy constructor
  ind(const ind & s){
    for(int i=0;i<n;i++)
      c[i]=s.c[i];
  }

  ~ind(){
    delete []c;
  }

  void operator = (const ind & s){
    for(int i=0;i<n;i++)
      c[i]=s.c[i];
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

  void calc_cost(city * city){
    //cost=cost(cities, c, n);
    cost=0.;
    for(int i=0;i<n-1;i++){
      cost+=sqrt(pow(city[c[i]].getx()-city[c[i]+1].getx(),2.)+pow(city[c[i]].gety()-city[c[i]+1].gety(),2.));
    }
    cost+=sqrt(pow(city[c[n-1]].getx()-city[0].getx(),2.)+pow(city[c[n-1]].gety()-city[c[0]].gety(),2.));
  }

  double get_cost(){
    return cost;
  }

private:
  int n = 30; //the last one is the same as the first
  int *c = new int [n]; //number of cities needs to be hard coded here (this is a city index array)
  double cost; //embed cost into individual class

};

int Pbc(int index){
  if(index>29){
    return index-30;
  }

  else return index;
}

void Input(ind* population){
  //Read Parameters from file
  ifstream in;
  in.open("input.dat");
  in >> pop;
  in >> n_matings;
  in >> n_children;
  in >> max_iterations;
  in >> mutation_number;
  in >> mutation_swaps;
  in >> n_cities;
  in >> city_dist;
  in.close();

  //Inizialization of the Random Number Generator
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

  ofstream c;
  c.open("cities.dat");
  for(int i=0;i<n_cities;i++){
    c << cities[i].getx() << " " << cities[i].gety() << endl;
  }
  c.close();

  //Initialize Population
  //ind population[pop];

  //Population Creation via random generation
  /*for(int i=0;i<pop;i++){
    int *c = new int [n_cities];
    for(int j=0;j<n_cities;j++){
      c[j]=int(rnd.Rannyu(1.,30.)+0.5);
      for(int k=0;k<j;k++){
        if(c[j]==c[k]){
          c[j]=int(rnd.Rannyu(1.,30.)+0.5);
          k=0;
        }
      }
    }
    population[i].setc(c);
    population[i].calc_cost(cities);
  }*/

  //Initialization via array shuffling
  for(int i=0;i<pop;i++){
    int *c = new int [n_cities];
    for(int j=0;j<n_cities;j++){
      c[j]=j;
    }

    for(int j=0;j<100;j++){
      int a = rnd.Rannyu(0.,29.)+0.5;
      int b = rnd.Rannyu(0.,29.)+0.5;
      swap(c[a],c[b]);
    }

    population[i].setc(c);
    population[i].calc_cost(cities);
  }

  /*for(int i=0;i<pop;i++){
    cout << "Individual " << i+1 << endl;
    for(int j=0;j<n_cities;j++){
      cout << population[i].getc(j) << endl;
    }
    cout << "Cost = " << population[i].get_cost() << endl;
    cout << endl;
  }*/
}

void swap(ind* population, int i, int j){
  ind t = population[i];
  population[i]=population[j];
  population[j]=t;
}

void quickSort(ind *array, int low, int high)
{
    //cout << "k" << endl;
    int i = low;
    int j = high;
    int pivot = array[i+(i - j) / 2].get_cost();
    i=low+1;
    while (i <= j)
    {
        while (i<=j and array[i].get_cost() <= pivot)
            i++;
        while (i<=j and array[j].get_cost() > pivot)
            j--;
        if (i <= j)
        {
            swap(array,i,j);
            i++;
            j--;
        }
    }
    if (j > low)
        quickSort(array, low, j);
    if (i < high)
        quickSort(array, i, high);
}
void generate(int n, ind* population){
  //Sort for fitness
  /*for(int i=0;i<pop-1;i++)
    for(int j=i+1;j<pop;j++)
      if(population[j].get_cost()<population[i].get_cost())
        swap(population, i, j);*/

  //Mutations
  int p = rnd.Rannyu();
  if(p<0.05){
    for(int i=0;i<pop;i++){
      int *c = new int [n_cities];
      for(int j=0;j<n_cities;j++)
        c[j]=population[i].getc(j);
      int h = rnd.Rannyu(0.,28.)+0.5;
      swap(c[h],c[h+1]);
      population[i].setc(c);
      population[i].calc_cost(cities);
    }
  }

  //Sort for fitness
  /*for(int i=0;i<pop-1;i++)
    for(int j=i+1;j<pop;j++)
      if(population[j].get_cost()<population[i].get_cost())
        swap(population, i, j);*/
  quickSort(population,0,pop);


  ofstream out;
  out.open("best.dat",ios::app);
  out << n << " " << population[0].get_cost() << endl;
  //cout << n << " " << population[0].get_cost() << endl;
  out.close();
  /*for(int i=0;i<n_cities;i++){
    cout << population[0].getc(i) << endl;
  }
  cout << endl;*/


  //Mating Season


}

#endif
