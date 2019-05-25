#ifndef _TSP_
#define _TSP_
#include <cmath>
#include <iostream>
#include <fstream>
#include "random.h"

using namespace std;

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
  int n = n_cities; //the last one is the same as the first
  int *c = NULL;
  double cost=0.; //embed cost into individual class
  double distance(city a, city b){
    return sqrt(pow(a.getx()-b.getx(), 2.)+pow(a.gety()-b.gety(), 2.));
  }

};

//---------------------------Global Variables--------------------------------------------------------------

int pop=0, n_children=0, max_iterations=0, mercy_number=0;
bool city_dist = false; //false for circumference, true for square
Random rnd;
ind *population = NULL;
city *cities = NULL;
ofstream out;

//-------------------------Functions------------------------------------------------------------------------

int Pbc(int index){
  if(index>n_cities-1){
    return index-n_cities;
  }

  else return index;
}

void check(ind s){
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

  //Sort for fitness
  for(int i=0;i<pop-1;i++)
    for(int j=i+1;j<pop;j++)
      if(population[j].get_cost()<population[i].get_cost())
        swap(population, i, j);

}

void Input(){
  //Read Parameters from file
  ifstream in;
  in.open("input.dat");
  in >> pop;
  in >> n_children;
  in >> max_iterations;
  in >> mercy_number;
  in >> n_cities;
  in >> city_dist;
  in.close();

  //Initialize Population and Cities
  population = new ind [pop];
  cities = new city [n_cities];

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

  out.open("cities.dat");
  for(int i=0;i<n_cities;i++){
    out << cities[i].getx() << " " << cities[i].gety() << endl;
  }
  out.close();

  creation();

  /*for(int i=0;i<pop;i++){
    cout << "Individual " << i+1 << endl;
    for(int j=0;j<n_cities;j++){
      cout << population[i].getc(j) << endl;
    }
    cout << "Cost = " << population[i].get_cost() << endl;
    cout << endl;
  }*/
}



void generate(int n){
  if(city_dist==false)
    out.open("best_circ.dat",ios::app);
  else
    out.open("best_sq.dat",ios::app);
  out << n << " " << population[0].get_cost() << endl;
  //cout << n << " " << population[0].get_cost() << endl;
  out.close();

  //Selection
  int *sel_ind = NULL;
  ind *lucky_boys = NULL;
  if(mercy_number!=0){
    sel_ind = new int [mercy_number];
    lucky_boys = new ind [mercy_number];
    for(int i=0;i<mercy_number;i++){
      sel_ind[i] = pop*pow(rnd.Rannyu(),10.);
      lucky_boys[i] = population[sel_ind[i]];
    }
  }

  //Mating Season
  for(int k=pop/2;k<pop;k++){
    ind mother, father;
    int cut_1=0;
    int sel_index=0;
    sel_index = pop*pow(rnd.Rannyu(),30.);
    mother = population[sel_index];
    check(mother);
    sel_index = pop*pow(rnd.Rannyu(),30.);
    father = population[sel_index];
    check(father);

    double p_mate = rnd.Rannyu();
    if(p_mate<0.9){
      ind *children = new ind [n_children];

      for(int i=0;i<n_children;i++){
        cut_1=rnd.Rannyu(0.,n_cities/2.)+0.5;
        //cut_1=15;

        for(int j=0;j<cut_1;j++)
          children[i].setc(mother.getc(j), j);

        for(int j=cut_1;j<n_cities;j++)
          for(int k=0;k<n_cities;k++)
            if(find(father.getc(k),children[i].getc(),n_cities)==false){
              children[i].setc(father.getc(k),j);
              break;
            }
        int index = rnd.Rannyu(n_children, n_cities-n_children)+0.5;
        check(children[i]);
        population[k-i-1]=children[i];
        check(population[k-i-1]);
        population[k-i-1].calc_cost(cities);

      }
      delete []children;
    }
  }

  //Mutations
  double p_mut = 0.;
  for(int i=0;i<pop;i++){
      p_mut = rnd.Rannyu();
      //Pair Swap
      if(p_mut<0.05){
        int *c = new int [n_cities];
        for(int j=0;j<n_cities;j++)
          c[j]=population[i].getc(j);
        int h = rnd.Rannyu(0.,n_cities-2.)+0.5;
        swap(c[h],c[h+1]);
        population[i].setc(c);
        check(population[i]);
        //population[i].calc_cost(cities);
        delete []c;
      }

      //Whole Translation
      if(n%30==0){
        int *c = new int [n_cities];
        for(int j=0;j<n_cities;j++)
          c[j]=population[i].getc(j);
        int h = rnd.Rannyu(0.,n_cities/2.)+0.5; //Max Half Positions Slide
        for(int j=0;j<n_cities;j++){
          swap(c[Pbc(j+h)],c[j]);
        }

        population[i].setc(c);
        check(population[i]);
        delete []c;
      }

      //Contiguous Cities Translation
      if(p_mut>0.05 and p_mut<0.1){
        int *c = new int [n_cities];
        for(int j=0;j<n_cities;j++)
          c[j]=population[i].getc(j);

        int h = rnd.Rannyu(0.,n_cities/2.-0.5)+0.5; //Start
        int l = rnd.Rannyu(n_cities/2.+0.5, n_cities-2.)+0.5; //Finish
        int g = rnd.Rannyu(0.,n_cities/2.)+0.5; //Lenght of Translation

        for(int j=h;j<l;j++)
          swap(c[Pbc(j+g)],c[j]);
        population[i].setc(c);
        check(population[i]);
        delete []c;
      }

      //Delete and reconstruct
      if(p_mut>0.1 and p_mut<0.11){
        int *c = new int [n_cities];
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
        delete []c;
      }

      population[i].calc_cost(cities);
  }

  for(int i=0;i<mercy_number;i++){
    int t = rnd.Rannyu(0.,pop-1.)+0.5;
    population[t]=lucky_boys[i];
    check(population[t]);
  }

  delete []lucky_boys;
  delete []sel_ind;
  //Sort for fitness
  for(int i=0;i<pop-1;i++){
    for(int j=i+1;j<pop;j++)
      if(population[j].get_cost()<population[i].get_cost())
        swap(population, i, j);
    check(population[i]);
  }


}

void save(){
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

#endif
