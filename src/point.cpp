
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
/*#include <utility>*/
#include "point.h"
#include "log.h"
using namespace std;

float dist(point *i,point *j){
  return sqrt((i->x - j->x)*(i->x - j->x) + (i->y - j->y) * (i->y - j->y));
}

instance* read_points(const char *filename){
  int i;
  instance *I;
  point *puntos;
  ifstream file(filename);
  
  if (!file) {
    logError(cerr << "ERROR: colud not open file '" << filename << "' for reading" << endl);
    return NULL;
  }

  I = new instance;
  file >> I->n >> I->S;
  puntos = new point[I->n];
  for (i = 0;i < I->n;i++)
    file >> puntos[i].x >> puntos[i].y >> puntos[i].demand;
  I->points = puntos;
  file.close();
  
  return I;
}

int unif(int a) {
  return floor(double(a) * rand() / RAND_MAX);
}

int comp(const void *a,const void *b) {
  std::pair<double,int> *x,*y;
  x = (std::pair<double,int>*)a;
  y = (std::pair<double,int>*)b;
  if (x->first > y->first) return 1;
  else if (x->first < y->first) return -1;
  else return 0;
}

void sort_dist (int n,double *d,int *c) {
  std::pair<double,int> *x;
  x = new std::pair<double,int>[n];
  for (int i = 0;i < n;i++) {
    x[i].first = d[i];
    x[i].second = i;
  }
  qsort(x,n,sizeof(std::pair<double,int>), comp);
  for (int i = 0;i < n;i++)
    c[i] = x[i].second;
  delete[] x;
}

