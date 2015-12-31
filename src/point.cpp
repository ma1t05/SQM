
#include "point.h"

float dist(point *i,point *j){
  return sqrt((i->x - j->x)*(i->x - j->x) + (i->y - j->y) * (i->y - j->y));
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

