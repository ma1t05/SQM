
#ifndef _POINT_H_
#define _POINT_H_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
/*#include <utility>*/
using namespace std;

struct point {
  float x,y;
  /* (x,y) coordinates of the point */
  float demand;
  /* demand of the point */
  bool *N_i;
  /* array that indicates whether the point j is in range S */
};

typedef struct point point;

float dist(point*,point*);
int comp(const void*,const void*);
void sort_dist (int,double*,int*);

#endif

/* eof */
