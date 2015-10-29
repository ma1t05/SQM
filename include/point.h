
#ifndef _POINT_H_
#define _POINT_H_

struct point {
  float x,y;
  /* (x,y) coordinates of the point */
  float demand;
  /* demand of the point */
  bool *N_i;
  /* array that indicates whether the point j is in range S */
};

struct instance {
  int n;
  /* numbers of nodes */
  float S;
  /* maximum distance */
  struct point *points;
  /* set of nodes */
};

typedef struct point point;
typedef struct instance instance;

float dist(point*,point*);
instance* read_points(const char*);
int unif(int);
int comp(const void*,const void*);
void sort_dist (int,double*,int*);

#endif

/* eof */
