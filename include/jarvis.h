
#ifndef _JARVIS_H
#define _JARVIS_H 1

#include <iostream>
#include <cmath>
#define epsilon 0.0000000001
using namespace std;

double** jarvis_hypercube_approximation
(int C, /* Number of types of customers*/
 int N, /* Number of servers */
 double *lambda, /* */
 double **Tao, /* */
 int **a /* */);

#endif

  /* eof */
