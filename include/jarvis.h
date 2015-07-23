
#ifndef _JARVIS_H
#define _JARVIS_H 1

#include <cmath>

double jarvis_hypercube_approximation
(int C, /* Number of types of customers*/
 int N, /* Number of servers */
 double *lambda, /* */
 double **Tao, /* */
 int **a /* */);

#endif

  /* eof */
