
#ifndef _JARVIS_H
#define _JARVIS_H 1

#include <iostream>
#include <cmath>
#include <gmp.h>
#define epsilon 0.00000001
using namespace std;
typedef long double num;

num** jarvis_hypercube_approximation
(int C, /* Number of types of customers*/
 int N, /* Number of servers */
 num *lambda, /* */
 num **Tao, /* */
 int **a, /* */
 num **f /* */);

#endif

  /* eof */
