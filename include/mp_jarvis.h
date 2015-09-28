

#ifndef _MP_JARVIS_H
#define _MP_JARVIS_H 1

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <gmp.h>

#define epsilon 0.000001
using namespace std;

void jarvis_hypercube_approximation
(int C, /* Number of types of customers*/
 int N, /* Number of servers */
 mpf_t *lambda, /* arrive rate according to a Poisson process per type */
 mpf_t **Tao, /* expected service time for service i and customer of node m */
 int **a, /* for customers of type m, the list of preferred servers */
 mpf_t **f);

#endif

/* eof*/
