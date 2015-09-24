

#ifndef _MP_JARVIS_H
#define _MP_JARVIS_H 1

#include <gmp.h>
typedef mpf_t num;

num** jarvis_hypercube_approximation
(int C, /* Number of types of customers*/
 int N, /* Number of servers */
 num *lambda, /* */
 num **Tao, /* */
 int **a, /* */
 num **f /* */);

#endif

/* eof*/
