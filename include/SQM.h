
#ifndef _SQM_H
#define _SQM_H

#include <math.h>
using namespace std;

#define EPSILON 0.001

struct response_unit {
  int location;
  double v;
  double beta;
};

struct network {
  int n;
};

typedef struct response_unit response_unit;
typedef struct network network;

#endif

/* eof */
