
#ifndef _SQM_H
#define _SQM_H

#include <fstream>
#include <math.h>
#include "point.h"

using namespace std;

#define            EPSILON 0.001
#define           TIME_MAX 1200.0
#define     MINS_PER_BLOCK 60
#define BLOCKS_PER_HORIZON 24

struct network {
  int n;
};
typedef struct network network;

struct SQM_instance {
  int M; /* Number of customers */
  int N; /* Number of potencial sites to locate a server */
  point* V; /* Set of demand points */
  point* W; /* Set of potencial locations sites */
};
typedef struct SQM_instance SQM_instance;

struct response_unit {
  int location;
  int past_location;
  double v;
  double beta;
};
typedef struct response_unit response_unit;

extern std::ofstream LogFile;
extern std::ofstream results;
extern std::ofstream dat;

#endif

/* eof */

