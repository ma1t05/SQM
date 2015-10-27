
#ifndef _SQM_H
#define _SQM_H

#include <math.h>
#include "point.h"
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

typedef IloArray<IloBoolVarArray> BoolVarMatrix;
typedef IloArray<BoolVarMatrix> BoolVarArrayMatrix;
typedef IloArray<IloNumArray> NumMatrix;
typedef IloArray<IloIntArray> IntMatrix;

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

#include "SQM_model.h"
#include "Goldberg.h"
#include "SQM_heuristic.h"
#include "config.h"
#include "gnuplot.h"
#include "SQM_GRASP.h"

#endif

/* eof */
