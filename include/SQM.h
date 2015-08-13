
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

#define EPSILON 0.001

struct network {
  int n;
};

struct SQM_instance {
  int M; /* Number of customers */
  int N; /* Number of potencial sites to locate a server */
  point* V; /* Set of demand points */
  point* W; /* Set of potencial locations sites */
};

typedef struct network network;
typedef struct SQM_instance SQM_instance;

#endif


/* eof */

