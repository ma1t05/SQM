
#ifndef _SQM_MODEL
#define _SQM_MODEL

#include "point.h"
#include <ilcplex/ilocplex.h>

ILOSTLBEGIN

typedef IloArray<IloBoolVarArray> BoolVarMatrix;
typedef IloArray<BoolVarMatrix> BoolVarArrayMatrix;
typedef IloArray<IloNumArray> NumMatrix;
typedef IloArray<IloIntArray> IntMatrix;
void SQM_model
(instance* I, // Set of points
 int p, // facilities
 float mu, // rate parameter
 float f, // portion of demand
 float speed); // speed

#endif

/* eof */
