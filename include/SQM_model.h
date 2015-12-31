
#ifndef _SQM_MODEL
#define _SQM_MODEL 1

#define EPSILON 0.00001
#define TIME_MAX 1200.0

#include <sstream>
#include "SQM_Instance.h"
#include "log.h"
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

typedef IloArray<IloBoolVarArray> BoolVarMatrix;
typedef IloArray<BoolVarMatrix> BoolVarArrayMatrix;
typedef IloArray<IloNumArray> NumMatrix;
typedef IloArray<IloIntArray> IntMatrix;


void SQM_model
(SQM_instance* I, // Set of points
 int p, // facilities
 float mu, // rate parameter
 float f, // portion of demand
 float speed); // speed

int* SQM_model
(SQM_instance* I, // Set of points
 int p, // facilities
 int k, // Number of facilities that care
 float mu, // rate parameter
 float f, // portion of demand
 float speed); // speed

#endif

/* eof */
