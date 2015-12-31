
#ifndef _CPLEX_H
#define _CPLEX_H 1

#include "log.h"
#include "SQM_Instance.h"
#include <ilcplex/ilocplex.h>

#define EPSILON 0.00001
#define TIME_MAX 1200.0

ILOSTLBEGIN

typedef IloArray<IloBoolVarArray> BoolVarMatrix;
typedef IloArray<BoolVarMatrix> BoolVarArrayMatrix;
typedef IloArray<IloNumArray> NumMatrix;
typedef IloArray<IloIntArray> IntMatrix;

#endif

/* eof */
