
#ifndef _CPLEX_H
#define _CPLEX_H 1

#include "log.h"
#include "SQM_Instance.h"
#include <ilcplex/ilocplex.h>

extern double EPSILON;
extern double TIME_MAX;

ILOSTLBEGIN

typedef IloArray<IloBoolVarArray> BoolVarMatrix;
typedef IloArray<BoolVarMatrix> BoolVarArrayMatrix;
typedef IloArray<IloNumArray> NumMatrix;
typedef IloArray<IloIntArray> IntMatrix;

#endif

/* eof */
