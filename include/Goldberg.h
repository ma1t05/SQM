#ifndef _GOLDBERG_H_
#define _GOLDBERG_H_

#include "cplex_matrix.h"
#include "SQM_instance.h"

void Goldberg(instance*, int, float, float);
void gnuplot_goldberg(instance*,int,IloCplex*,IloBoolVarArray*,BoolVarArrayMatrix*);

#endif
