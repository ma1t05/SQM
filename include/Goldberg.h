
#ifndef _GOLDBERG_H_
#define _GOLDBERG_H_

#include "SQM_model.h"

void Goldberg(SQM_instance*, int, float, float);
void gnuplot_goldberg(SQM_instance*,int,IloCplex*,IloBoolVarArray*,BoolVarArrayMatrix*);

#endif
