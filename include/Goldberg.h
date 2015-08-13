#ifndef _GOLDBERG_H_
#define _GOLDBERG_H_

#include "SQM.h"
#include <sstream>

void Goldberg(SQM_instance*, int, float, float);
void gnuplot_goldberg(SQM_instance*,int,IloCplex*,IloBoolVarArray*,BoolVarArrayMatrix*);

#endif
