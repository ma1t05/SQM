
#ifndef _SQM_MODEL
#define _SQM_MODEL 1

#include <sstream>
#include "cplex.h"

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
