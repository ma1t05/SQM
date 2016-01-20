
#ifndef _SQM_MODEL
#define _SQM_MODEL 1

#include <sstream>
#include "cplex.h"

void SQM_model
(SQM_instance&I, // Set of points
 int,            // facilities
 float         // speed
 );

int* SQM_model
(SQM_instance&, // Set of points
 int,           // facilities
 int,           // Number of facilities that care
 float          // speed
 );

#endif

/* eof */
