
#ifndef _SQM_MODEL
#define _SQM_MODEL 1

#include "point.h"
#include "instance-creator.h"
#include <sstream>

void SQM_model
(instance* I, // Set of points
 int p, // facilities
 float mu, // rate parameter
 float f, // portion of demand
 float speed); // speed
void SQM_model
(SQM_instance* I, // Set of points
 int p, // facilities
 int k, // Number of facilities that care
 float mu, // rate parameter
 float f, // portion of demand
 float speed); // speed

#endif

/* eof */
