
#ifndef _SQM_HEURISTIC_H
#define _SQM_HEURISTIC_H

#include "SQM_Solution.h"

void SQM_heuristic(SQM_solution*,double,double);
double SQM_response_time
(SQM_solution*,
 double lambda, /* mean arrival rate */
 double Mu_NT // mean of non-travel time component of the service time
 );

#endif

/* eof */
