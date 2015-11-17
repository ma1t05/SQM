
#ifndef _SQM_HEURISTIC_H
#define _SQM_HEURISTIC_H

#include <gmp.h>
#include "SQM.h"

typedef mpf_t num;
void SQM_heuristic(SQM_instance*,int,double,double,response_unit*);
double SQM_response_time
(SQM_solution*,
 double lambda, /* mean arrival rate */
 double Mu_NT // mean of non-travel time component of the service time
 );

#endif

/* eof */
