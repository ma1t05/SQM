
#ifndef _SQM_HEURISTIC_H
#define _SQM_HEURISTIC_H

#include <utility>
#include <gmp.h>
#include "point.h"
#include "SQM.h"
#include "mp_jarvis.h"

typedef mpf_t num;
void SQM_heuristic(SQM_instance*,int,double,double,response_unit*);
response_unit* guess_a_location_01(int,int,point*); // Returns the first p
response_unit* guess_a_location_02(int,int,point*); // Returns p random with replace
response_unit* guess_a_location_03(int,int,point*); // Returns p random
double SQM_response_time
(SQM_instance*,
 int,/* Number of adjusters */
 response_unit* X,
 double lambda, /* mean arrival rate */
 double Mu_NT // mean of non-travel time component of the service time
 );

#endif

/* eof */
