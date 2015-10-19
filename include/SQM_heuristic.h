
#ifndef _SQM_HEURISTIC_H
#define _SQM_HEURISTIC_H

#include <utility>
#include <cstdlib>
#include <gmp.h>
#include "point.h"
#include "SQM.h"
#include "mp_jarvis.h"

struct response_unit {
  int location;
  int past_location;
  double v;
  double beta;
};

typedef mpf_t num;
typedef struct response_unit response_unit;
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
response_unit* SQM_GRASP
(SQM_instance *I,
 int p, // Number of adjusters
 double lambda, // mean rate per unit of time within service calls are generated in Poisson manner
 double Mu_NT, // mean of non-travel time component of the service time
 double v // Speed
 );

#endif

/* eof */
