
#ifndef _SQM_GRASP
#define _SQM_GRASP 1

#include "SQM.h"

response_unit* GRASP
(SQM_instance *I, // Instance
 int p, // Number of adjusters
 double lambda, // mean rate per unit of time within service calls are generated in Poisson manner
 double Mu_NT, // mean of non-travel time component of the service time
 double v, // Speed
 double alpha
 );

double GRASP_func_NN
(SQM_instance *I, // Instance
 int p, // Number of adjusters
 response_unit *X, // Current solution
 double lambda, // mean rate per unit of time within service calls are generated in Poisson manner
 double Mu_NT // mean of non-travel time component of the service time
 );

double GRASP_func_kNN
(SQM_instance *I, // Instance
 int p, // Number of adjusters
 response_unit *X, // Current solution
 double lambda, // mean rate per unit of time within service calls are generated in Poisson manner
 double Mu_NT, // mean of non-travel time component of the service time
 int K // Number of servers to consider
 );

#endif

/* eof */
