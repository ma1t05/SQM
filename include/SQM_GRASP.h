
#ifndef _SQM_GRASP
#define _SQM_GRASP 1

#include "SQM_Solution.h"

SQM_solution* GRASP
(SQM_instance *I, // Instance
 int p, // Number of adjusters
 double lambda, // mean rate per unit of time within service calls are generated in Poisson manner
 double Mu_NT, // mean of non-travel time component of the service time
 double v, // Speed
 double beta,
 double alpha
 );

double GRASP_func_NN (SQM_solution *Sol /* Current solution*/);
double GRASP_func_kNN (SQM_solution *Sol /* Current solution */, int K /* Number of servers to consider */);
extern int GRASP_kNN_param;

#endif

/* eof */
