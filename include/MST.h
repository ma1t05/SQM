
#ifndef _MST_H
#define _MST_H 1

#include "SQM.h"
#include <gmp.h>

double MST_response_time
(SQM_instance *I,
 int p, /* Number of adjusters */
 response_unit* X, /* Current configuration */
 double lambda, // mean rate per unit of time within service calls are generated in Poisson manner
 double Mu_NT // mean of non-travel time component of the service time
 );

void MST_update_mst
(mpf_t *mst,       /* Stores the response time */
 int m,            /* number of clients */
 int p,            /* number of servers */
 double Mu_NT,     /* Mean of Non-Travel time */
 double **Dist,    /* Matrix of distances between clients & servers */
 response_unit* X, /* Current configuration */
 mpf_t **f         /* Curren matrix of assignations */
);

void MST_expected_travel_time
(mpf_t t_r, /* stores the response time */
 int m, /* number of clients */
 int p, /* number of servers */
 double **Dist,
 response_unit* X,
 mpf_t **f
 );

void MST_mean_queue_delay
(mpf_t t_r,
 int m,
 int p,
 mpf_t *Lambda,
 mpf_t *MST,
 mpf_t **Tao,
 mpf_t **f
 );

#endif

/* eof */
