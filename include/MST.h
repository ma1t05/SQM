
#ifndef _MST_H
#define _MST_H 1

#include "SQM_Solution.h"
#include <gmp.h>

double MST_response_time (SQM_solution* Sol /* Current configuration */);

void MST_update_mst
(mpf_t *mst,       /* Stores the response time */
 SQM_solution *Sol,/* Current configuration */
 mpf_t **f         /* Curren matrix of assignations */
);

void MST_expected_travel_time
(mpf_t t_r,        /* stores the response time */
 SQM_solution *Sol,/* Current configuration */
 mpf_t **f         /* Curren matrix of assignations */
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
