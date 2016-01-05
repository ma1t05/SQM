
#ifndef _MST_H
#define _MST_H 1

#include "SQM_Server.h"
#include "mp_jarvis.h"
#include "log.h"

extern int MINS_PER_BLOCK;
extern int BLOCKS_PER_HORIZON;
/* */

double MST_response_time 
(SQM_instance* Inst, /* Current configuration */
 int p,
 server *Servers,
 int **preferred_servers
);

double* MST_workload
(SQM_instance* Inst, /* Current configuration */
 int p,
 server *Servers,
 int **preferred_servers
);

void MST_Calibration
(mpf_t **f,      /* Stores the assignations */
 mpf_t *mst,     /* Stores the response time */
 mpf_t **Tao,    /* Stores the service times */
 mpf_t *Lambda,  /* Stores the arrival rate */
 SQM_instance *Inst, /* Current configuration */
 int p,
 server *Servers,
 int ** preferred_servers
 );

void MST_update_mst
(mpf_t *mst,       /* Stores the response time */
 SQM_instance *Inst,/* Current configuration */
 int p,
 server *Servers,
 mpf_t **f         /* Current matrix of assignations */
 );

void MST_expected_travel_time
(mpf_t t_r,        /* stores the response time */
 SQM_instance *Inst,/* Current configuration */
 int p,
 server *Servers,
 mpf_t **f         /* Current matrix of assignations */
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

void _MST_mpf_init(mpf_t**,int);
void _MST_mpf_init(mpf_t***,int,int);
void _MST_mpf_clear(mpf_t**,int);
void _MST_mpf_clear(mpf_t***,int,int);

#endif

/* eof */
