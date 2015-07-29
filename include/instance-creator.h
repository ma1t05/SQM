
#ifndef _INSTANCE_CREATOR_H
#define _INSTANCE_CREATOR_H 1

#include "point.h"

struct SQM_instance {
  int M; /* Number of customers */
  int N; /* Number of potencial sites to locate a server */
  point* V; /* Set of demand points */
  point* W; /* Set of potencial locations sites */
};

typedef struct SQM_instance instance;



#endif

/* eof */
