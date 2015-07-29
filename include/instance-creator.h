
#ifndef _INSTANCE_CREATOR_H
#define _INSTANCE_CREATOR_H 1

#include <cstdlib>
#include "point.h"

#define MIN_X 0.0
#define MAX_X 1024.0
#define MIN_Y 0.0
#define MAX_Y 1024.0

struct SQM_instance {
  int M; /* Number of customers */
  int N; /* Number of potencial sites to locate a server */
  point* V; /* Set of demand points */
  point* W; /* Set of potencial locations sites */
};

typedef struct SQM_instance SQM_instance;

SQM_instance* IC_create_instance (int n,int m);
SQM_instance* IC_read_instance (string,string);
void IC_write_instance (SQM_instance*,string);

#endif

/* eof */
