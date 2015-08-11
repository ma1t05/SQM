
#ifndef _SQM_HEURISTIC_H
#define _SQM_HEURISTIC_H

#include <utility>
#include <cstdlib>
#include "point.h"
#include "jarvis.h"

struct response_unit {
  int location;
  double v;
  double beta;
};

struct SQM_instance {
  int M; /* Number of customers */
  int N; /* Number of potencial sites to locate a server */
  point* V; /* Set of demand points */
  point* W; /* Set of potencial locations sites */
};

typedef struct response_unit response_unit;
typedef struct SQM_instance SQM_instance;

#endif

/* eof */
