
#ifndef _SQM_HEURISTIC_H
#define _SQM_HEURISTIC_H

#include <utility>
#include <cstdlib>
#include "point.h"
#include "SQM.h"
#include "jarvis.h"
typedef long double num;

struct response_unit {
  int location;
  int past_location;
  double v;
  double beta;
};

typedef struct response_unit response_unit;
response_unit* SQM_heuristic(SQM_instance*,int,double,double);

#endif

/* eof */
