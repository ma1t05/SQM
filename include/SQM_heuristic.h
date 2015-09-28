
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
response_unit* SQM_heuristic(SQM_instance*,int,double,double);

#endif

/* eof */
