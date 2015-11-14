
#ifndef _SQM_H
#define _SQM_H

#include <fstream>
#include "point.h"

#define EPSILON 0.001
#define TIME_MAX 1200.0
/* */
#define MINS_PER_BLOCK 60
#define BLOCKS_PER_HORIZON 24

extern std::ofstream LogFile;
extern std::ofstream results;
extern std::ofstream dat;

#endif

/* eof */

