
#ifndef _INSTANCE_CREATOR_H
#define _INSTANCE_CREATOR_H 1

#include <cstdlib>
#include <ctime>
#include "SQM.h"

#define MIN_X 0.0
#define MAX_X 1024.0
#define MIN_Y 0.0
#define MAX_Y 1024.0

SQM_instance* IC_create_instance (int n,int m);
SQM_instance* IC_read_instance (string,string);
void IC_write_instance (SQM_instance*,string,string);
void IC_plot_instance(SQM_instance*I,int*,string);

#endif

/* eof */
