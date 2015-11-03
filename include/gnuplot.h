
#ifndef _GNUPLOT_H
#define _GNUPLOT_H 1

#include <cstdlib>
#include "SQM.h"

void plot_instance_solution(SQM_instance*,int*,string);
void plot_solution_allocation(SQM_instance*,int,response_unit*,double**,string,string);

#endif

/* eof */
