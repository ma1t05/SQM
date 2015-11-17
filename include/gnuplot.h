
#ifndef _GNUPLOT_H
#define _GNUPLOT_H 1

#include "SQM.h"

void plot_instance_solution(SQM_instance*,int*,string);
void plot_solution_allocation(SQM_instance*,int,response_unit*,double**,string,string);
void gnuplot_GRASP(char *filename);
void IC_plot_instance(SQM_instance*,int*,string);

#endif

/* eof */
