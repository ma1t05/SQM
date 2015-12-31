
#ifndef _GNUPLOT_H
#define _GNUPLOT_H 1

#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "SQM_Solution.h"

void plot_instance_solution(SQM_instance*,int*,string);
void plot_solution_allocation(SQM_solution*,double**,string,string);
void gnuplot_GRASP(char *filename);
void IC_plot_instance(SQM_instance*,int*,string);
void gnuplot_solution(SQM_solution*,int);

#endif

/* eof */
