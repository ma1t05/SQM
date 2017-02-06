
#ifndef _GNUPLOT_H
#define _GNUPLOT_H 1

#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "SQM_Solution.h"

void plot_instance_solution(SQM_instance&,int*,string);
void plot_solution_allocation(SQM_solution*,double**,string,string);
void gnuplot_GRASP(char *filename);
void gnuplot_solution(SQM_solution*,int);
void gnuplot_Path_Relinking(SQM_solution&,int*,SQM_solution&,string);
void gnuplot_write_points_file(char*,point*,int);
void gnuplot_sets(FILE*);
void gnuplot_unsets(FILE*);

#endif

/* eof */
