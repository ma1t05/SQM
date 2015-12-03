
#ifndef _PATH_RELINKING_H
#define _PATH_RELINKING_H

#include <list>
#include "SQM_Solution.h"

list<SQM_solution*>* Path_Relinking (SQM_solution*,SQM_solution*);
SQM_solution* SQM_path_relinking(list<SQM_solution*>*);
SQM_solution* SQM_best_solution(list<SQM_solution*>* Solutions);
SQM_solution* SQM_leave_only_the_best(list<SQM_solution*>* Solutions);
void SQM_delete_sols(list<SQM_solution*>* Solutions);
double PR_perfect_matching_cost(SQM_solution *X,SQM_solution *Y);
double SQM_min_cost_pm(list<SQM_solution*>*,SQM_solution*);

#endif

/* eof */
