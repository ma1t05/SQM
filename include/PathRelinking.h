
#ifndef _PATH_RELINKING_H
#define _PATH_RELINKING_H

#include "Reference_Set.h"
#include "PerfectMatching.h"
#include "SQM_heuristic.h"

list<SQM_solution*>* Path_Relinking (SQM_solution&,SQM_solution&);
void SQM_path_relinking(RefSet&,list<SQM_solution*>&);
SQM_solution* SQM_best_solution(list<SQM_solution*>* Solutions);
SQM_solution* SQM_leave_only_the_best(list<SQM_solution*>* Solutions);
void SQM_delete_sols(list<SQM_solution*>* Solutions);
double PR_perfect_matching_cost(SQM_solution&,SQM_solution&);
double SQM_min_cost_pm(RefSet&,SQM_solution&);

enum matching_type {perfect_matching,workload_matching,random_matching,
		    invalid_matching};
matching_type& operator++(matching_type&);
extern int* (*matching_function)(SQM_solution&,SQM_solution&);
matching_type& operator++(matching_type& target);
int* PR_perfect_matching(SQM_solution&,SQM_solution&);
int* PR_random_matching(SQM_solution&,SQM_solution&);
int* PR_workload_matching(SQM_solution&,SQM_solution&);

enum order_type {nearest_first,farthest_first,random_order,invalid_order};
order_type& operator++(order_type&);
extern int* (*order_function)(SQM_solution&,int*,SQM_solution&);
int* PR_processing_order_random(SQM_solution&,int*,SQM_solution&);
int* PR_processing_order_nf(SQM_solution&,int*,SQM_solution&);
int* PR_processing_order_ff(SQM_solution&,int*,SQM_solution&);
extern void (*Improvement_Method)(SQM_solution&);

#endif

/* eof */
