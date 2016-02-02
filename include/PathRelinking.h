
#ifndef _PATH_RELINKING_H
#define _PATH_RELINKING_H

#include <list>
#include "SQM_Solution.h"
#include "PerfectMatching.h"
#include "SQM_heuristic.h"

enum Subset {two_element,   /* all 2-element subsets. */
	     three_element, /* 3-element subsets derived from the 2-element 
			       subsets by augmenting each 2-element subset to 
			       include the best solution not in this subset */
	     four_element,  /* 4-element subsets derived from the 3-element 
			       subsets by augmenting each 3-element subset to
			       include the best solution not in this subset */
	     best_i,        /* the subsets consisting in the best i elements,
			       for i = 5 to bNow */
	     invalid_subset};
Subset& operator++(Subset&);

class RefSet {
 private:
  int bMax;
  int bNow;
  int NewRank;
  int RefSetCall;
  int RefSetAdd;
  int DupCheck;
  int FullDupCheck;
  int FullDupFound;
  int *loc;
  double E0,Hash0;
  double *E,*Hash;
  SQM_solution **Solutions;
  void Add(SQM_solution&);
  void SubsetControl ();
  void algorithm_for_SubsetType1 ();
  void algorithm_for_SubsetType2 ();
  void algorithm_for_SubsetType3 ();
  void algorithm_for_SubsetType4 ();
protected:
  int NowTime;
  Subset SubsetType;
  int StopCondition;
  int *LastChange;
  int *LastRunTime;
  int iNew,jOld;
  int *LocNew;
  int *LocOld;
 public:
  RefSet (int);
  ~RefSet ();
  void Update (SQM_solution&);
  double best () const {return E[loc[0]];};
  double worst () const {return E[loc[bNow-1]];};
  SQM_solution* best_sol () const {return Solutions[loc[0]];};
  SQM_solution* worst_sol () const {return Solutions[loc[bNow-1]];};
  int Calls () const {return RefSetCall;};
  int Adds () const {return RefSetAdd;};
  int Checks () const {return DupCheck;};
  int FullCheck () const {return FullDupCheck;};
  int DupFound () const {return FullDupFound;};
};

list<SQM_solution*>* Path_Relinking (SQM_solution*,SQM_solution*);
SQM_solution* SQM_path_relinking(list<SQM_solution*>*);
SQM_solution* SQM_best_solution(list<SQM_solution*>* Solutions);
SQM_solution* SQM_leave_only_the_best(list<SQM_solution*>* Solutions);
void SQM_delete_sols(list<SQM_solution*>* Solutions);
double PR_perfect_matching_cost(SQM_solution *X,SQM_solution *Y);
double SQM_min_cost_pm(list<SQM_solution*>*,SQM_solution*);

enum matching_type {perfect_matching,workload_matching,random_matching,
		    invalid_matching};
matching_type& operator++(matching_type&);
extern int* (*matching_function)(SQM_solution*,SQM_solution*);
matching_type& operator++(matching_type& target);
int* PR_run_perfect_matching(SQM_solution*,SQM_solution*);
int* PR_random_matching(SQM_solution*,SQM_solution*);
int* PR_workload_matching(SQM_solution*,SQM_solution*);

enum order_type {nearest_first,farthest_first,random_order,invalid_order};
order_type& operator++(order_type&);
extern int* (*order_function)(SQM_solution*,int*,SQM_solution*);
int* PR_processing_order_random(SQM_solution*,int*,SQM_solution*);
int* PR_processing_order_nf(SQM_solution*,int*,SQM_solution*);
int* PR_processing_order_ff(SQM_solution*,int*,SQM_solution*);

#endif

/* eof */
