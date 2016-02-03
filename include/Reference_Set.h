#ifndef _REFERENCE_SET_H
#define _REFERENCE_SET_H 1

#include <list>
#include "SQM_Solution.h"

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
  bool Update (SQM_solution&);
  void SubsetControl ();
  double best () const {return E[loc[0]];};
  double worst () const {return E[loc[bNow-1]];};
  double evaluation (int i) const {return ((i >= 0 && i < bNow) ? E[loc[i]] : -1);};
  SQM_solution* best_sol () const {return Solutions[loc[0]];};
  SQM_solution* worst_sol () const {return Solutions[loc[bNow-1]];};
  SQM_solution* get_sol (int i) const {
    return ((i >= 0 && i < bNow) ? Solutions[loc[i]] : NULL);
  };
  int get_elements () const {return bNow;};
  int get_Calls () const {return RefSetCall;};
  int get_Adds () const {return RefSetAdd;};
  int get_Checks () const {return DupCheck;};
  int get_FullCheck () const {return FullDupCheck;};
  int get_DupFound () const {return FullDupFound;};
  int get_bNow () const {return bNow;};
};

#endif

/* eof */
