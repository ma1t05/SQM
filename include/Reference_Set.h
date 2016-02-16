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

typedef double (*EvaluationMethod)(SQM_solution&);
typedef void (*ImprovementMethod)(SQM_solution&);
typedef std::list<SQM_solution*> SolList;

class Reference_Set {
private:
  int bMax;
  int bNow;
  int Adds; /* Count the additions to Reference Set */
  /* Count the Checks */
  int DupCheck;
  int FullDupCheck;
  int FullDupFound;
  int *loc;
  double *E,*Hash;
  SolList garbage;
protected:
  int NewRank;
  int RefSetCall;
  double E0,Hash0;
  SQM_solution **Solutions;
  int Add(SQM_solution&);
  bool EqualSol(SQM_solution&,int);
public:
  Reference_Set (int);
  ~Reference_Set ();
  virtual bool Update (SQM_solution&) = 0;
  void clean_garbage ();
  bool it_is_full () {return bNow == bMax;};
  int get_Adds () const {return Adds;};
  int get_Checks () const {return DupCheck;};
  int get_FullCheck () const {return FullDupCheck;};
  int get_DupFound () const {return FullDupFound;};
  int get_elements () const {return bNow;};
  int location (int i) const {
    return ((i >= 0 && i < bNow) ? loc[i] : -1);
  };
  double evaluation (int i) const {
    return ((i >= 0 && i < bNow) ? E[loc[i]] : -1);
  };
};

class RefSet : public Reference_Set {
private:
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
  ImprovementMethod Improvement;
  EvaluationMethod Evaluation;
public:
  RefSet (int,ImprovementMethod,EvaluationMethod);
  ~RefSet ();
  bool Update (SQM_solution&);
  void Call_Improvement(SQM_solution&);
  void SubsetControl ();
  double best () const {return evaluation(location(0));};
  double worst () const {return evaluation(location(get_elements()-1));};
  SQM_solution* best_sol () const {return Solutions[location(0)];};
  SQM_solution* worst_sol () const {
    return Solutions[location(get_elements()-1)];
  };
  SQM_solution* get_sol (int i) const {
    return ((i >= 0 && i < get_elements()) ? Solutions[location(i)] : NULL);
  };
  int get_Calls () const {return RefSetCall;};
  /*int get_bNow () const {return bNow;};*/
};

typedef double (*DiversificationMethod)(RefSet&,SQM_solution&);

class RefSet_2 : public RefSet {
private:
  SQM_solution **quality;
  SQM_solution **diverse;
  DiversificationMethod Diversification;
public:
  RefSet_2 (int,int,ImprovementMethod,EvaluationMethod,DiversificationMethod);
  ~RefSet_2 ();
  
};

double get_response_time (SQM_solution&);
double get_perfect_matching_cost (SQM_solution&);

#endif

/* eof */
