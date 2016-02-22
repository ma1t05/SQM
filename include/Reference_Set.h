#ifndef _REFERENCE_SET_H
#define _REFERENCE_SET_H 1

#include <list>
#include "SQM_Solution.h"

typedef std::list<SQM_solution*> SolList;

class RefSet {
private:
  int bMax;
  int bNow;
  SQM_solution **Solutions;
  int *loc;
  double *ObjVal;
  double *DivVal;
  double *Hash;
  SolList garbage;
  bool NewSolutions;
protected:
  int Calls;
  int Adds;
  int DupCheck;
  int FullDupCheck;
  int FullDupFound;

  int NewRank;
  double NewObjVal;
  double NewDivVal;
  double NewHash;
  /* Pass evaluation and hash in NewObjVal & NewHash */
  int Add(SQM_solution&);
  bool EqualSol(SQM_solution&,int);
public:
  RefSet (int);
  ~RefSet ();
  int TryAdd (SQM_solution&,double);
  void clean_garbage ();
  SQM_solution* operator[](int) const;
  double best () const;
  double worst () const;
  bool is_full () const;
};

typedef double (*EvaluationMethod)(SQM_solution&);
typedef void (*ImprovementMethod)(SQM_solution&);
extern SolList* (*Combine_Solutions)(SQM_solution&,SQM_solution&);
class Reference_Set {
public:
  void Update (SolList*);
  virtual void SubsetControl () = 0;
  virtual double min_cost_pm (SQM_solution&);
  bool it_is_full () {return bNow == bMax;};
  int get_Calls () const {return Calls;};
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
  double best () const {return evaluation(0);};
  double worst () const {return evaluation(get_elements()-1);};
};

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

class RefSet_dynamic : public Reference_Set {
private:
  void algorithm_for_SubsetType1 ();
  void algorithm_for_SubsetType2 ();
  void algorithm_for_SubsetType3 ();
  void algorithm_for_SubsetType4 ();
  ImprovementMethod Improvement;
  EvaluationMethod Evaluation;
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
  RefSet_dynamic (int,ImprovementMethod,EvaluationMethod);
  ~RefSet_dynamic ();
  bool Update (SQM_solution&);
  void SubsetControl ();
  double min_cost_pm (SQM_Solution&);
  void Call_Improvement(SQM_solution&); /* Is required? */
  SQM_solution* best_sol () const {return Solutions[location(0)];};
  SQM_solution* worst_sol () const {
    return Solutions[location(get_elements()-1)];
  };
  SQM_solution* get_sol (int i) const {
    return ((i >= 0 && i < get_elements()) ? Solutions[location(i)] : NULL);
  };
  /*int get_bNow () const {return bNow;};*/
};

class RefSet_2Tier : public Reference_Set {
private:
  int bMax2;
  int bNow2;
  int Adds2;
  int *loc;
  double *D,*Hash2;
protected:
  SQM_solution **Diverse;
  void reorder_diverse_set ();
  bool quality_update (SQM_solution&);
  bool diverse_update (SQM_solution&)
public:
  RefSet_2Tier (int,int,ImprovementMethod,EvaluationMethod)
  ~RefSet_2Tier ();
  bool Update (SQM_solution&);
  void SubsetControl ();
  double min_cost_pm (SQM_solution&);
};

double get_response_time (SQM_solution&);
double get_perfect_matching_cost (SQM_solution&);

#endif

/* eof */
