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
  double Obj (int) const;
  int size () const;
  int elements () const;
  int location (int) const;
  bool is_full () const;
  
  int get_Calls () const;
  int get_Adds () const;
  int get_Checks () const;
  int get_FullCheck () const;
  int get_DupFound () const;
};

extern void (*Improvement_Method)(SQM_solution&);
extern SolList* (*Combine_Solutions)(SQM_solution&,SQM_solution&);
double min_cost_pm (RefSet&,SQM_solution&);

class Dynamic_SubsetControl {
private:
  void algorithm_for_SubsetType1 (RefSet&);
  void Update(SolList*,RefSet&);
public:
  Dynamic_SubsetControl(RefSet&);
  ~Dynamic_SubsetControl();
  int NowTime;
  int StopCondition;
  int LastRunTime;
  int *LastChange;
  int iNew,jOld;
  int *LocNew,*LocOld;
};

#endif

/* eof */
