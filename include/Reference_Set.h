#ifndef _REFERENCE_SET_H
#define _REFERENCE_SET_H 1

#include "SQM_Solution.h"
#include "PathRelinking.h"

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
  void recover_garbage (SolList&);
  void sort_by_diversity ();
  bool is_not_in(SQM_solution&);
  SQM_solution* remove(int,int*);
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

#define MAX_ITER 1000

class SubsetControl {
protected:
  int CurrentIter;
  RefSet *rs;
  SolList *pool;
  virtual void Generate_Subsets () = 0;
  virtual void Update(SolList*) = 0;
public:
  RefSet* get_RefSet ();
};

class Static_SC : public SubsetControl {
private:
  int LastRunTime;
  int *LastChange;
  int iNew,jOld;
  int *LocNew,*LocOld;
  void Generate_Subsets ();
  void Update(SolList*);
public:
  Static_SC (int,SolList&);
  ~Static_SC ();
};

class Dynamic_SC : public SubsetControl {
private:
  int StopCondition;
  int LastRunTime;
  int *LastChange;
  int iNew,jOld;
  int *LocNew,*LocOld;
  void Generate_Subsets ();
  void Update(SolList*);
public:
  Dynamic_SC (int,SolList&);
  ~Dynamic_SC ();
};

bool compare_SQMSols(SQM_solution*,SQM_solution*);

class TwoTier_SC : public SubsetControl {
private:
  int b1,b2;
  RefSet *rs2;
  int LastRunTime;
  int *LastChange;
  int iNew,jOld;
  int *LocNew,*LocOld;
  void Generate_Subsets ();
  void Update(SolList*);
  int location(int);
  SQM_solution* Solution(int);
  void Update_diversity ();
public:
  TwoTier_SC (int,int,SolList&);
  ~TwoTier_SC ();
};

void No_Improvement (SQM_solution&);

#endif

/* eof */
