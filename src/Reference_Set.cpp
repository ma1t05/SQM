
#include "Reference_Set.h"
#include "PathRelinking.h"
#define UNAGGREGATED -1

RefSet::RefSet (int Max) {
  bMax = Max;
  bNow = 0;

  Solutions = new SQM_solution*[bMax];
  for (int i = 0;i < bMax;i++)
    Solutions[i] = NULL;
  loc = new int [bMax];
  ObjVal = new double [bMax];
  DivVal = new double [bMax];
  Hash = new double [bMax];
  NewSolutions = false;

  Calls = 0;
  Adds = 0;
  DupCheck = 0;
  FullDupCheck = 0;
  FullDupFound = 0;
}

RefSet::~RefSet () {
  delete [] Hash;
  delete [] DivVal;
  delete [] ObjVal;
  clean_garbage ();
  delete [] loc;
  for (int i = 0;i < bMax;i++)
    if (Solutions[i] != NULL)
      delete Solutions[i];
  delete [] Solutions;
}

int RefSet::Add (SQM_solution &Sol) {
  log_depth++;
  string tag = log_tag("Reference_Set::Add: ");
  logDebug(cout << tag << "Start" << endl);

  int loc0;
  Adds++;
  if (bNow < bMax) {
    bNow++;
    loc[bNow-1] = bNow-1;
  }
  else {
    /* To avoid segmentation default, DO NOT DELETE the sol, send to grabage */
    garbage.push_back(Solutions[loc[bNow-1]]);
  }
  loc0 = loc[bNow-1];
  Solutions[loc0] = &Sol;

  if (NewRank < bNow)
    for (int i = bNow - 1;i >= NewRank;i--)
      loc[i] = loc[i-1];
  loc[NewRank-1] = loc0;

  ObjVal[loc0] = NewObjVal;
  DivVal[loc0] = NewDivVal;
  Hash[loc0] = NewHash;

  logDebug(cout << tag << "Finish" << endl);
  log_depth--;
  return loc0;
}

bool RefSet::EqualSol (SQM_solution &Sol,int i) {
  int iLoc = loc[i];
  /* Objective value check */
  if (ObjVal[iLoc] != NewObjVal)
    return false;

  /* Hash value check */
  DupCheck++;
  if (Hash[iLoc] != NewHash)
    return false;

  /* Full duplication check */
  FullDupCheck++;
  if (Sol != *Solutions[iLoc])
    return false;
  FullDupFound++;
  return true;
}

int RefSet::TryAdd (SQM_solution &Sol,double value) {
  log_depth++;
  string tag = log_tag("RefSet_Dynamic::Update: ");
  logDebug(cout << tag << "Start" << endl);

  int loc0;
  Calls++;
  NewRank = 0;
  NewObjVal = value;
  /*NewDivVal = Evaluation(Sol); /* Diversification Method */
  logDebug(cout << tag << "bNow = " << bNow << endl);
  if (bNow == 0) {
    NewRank = 1;
    logDebug(cout << tag << "Calculate Hash (bNow = 0)" << endl);
    NewDivVal = -1;
    NewHash = Sol.Hash();
    loc0 = Add (Sol);
    
    logDebug(cout << tag << "Finish" << endl);
    log_depth--;
    return loc0;
  }
  else {
    if (NewObjVal >= worst() && is_full ()) {
      logDebug(cout << tag << "Finish" << endl);
      log_depth--;
      return UNAGGREGATED;
    }
    logDebug(cout << tag << "Calculate Hash (bNow > 0)" << endl);
    NewHash = Sol.Hash();

    if (NewObjVal < best()) {
      NewRank = 1;
    }
    else {
      for (int i = bNow-1;i >= 0; i--) {
	if (EqualSol(Sol,i)) {
	  logDebug(cout << tag << "Finish" << endl);
	  log_depth--;
	  return UNAGGREGATED;
	}
	else if (Obj(i) < NewObjVal){
	  NewRank = i+2;
	  loc0 = Add (Sol);
	  logDebug(cout << tag << "Finish" << endl);
	  log_depth--;
	  return loc0;
	}
      }
      NewRank = 1; /* This case is when the solution gives the same evaluation 
		      as the best solution */
    }
    loc0 = Add (Sol);
  }
  logDebug(cout << tag << "Finish" << endl);
  log_depth--;
  return loc0;
}

void RefSet::clean_garbage () {
  logDebug(cout << "Clean garbage: " << garbage.size() << endl);
  while (!garbage.empty()) {
    delete garbage.back();
    garbage.pop_back();
  }
}

SQM_solution* RefSet::operator[] (int idx) const {
  if (idx < 0 || idx >= bNow) return NULL;
  return Solutions[idx];
}

double RefSet::best () const {
  return ObjVal[loc[0]];
}

double RefSet::worst () const {
  return ObjVal[loc[bNow-1]];
}

double RefSet::Obj (int i) const {
  if (i < 0 || i >= bNow) return -1;
  return ObjVal[loc[i]];
}

int RefSet::size () const {
  return bMax;
}

int RefSet::elements () const {
  return bNow;
}

int RefSet::location (int i) const {
  if (i < 0 || i >= bNow) return -1;
  return loc[i];
}

bool RefSet::is_full () const {
  return (bNow == bMax);
}

int RefSet::get_Calls () const {
  return Calls;
}

int RefSet::get_Adds () const {
  return Adds;
}

int RefSet::get_Checks () const {
  return DupCheck;
}

int RefSet::get_FullCheck () const {
  return FullDupCheck;
}

int RefSet::get_DupFound () const {
  return FullDupFound;
}

/*******************************************************************************
 * Static_SubsetControl
 ******************************************************************************/
RefSet* SubsetControl::get_RefSet () {
  return rs;
}
/*******************************************************************************
 * Static_SC
 ******************************************************************************/
Static_SC::Static_SC (int size,SolList &P) {
  SQM_solution *X;
  int iLoc;

  /* Pool list is sorted */
  pool = &P;

  LastChange = new int [size];
  for (iLoc = 0;iLoc < size;iLoc++)
    LastChange[iLoc] = 0;
  LocNew = new int [size];
  LocOld = new int [size];
  CurrentIter = 0;

  rs = new RefSet (size);

  do {

    /* Transfer elements of the pool to the RefSet */
    while (!pool->empty()) {
      X = pool->front();
      pool->pop_front();
      iLoc = rs->TryAdd(*X,X->get_response_time());
      if (iLoc > 0)
	LastChange[iLoc] = CurrentIter;
      else delete X;
    } while (!rs->is_full());

    /* Divide RefSet in New and Old Solutions */
    iNew = 0;
    jOld = 0;
    for (int i = 0;i < size;i++) {
      iLoc = rs->location(i);
      if (LastChange[iLoc] >= CurrentIter)
	LocNew[iNew++] = iLoc;
      else
	LocOld[jOld++] = iLoc;
    }

    if (iNew == 0) break;
    Generate_Subsets ();
    
    
    CurrentIter++;
  } while (CurrentIter < MAX_ITER);
  if (CurrentIter == MAX_ITER)
    logInfo(cout << "Static_SC: Finish by MAX_ITER" << endl);
  else
    logInfo(cout << "Static_SC: Finish by no new solutions" << endl);

  delete [] LocNew;
  delete [] LocOld;
  delete [] LastChange;
}

void Static_SC::Generate_Subsets () {
  int iLoc,jLoc;
  SQM_solution *X,*Y;
  SolList *Combined_Solutions;
  logDebug(cout << "Start algorithm for Subsets" << endl);

  if (iNew > 1) 
    for (int i = 0;i < iNew;i++) {
      iLoc = LocNew[i];
      X = (*rs)[iLoc];
      for (int j = i+1;j < iNew;j++) {
	jLoc = LocNew[j];
	Y = (*rs)[jLoc];
	Combined_Solutions = Combine_Solutions(*X,*Y);
	Update(Combined_Solutions);
      }
    }

  if (jOld > 0)
    for (int i = 0;i < iNew;i++) {
      iLoc = LocNew[i];
      X = (*rs)[iLoc];
      for (int j = 0;j < jOld;j++) {
	jLoc = LocOld[j];
	Y = (*rs)[jLoc];
	Combined_Solutions = Combine_Solutions(*X,*Y);
	Update(Combined_Solutions);
      }
    }
  
}

void Static_SC::Update (SolList *Sols) {
  SQM_solution *X;
  if (Sols == NULL) return;
  while (!Sols->empty()) {
    X = Sols->front();
    Sols->pop_front();
    Improvement_Method(*X);
    if (X->get_response_time() < rs->worst())
      pool->push_back(X);
    else
      delete X;
  }
  delete Sols;
}

/*******************************************************************************
 * Dynamic_SC
 ******************************************************************************/

Dynamic_SC::Dynamic_SC (int size,SolList &P) {
  logDebug(cout << "RefSet_Dynamic::Dynamic_SC: Start" << endl);
  /* Initialization */
  int iLoc;
  SolList::iterator Z,it;
  SQM_solution *X,*Div;
  pool = &P;

  rs = new RefSet (size);
  CurrentIter = 0;
  LastChange = new int [size];
  for (iLoc = 0;iLoc < size;iLoc++)
    LastChange[iLoc] = 0;

  /* Insert quality solutions */
  do {
    X = pool->front();
    pool->pop_front();
    rs->TryAdd(*X,X->get_response_time());
  } while (rs->elements() < size/2);
 
  /*rs->drop_garbage();*/
  
  /* Insert diverse solutions */
  /* Determine min cost of pm */
  for (Z = pool->begin();Z != pool->end();Z++) {
    X = *Z;
    X->pm_cost = min_cost_pm(*rs,*X);
  }
  
  do {
    Div = pool->front();
    it = pool->begin();
    for (Z = pool->begin();Z != pool->end();Z++)
      if (Div->pm_cost < (**Z).pm_cost) {
	Div = *Z;
	it = Z;
      }
    pool->erase(it);
    rs->TryAdd(*Div,Div->get_response_time());

    if (rs->is_full()) break;

    /* Update Diversity value */
    double min_cost;
    for (Z = pool->begin();Z != pool->end();Z++) {
      X = *Z;
      min_cost = PR_perfect_matching_cost(*Div,*X);
      if (min_cost < X->pm_cost)
	X->pm_cost = min_cost;
    }

  } while (true);
  
  LocNew = new int [rs->elements()];
  LocOld = new int [rs->elements()];
  StopCondition = 0;
  LastRunTime = 0;
  while (StopCondition == 0) {

    /*++SubsetType;*/
    CurrentIter++;

    iNew = 0;
    jOld = 0;
    for (int i = 0;i < rs->elements();i++) {
      iLoc = rs->location(i);
      if (LastChange[iLoc] >= LastRunTime)
	LocNew[iNew++] = iLoc;
      else
	LocOld[jOld++] = iLoc;
    }

    if (iNew == 0)
      break;

    logDebug(cout << "Continue Dynamic_SC with "
	     << 100 * iNew / rs->elements()
	     << "% of new solutions" << endl);

    Generate_Subsets ();

    if (StopCondition > 0) {
      /* Actually no StopCondition while new solutions were found */
    }
    LastRunTime = CurrentIter;
  }
}

Dynamic_SC::~Dynamic_SC () {
  delete [] LocOld;
  delete [] LocNew;
  delete [] LastChange;
  delete rs;
}

double min_cost_pm (RefSet &Sols,SQM_solution &Sol) {
  double cost,min_cost;
  
  min_cost = PR_perfect_matching_cost(*Sols[0],Sol);
  for (int i = 1;i < Sols.elements();i++) {
    cost = PR_perfect_matching_cost(*Sols[i],Sol);
    if (min_cost > cost) min_cost = cost;
  }

  return min_cost;
}

void Dynamic_SC::Generate_Subsets () {
  int iLoc,jLoc;
  SQM_solution *X,*Y;
  SolList *Combined_Solutions;
  logDebug(cout << "Start algorithm for SubsetType1" << endl);
  if (LogDebug) {
    cout << "LocNew :";
    for (int i = 0;i < iNew;i++)
      cout << " " << LocNew[i];
    cout << endl;
    cout << "LocOlc :";
    for (int j = 0;j < jOld;j++)
      cout << " " << LocOld[j];
    cout << endl;
  }
  if (iNew > 1) 
    for (int i = 0;i < iNew;i++) {
      iLoc = LocNew[i];
      if (LastChange[iLoc] < CurrentIter) {
	X = (*rs)[iLoc];
	logDebug(cout << "Combine " << iLoc << " :");
	for (int j = i+1;j < iNew;j++) {
	  jLoc = LocNew[j];
	  if (LastChange[jLoc] < CurrentIter) {
	    Y = (*rs)[jLoc];
	    logDebug(cout << " " << jLoc);
	    /* Create C(X) and execute improvement method */
	    Combined_Solutions = Combine_Solutions(*X,*Y);
	    Update(Combined_Solutions);
	    /* Optional Check: if LastChange[iLoc] == CurrentIter, 
	       then can jump to the end of "i loop" to pick up the next i,
	       and generate fewer solutions. */
	    if (LastChange[iLoc] == CurrentIter)
	      break;
	  }
	}
	logDebug(cout << endl);
	rs->clean_garbage();
      }
    }
  if (jOld > 0)
    for (int i = 0;i < iNew;i++) {
      iLoc = LocNew[i];
      if (LastChange[iLoc] < CurrentIter) {
	X = (*rs)[iLoc];
	logDebug(cout << "Combine " << iLoc << " :");
	for (int j = 0;j < jOld;j++) {
	  jLoc = LocOld[j];
	  if (LastChange[jLoc] < CurrentIter) {
	    Y = (*rs)[jLoc];
	    logDebug(cout << " " << jLoc);
	    /* Create C(X) and execute improvement method */
	    Combined_Solutions = Combine_Solutions(*X,*Y);
	    Update(Combined_Solutions);
	    /* Optional Check: if LastChange[iLoc] == CurrentIter, 
	       then can jump to the end of "i loop" to pick up the next i,
	       and generate fewer solutions. */
	    if (LastChange[iLoc] == CurrentIter)
	      break;
	  }
	}
	logDebug(cout << endl);
	rs->clean_garbage();
      }
    }
}

void Dynamic_SC::Update (SolList *NewSols) {
  int iLoc;
  SolList::iterator Z;
  SQM_solution *X;
  if (NewSols == NULL) return;
  while (!NewSols->empty()) {
    X = NewSols->front();
    NewSols->pop_front();
    Improvement_Method(*X);
    iLoc = rs->TryAdd(*X,X->get_response_time());
    if (iLoc == UNAGGREGATED)
      delete X;
    else
      LastChange[iLoc] = CurrentIter;
  }
}

bool compare_SQMSols(SQM_solution *first,SQM_solution *second) {
  return (first->get_response_time() < second->get_response_time());
}
