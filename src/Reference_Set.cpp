
#include "Reference_Set.h"
#include "PathRelinking.h"
#define UNAGGREGATED -1
#define INVALID_POSITION -1

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

void RefSet::recover_garbage (SolList &pool) {
  logDebug(cout << "Recover garbage: " << garbage.size() << endl);
  while (!garbage.empty()) {
    pool.push_front(garbage.back());
    garbage.pop_back();
  }
}

void RefSet::sort_by_diversity () {
  double *Div;

  /* Update DivVal */
  Div = new double [bNow];
  for (int i = 0;i < bNow;i++)
    Div[i] = -Solutions[i]->pm_cost;

  /* Sort by DivVal */
  sort_dist(bNow,Div,loc);
  delete [] Div;

  for (int i = 0;i < bNow;i++)
    ObjVal[i] = -Solutions[i]->pm_cost;

}

bool RefSet::is_not_in (SQM_solution& Sol) {
  log_depth++;
  string tag = log_tag("RefSet::is_not_in: ");
  int iLoc;
  
  NewDivVal = Sol.pm_cost;
  if (NewDivVal < worst()) {
    logDebug(cout << tag << "Finish" << endl);
    log_depth--;
    return true;
  }
  logDebug(cout << tag << "Calculate Hash" << endl);
  NewHash = Sol.Hash();

  for (int i = 0;i < bNow;i++) {
    iLoc = loc[i];
    if (&Sol == Solutions[iLoc]) {
      logDebug(cout << tag << "Finish - Solution in RefSet at location " << i
	      << endl);
      log_depth--;
      return false;
    }
  }

  SolList::iterator Z;
  for (Z = garbage.begin();Z != garbage.end();Z++)
    if (&Sol == *Z) {
      logDebug(cout << tag << "Solution was in garbage (size:" << garbage.size()
	      << ") " << Sol.pm_cost << " : " << (**Z).pm_cost << endl);
      garbage.erase(Z);
      log_depth--;
      return true;
    }
  log_depth--;
  return true;
  
}

SQM_solution* RefSet::remove (int idx,int *LastChange) { 
  int iLoc,lastLoc;
  SQM_solution *Sol;
  if (idx < 0 || idx >= bNow) return NULL;
  iLoc = loc[idx];
  lastLoc = bNow-1;
  Sol = Solutions[iLoc];
  if (iLoc != lastLoc) {
    /* Move information from lastLoc to iLoc */
    Solutions[iLoc] = Solutions[lastLoc];
    ObjVal[iLoc] = ObjVal[lastLoc];
    DivVal[iLoc] = DivVal[lastLoc];
    Hash[iLoc] = Hash[lastLoc];
    if (LastChange != NULL)
      LastChange[iLoc] = LastChange[lastLoc];
    Solutions[lastLoc] = NULL;
    int i = 0;
    while (loc[i] != lastLoc) i++;
    loc[i] = iLoc;
    bNow--;
  }

  for (int i = idx;i < bNow;i++)
    loc[i] = loc[i+1];

  return Sol;
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
 * SubsetControl
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
      if (iLoc != UNAGGREGATED)
	LastChange[iLoc] = CurrentIter;
      else delete X;
    }
    rs->clean_garbage();

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
  } while (++CurrentIter < MAX_ITER);

  if (CurrentIter == MAX_ITER)
    logInfo(cout << "Static_SC: Finish by MAX_ITER" << endl);
  else
    logInfo (cout << "Static_SC: Finish by no new solutions at " << CurrentIter
	     << " iterations"  << endl);

}

Static_SC::~Static_SC () { 
  delete rs;
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
  
  /* Clean pool */
  while (!pool->empty()) {
    delete pool->front();
    pool->pop_front();
  }

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

/*******************************************************************************
 * SubsetControl
 ******************************************************************************/

TwoTier_SC::TwoTier_SC (int B1,int B2,SolList &P) {
  log_depth++;
  string tag = log_tag("TwoTier_SC::constructor: ");
  logInfo(cout << tag << "Start" << endl);

  /* Initialization */
  int iLoc;
  SolList::iterator Z,it;
  SQM_solution *X,*Div;
  b1 = B1;
  b2 = B2;
  pool = &P;

  logInfo(cout << tag << "Create RefSet with quality solutions" << endl);
  rs = new RefSet (b1);
  CurrentIter = 0;
  LastChange = new int [b1+b2];
  for (iLoc = 0;iLoc < b1+b2;iLoc++)
    LastChange[iLoc] = 0;

  /* Insert quality solutions */
  do {
    X = pool->front();
    pool->pop_front();
    rs->TryAdd(*X,X->get_response_time());
  } while (rs->elements() < b1);
  if (LogInfo) {
    cout << "The best " << b1 << " response times" << endl;
    for (int i = 0;i < b1;i++)
      cout << (*rs)[i]->get_response_time() << endl;
  }
  
  logInfo(cout << tag << "Create RefSet with diverity solutions" << endl);
  /* Insert diverse solutions */
  rs2 = new RefSet (b2);
  Update_diversity ();
  if (LogInfo) {
    cout << "The " << b2 << " diverse values" << endl;
    for (int i = 0;i < b2;i++)
      cout << (*rs2)[i]->pm_cost << endl;
  }

  LocNew = new int [b1+b2];
  LocOld = new int [b1+b2];
  LastRunTime = 0;
  while (true) {

    /*++SubsetType;*/
    CurrentIter++;

    iNew = 0;
    jOld = 0;
    for (int i = 0;i < b1+b2;i++) {
      iLoc = location(i);
      if (LastChange[iLoc] >= LastRunTime)
	LocNew[iNew++] = iLoc;
      else
	LocOld[jOld++] = iLoc;
    }

    if (iNew == 0)
      break;

    logInfo(cout << tag << "Continue TwoTier_SC with "
	     << 100 * iNew / rs->elements()
	     << "% of new solutions" << endl);

    Generate_Subsets ();

    LastRunTime = CurrentIter;
  }
  logInfo(cout << tag << "Finish" << endl);
  log_depth--;
}

TwoTier_SC::~TwoTier_SC ( ) {
  delete [] LocOld;
  delete [] LocNew;
  delete rs2;
  delete [] LastChange;
  delete rs;
}

void TwoTier_SC::Generate_Subsets () {
  int iLoc,jLoc;
  SQM_solution *X,*Y;
  SolList *Combined_Solutions;

  if (LogInfo) {
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
	X = Solution(iLoc);
	for (int j = i+1;j < iNew;j++) {
	  jLoc = LocNew[j];
	  if (LastChange[jLoc] < CurrentIter) {
	    Y = Solution(jLoc);
	    logInfo(cout << "Combine " << iLoc << " : " << jLoc << endl);
	    cout << "with response times" << endl;
	    cout << X->get_response_time() << endl
		 << Y->get_response_time() << endl;
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
	rs->recover_garbage(*pool); 
	Update_diversity ();
      }
    }
  if (jOld > 0)
    for (int i = 0;i < iNew;i++) {
      iLoc = LocNew[i];
      if (LastChange[iLoc] < CurrentIter) {
	X = Solution(iLoc);
	for (int j = 0;j < jOld;j++) {
	  jLoc = LocOld[j];
	  if (LastChange[jLoc] < CurrentIter) {
	    Y = Solution(jLoc);
	    logInfo(cout << "Combine " << iLoc << " : " << jLoc << endl);
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
	rs->recover_garbage(*pool);
	Update_diversity ();
      }
    }
}

void TwoTier_SC::Update (SolList *NewSols) {
  int iLoc;
  SolList::iterator Z;
  SQM_solution *X;

  if (NewSols == NULL) return;
  log_depth++;
  string tag = log_tag("TwoTier_SC::Update: ");
  logInfo(cout << tag << "Start" << endl);

  while (!NewSols->empty()) {
    X = NewSols->front();
    NewSols->pop_front();
    Improvement_Method(*X);
    iLoc = rs->TryAdd(*X,X->get_response_time());
    if (iLoc == UNAGGREGATED)
      pool->push_front(X);
    else
      LastChange[iLoc] = CurrentIter;
  }
  delete NewSols;
  logInfo(cout << tag << "Finish" << endl);
  log_depth--;
}

int TwoTier_SC::location (int i) {
  if (i < 0 || i >= b1+b2) return INVALID_POSITION;
  if (i < b1)
    return rs->location(i);
  else 
    return b1 + rs2->location(i-b1);
}

SQM_solution* TwoTier_SC::Solution (int iLoc) {
  if (iLoc < 0 || iLoc >= b1+b2) {
    cout << "Location " << iLoc << " out of range" << endl;
    return NULL;
  }
  if (iLoc < b1) {
    cout << "Solution in RefSet1 in location " << iLoc << endl;
    return (*rs)[iLoc];
  }
  else {
    cout << "Solution in RefSet2 in location " << iLoc << endl;
    return (*rs2)[iLoc-b1];
  }
}

void TwoTier_SC::Update_diversity () {
  if (pool->empty()) return;
  log_depth++;
  string tag = log_tag("TwoTier_SC::Update_diversity: ");
  logInfo(cout << tag << "Start" << endl);

  int iLoc;
  int adds;
  SQM_solution *DivSol,*Sol;
  SolList::iterator DivIt,Z;

  /* Inser diverse Solution to pool */
  for (int i = b2;i > 0;i--) {
    iLoc = rs2->location(i-1);
    if (iLoc != -1) {
      Sol = (*rs2)[iLoc];
      pool->push_front(Sol);
    }
  }

  /* Update min cost of pm */
  for (Z = pool->begin();Z != pool->end();Z++)
    (**Z).pm_cost = min_cost_pm(*rs,**Z);

  adds = 0;
  do {
    /* sort rs2 */
    rs2->sort_by_diversity();

    DivSol = pool->front();
    DivIt = pool->begin();
    for (Z = pool->begin();Z != pool->end();Z++)
      if (DivSol->pm_cost < (**Z).pm_cost) {
	DivSol = *Z;
	DivIt = Z;
      }
    pool->erase(DivIt);

    adds++;
    if (rs2->is_not_in(*DivSol)) {
      iLoc = rs2->TryAdd(*DivSol,-DivSol->pm_cost);
      LastChange[b1+iLoc] = CurrentIter;
      logDebug(cout << "Add diverse sol " << DivSol->pm_cost << " location " 
	       << b1+iLoc << endl);
    }
      
    if (adds == b2) break;

    /* Update Diversity value */
    double min_cost;
    for (Z = pool->begin();Z != pool->end();Z++) {
      Sol = *Z;
      min_cost = PR_perfect_matching_cost(*DivSol,*Sol);
      if (min_cost < Sol->pm_cost)
	Sol->pm_cost = min_cost;
    }

  } while (!pool->empty());
  SolList *tmp_list = new SolList;
  rs2->recover_garbage (*tmp_list); /* Double free or corruption */
  delete tmp_list;

  /* Optional: Clean pool */
  while (!pool->empty()) {
    delete pool->front();
    pool->pop_front();
  }
  
  logInfo(cout << tag << "Finish" << endl);
  log_depth--;
}

bool compare_SQMSols(SQM_solution *first,SQM_solution *second) {
  return (first->get_response_time() < second->get_response_time());
}

void No_Improvement (SQM_solution &Sol) {};
