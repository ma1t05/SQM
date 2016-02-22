
#include "Reference_Set.h"
#include "PathRelinking.h"


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
    return 1;
  }
  else {
    if (NewObjVal >= worst() && is_full ()) {
      logDebug(cout << tag << "Finish" << endl);
      log_depth--;
      return 0;
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
	  return 0;
	}
	else if (Obj(i) < NewObjVal){
	  NewRank = i+2;
	  loc0 = Add (Sol);
	  logDebug(cout << tag << "Finish" << endl);
	  log_depth--;
	  return NewRank;
	}
      }
      NewRank = 1; /* This case is when the solution gives the same evaluation 
		      as the best solution */
    }
    loc0 = Add (Sol);
    /*LastChange[loc0] = NowTime;*/
  }
  logDebug(cout << tag << "Finish" << endl);
  log_depth--;
  return NewRank;
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
 * Dynamic_SubsetControl
 ******************************************************************************/

Dynamic_SubsetControl::Dynamic_SubsetControl (RefSet &EliteSols) {
  /* Initialization */
  int iLoc;
  logDebug(cout << "RefSet_Dynamic::Dynamic_SubsetControl: Start" << endl);
  LastChange = new int [EliteSols.size ()];
  for (iLoc = 0;iLoc < EliteSols.size();iLoc++)
    LastChange[iLoc] = 0;
  NowTime = 0;
  StopCondition = 0;
  LocNew = new int [EliteSols.elements()];
  LocOld = new int [EliteSols.elements()];
  LastRunTime = 0;
  while (StopCondition == 0) {

    /*++SubsetType;*/
    NowTime++;

    iNew = 0;
    jOld = 0;
    for (int i = 0;i < EliteSols.elements();i++) {
      iLoc = EliteSols.location(i);
      if (LastChange[iLoc] >= LastRunTime)
	LocNew[iNew++] = iLoc;
      else
	LocOld[jOld++] = iLoc;
    }

    if (iNew == 0) {
      delete [] LocNew;
      delete [] LocOld;
      return;
    }
    logDebug(cout << "Continue Dynamic_SubsetControl with "
	     << 100 * iNew / EliteSols.elements()
	     << "% of new solutions" << endl);

    algorithm_for_SubsetType1 (EliteSols);

    if (StopCondition > 0) {
      /* Actually no StopCondition while new solutions were found */
    }
    LastRunTime = NowTime;
  }
}

Dynamic_SubsetControl::~Dynamic_SubsetControl () {
  delete [] LastChange;
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

void Dynamic_SubsetControl::algorithm_for_SubsetType1 (RefSet &Solutions) {
  int iLoc,jLoc;
  SQM_solution *X,*Y;
  SolList *Combined_Solutions;
  logDebug(cout << "Start algorithm for SubsetType1" << endl);
  if (iNew > 1) 
    for (int i = 0;i < iNew;i++) {
      iLoc = LocNew[i];
      if (LastChange[iLoc] < NowTime) {
	X = Solutions[iLoc];
	for (int j = i+1;j < iNew;j++) {
	  jLoc = LocNew[j];
	  if (LastChange[jLoc] < NowTime) {
	    Y = Solutions[jLoc];
	    /* Create C(X) and execute improvement method */
	    Combined_Solutions = Combine_Solutions(*X,*Y);
	    Update(Combined_Solutions,Solutions);
	    /* Optional Check: if LastChange[iLoc] == NowTime, 
	       then can jump to the end of "i loop" to pick up the next i,
	       and generate fewer solutions. */
	    if (LastChange[iLoc] == NowTime)
	      break;
	  }
	}
	Solutions.clean_garbage();
      }
    }
  if (jOld > 0)
    for (int i = 0;i < iNew;i++) {
      iLoc = LocNew[i];
      if (LastChange[iLoc] < NowTime) {
	X = Solutions[iLoc];
	for (int j = 0;j < jOld;j++) {
	  jLoc = LocOld[j];
	  if (LastChange[jLoc] < NowTime) {
	    Y = Solutions[jLoc];
	    /* Create C(X) and execute improvement method */
	    Combined_Solutions = Combine_Solutions(*X,*Y);
	    Update(Combined_Solutions,Solutions);
	    /* Optional Check: if LastChange[iLoc] == NowTime, 
	       then can jump to the end of "i loop" to pick up the next i,
	       and generate fewer solutions. */
	    if (LastChange[iLoc] == NowTime)
	      break;
	  }
	}
	Solutions.clean_garbage();
      }
    }
}

void Dynamic_SubsetControl::Update (SolList *NewSols,RefSet &Sols) {
  SolList::iterator Z;
  SQM_solution *X;
  while (!NewSols->empty()) {
    X = NewSols->front();
    NewSols->pop_front();
    Improvement_Method(*X);
    if (!Sols.TryAdd(*X,X->get_response_time())) delete X;
  }
}
