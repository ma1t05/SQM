
#include "Reference_Set.h"
#include "PathRelinking.h"

Reference_Set::Reference_Set (int Max) : bMax(Max) {
  bNow = 0;

  DupCheck = 0;
  FullDupCheck = 0;
  FullDupFound = 0;

  RefSetCall = 0;
  Adds = 0;

  loc = new int [bMax];
  E = new double [bMax];
  Hash = new double [bMax];
  Solutions = new SQM_solution*[bMax];
}

Reference_Set::~Reference_Set () {
  cout << "Delete Reference_Set arrays!" << endl;
  for (int i = 0;i < bNow;i++) {
    delete Solutions[loc[i]];
  }
  delete [] Solutions;
  delete [] Hash;
  delete [] E;
  delete [] loc;
}

bool Reference_Set::EqualSol (SQM_solution &Sol,int i) {
  int iLoc = location(i);
  /* Objective value check */
  if (E[iLoc] != E0)
    return false;

  /* Hash value check */
  DupCheck++;
  if (Hash[iLoc] != Hash0)
    return false;

  /* Full duplication check */
  FullDupCheck++;
  if (Sol != *Solutions[iLoc])
    return false;
  FullDupFound++;
  return true;
}

int Reference_Set::Add (SQM_solution &Sol) {
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
  if (NewRank < bNow)
    for (int i = bNow - 1;i >= NewRank;i--)
      loc[i] = loc[i-1];
  loc[NewRank-1] = loc0;

  Solutions[loc0] = &Sol;
  Hash[loc0] = Hash0;
  E[loc0] = E0;
  /* LastChange[loc0] = NowTime; */
  logDebug(cout << tag << "Finish" << endl);
  log_depth--;
  return loc0;
}

void Reference_Set::clean_garbage () {
  logDebug(cout << "Clean garbage: " << garbage.size() << endl);
  while (!garbage.empty()) {
    delete garbage.back();
    garbage.pop_back();
  }
}

RefSet::RefSet (int Max,ImprovementMethod ImprovementFunction,
		EvaluationMethod EvaluationFunction) : Reference_Set(Max) {
  LastChange = new int [Max];
  LastRunTime = new int [invalid_subset];
  Improvement = ImprovementFunction;
  Evaluation = EvaluationFunction;
}

RefSet::~RefSet () {
  clean_garbage();
  delete [] LastRunTime;
  delete [] LastChange;
}

bool RefSet::Update (SQM_solution &Sol) {
  int loc0;
  log_depth++;
  string tag = log_tag("RefSet::Update: ");
  logDebug(cout << tag << "Start" << endl);
  RefSetCall++;
  NewRank = 0;
  E0 = Evaluation(Sol); /* Evaluation Method */
  logDebug(cout << tag << "bNow = " << get_elements() << endl);
  if (get_elements() == 0) {
    NewRank = 1;
    logDebug(cout << tag << "Calculate Hash (bNow = 0)" << endl);
    Hash0 = Sol.Hash();
    loc0 = Add (Sol);
    LastChange[loc0] = NowTime;
    logDebug(cout << tag << "Finish" << endl);
    log_depth--;
    return true;
  }
  else {
    if (E0 >= worst() && it_is_full ()) {
      logDebug(cout << tag << "Finish" << endl);
      log_depth--;
      return false;
    }
    logDebug(cout << tag << "Calculate Hash (bNow > 0)" << endl);
    Hash0 = Sol.Hash();

    if (E0 < best()) {
      NewRank = 1;
    }
    else {
      for (int i = get_elements()-1;i >= 0; i--) {
	if (EqualSol(Sol,i)) {
	  logDebug(cout << tag << "Finish" << endl);
	  log_depth--;
	  return false;
	}
	else if (evaluation(i) < E0){
	  NewRank = i+2;
	  loc0 = Add (Sol);
	  LastChange[loc0] = NowTime;
	  logDebug(cout << tag << "Finish" << endl);
	  log_depth--;
	  return true;
	}
      }
      NewRank = 1; /* This case is when the solution gives the same evaluation 
		      as the best solution */
    }
    loc0 = Add (Sol);
    LastChange[loc0] = NowTime;
  }
  logDebug(cout << tag << "Finish" << endl);
  log_depth--;
  return true;
}

void RefSet::Call_Improvement (SQM_solution &Sol) {
  Improvement(Sol);
}

void RefSet::SubsetControl () {
  /* Initialization */
  int iLoc;
  logDebug(cout << "RefSet::SubsetControl: Start" << endl);
  for (iLoc = 0;iLoc < get_elements();iLoc++)
    LastChange[iLoc] = 0;
  SubsetType = invalid_subset;
  for (int i = 0;i < invalid_subset;i++)
    LastRunTime[i] = 0;
  NowTime = 0;
  StopCondition = 0;
  SubsetType = two_element;
  LocNew = new int [get_elements()];
  LocOld = new int [get_elements()];
  while (StopCondition == 0) {

    /*++SubsetType;*/
    NowTime++;

    iNew = 0;
    jOld = 0;
    for (int i = 0;i < get_elements();i++) {
      iLoc = location(i);
      if (LastChange[iLoc] >= LastRunTime[SubsetType])
	LocNew[iNew++] = iLoc;
      else
	LocOld[jOld++] = iLoc;
    }

    if (iNew == 0) {
      delete [] LocNew;
      delete [] LocOld;
      return;
    }
    logDebug(cout << "Continue SubsetControl with "
	     << 100 * iNew / get_elements()
	     << "% of new solutions" << endl);

    switch (SubsetType) {
    case two_element: /* Call Algorithm 1 Subrutine */
      algorithm_for_SubsetType1 ();
      break;
    case three_element: /* Call Algorithm 2 Subrutine */
      algorithm_for_SubsetType2 ();
      break;
    case four_element: /* Call Algorithm 3 Subrutine */
      algorithm_for_SubsetType3 ();
      break;
    case best_i: /* Call Algorithm 4 Subrutine */
      algorithm_for_SubsetType4 ();
      break;
    }

    if (StopCondition > 0) {
      /* Actually no StopCondition while new solutions were found */
    }
    LastRunTime[SubsetType] = NowTime;
  }
}

void RefSet::algorithm_for_SubsetType1 () {
  int iLoc,jLoc;
  SQM_solution *X,*Y;
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
	    /* SQM_path_relinking(*this,X,Y); */
	    /* Optional Check: if LastChange[iLoc] == NowTime, 
	       then can jump to the end of "i loop" to pick up the next i,
	       and generate fewer solutions. */
	    if (LastChange[iLoc] == NowTime)
	      break;
	  }
	}
	clean_garbage();
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
	    /* SQM_path_relinking(*this,X); */
	    /* Optional Check: if LastChange[iLoc] == NowTime, 
	       then can jump to the end of "i loop" to pick up the next i,
	       and generate fewer solutions. */
	    if (LastChange[iLoc] == NowTime)
	      break;
	  }
	}
	clean_garbage();
      }
    }
}


void RefSet::algorithm_for_SubsetType2 () {
  int loc1,iLoc,jLoc;
  list<SQM_solution*> X;
  logDebug(cout << "Start algorithm for SubsetType2" << endl);
  loc1 = location(0);
  X.push_back(Solutions[loc1]);
  if (LastChange[loc1] >= LastRunTime[SubsetType]) {
    for (int i = 1;i < get_elements();i++) {
      iLoc = location(i);
      if (LastChange[iLoc] < NowTime) {
	X.push_back(Solutions[iLoc]);
	for (int j = i+1;j < get_elements();j++) {
	  jLoc = location(j);
	  if (LastChange[jLoc] < NowTime) {
	    X.push_back(Solutions[jLoc]);
	    /* Create C(X) and execute improvement method */
	    SQM_path_relinking(*this,X);
	    clean_garbage();
	    X.pop_back();
	    /* Optional Check: if LastChange[iLoc] == NowTime, 
	       then can jump to the end of "i loop" to pick up the next i,
	       and generate fewer solutions. */
	  }
	}
	X.pop_back();
      }
    }
  }
  else {
    if (iNew > 1) 
      for (int i = 0;i < iNew;i++) {
	iLoc = LocNew[i];
	if (LastChange[iLoc] < NowTime) {
	  X.push_back(Solutions[iLoc]);
	  for (int j = i+1;j < iNew;j++) {
	    jLoc = LocNew[j];
	    if (LastChange[jLoc] < NowTime) {
	      X.push_back(Solutions[jLoc]);
	      /* Create C(X) and execute improvement method */
	      SQM_path_relinking(*this,X);
	      clean_garbage();
	      X.pop_back();
	      /* Optional Check: if LastChange[iLoc] == NowTime, 
		 then can jump to the end of "i loop" to pick up the next i,
		 and generate fewer solutions. */
	    }
	  }
	  X.pop_back();
	}
      }
    if (jOld > 0)
      for (int i = 0;i < iNew;i++) {
	iLoc = LocNew[i];
	if (LastChange[iLoc] < NowTime) {
	  X.push_back(Solutions[iLoc]);
	  for (int j = 0;j < jOld;j++) {
	    jLoc = LocOld[j];
	    if (LastChange[jLoc] < NowTime) {
	      X.push_back(Solutions[jLoc]);
	      /* Create C(X) and execute improvement method */
	      SQM_path_relinking(*this,X);
	      clean_garbage();
	      X.pop_back();
	      /* Optional Check: if LastChange[iLoc] == NowTime, 
		 then can jump to the end of "I loop" to pick up the next I,
		 and generate fewer solutions. */
	    }
	  }
	  X.pop_back();
	}
      }
  }
}

void RefSet::algorithm_for_SubsetType3 () {
  int loc1,loc2,iLoc,jLoc;
  list<SQM_solution*> X;
  logDebug(cout << "Start algorithm for SubsetType3" << endl);
  loc1 = location(0);
  loc2 = location(1);
  X.push_back(Solutions[loc1]);
  X.push_back(Solutions[loc2]);
  if (LastChange[loc1] >= LastRunTime[SubsetType] || 
      LastChange[loc2] >= LastRunTime[SubsetType]) {
    for (int i = 2;i < get_elements();i++) {
      iLoc = location(i);
      if (LastChange[iLoc] < NowTime) {
	X.push_back(Solutions[iLoc]);
	for (int j = i+1;j < get_elements();j++) {
	  jLoc = location(j);
	  if (LastChange[jLoc] < NowTime) {
	    X.push_back(Solutions[jLoc]);
	    /* Create C(X) and execute improvement method */
	    SQM_path_relinking(*this,X);
	    clean_garbage();
	    X.pop_back();
	  }
	}
	X.pop_back();
      }
    }
  }
  else {
    if (iNew > 1) 
      for (int i = 0;i < iNew;i++) {
	iLoc = LocNew[i];
	if (LastChange[iLoc] < NowTime) {
	  X.push_back(Solutions[iLoc]);
	  for (int j = i+1;j < iNew;j++) {
	    jLoc = LocNew[j];
	    if (LastChange[jLoc] < NowTime) {
	      X.push_back(Solutions[jLoc]);
	      /* Create C(X) and execute improvement method */
	      SQM_path_relinking(*this,X);
	      clean_garbage();
	      X.pop_back();
	    }
	  }
	  X.pop_back();
	}
      }
    if (jOld > 0)
      for (int i = 0;i < iNew;i++) {
	iLoc = LocNew[i];
	if (LastChange[iLoc] < NowTime) {
	  X.push_back(Solutions[iLoc]);
	  for (int j = 0;j < jOld;j++) {
	    jLoc = LocOld[j];
	    if (LastChange[jLoc] < NowTime) {
	      X.push_back(Solutions[jLoc]);
	      /* Create C(X) and execute improvement method */
	      SQM_path_relinking(*this,X);
	      clean_garbage();
	      X.pop_back();
	    }
	  }
	  X.pop_back();
	}
      }
  }
}

void RefSet::algorithm_for_SubsetType4 () {
  int iLoc;
  bool new_subset = false;
  list<SQM_solution*> X;
  logDebug(cout << "Start algorithm for SubsetType4" << endl);
  for (int i = 0;i < 4;i++) {
    iLoc = location(i);
    X.push_back(Solutions[iLoc]);
    if (LastChange[iLoc] >= LastRunTime[SubsetType]) new_subset = true;
  }
  for (int i = 4;i < get_elements();i++) {
    iLoc = location(i);
    if (LastChange[iLoc] >= LastRunTime[SubsetType]) new_subset = true;
    if (LastChange[iLoc] < NowTime) {
      X.push_back(Solutions[iLoc]);
      if (new_subset) {
	/* Create C(X) and execute improvement method */
	SQM_path_relinking(*this,X);
	clean_garbage();
      }
    }
  }
}

Subset& operator++(Subset &target) {
  target = (target == invalid_subset ?
	    static_cast<Subset>(0) : static_cast<Subset>(target + 1));
  if (target == invalid_subset)
    target = static_cast<Subset>(0);
  return target;
}

double get_response_time (SQM_solution &Sol) {
  return Sol.get_response_time();
}

double get_perfect_matching_cost (SQM_solution &Sol) {
  return Sol.pm_cost;
}
