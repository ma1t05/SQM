
#include "Reference_Set.h"
#include "PathRelinking.h"

RefSet::RefSet (int Max) {
  bMax = Max;
  bNow = 0;
  RefSetCall = 0;
  RefSetAdd = 0;
  DupCheck = 0;
  FullDupCheck = 0;
  FullDupFound = 0;
  loc = new int [bMax];
  E = new double [bMax];
  Hash = new double [bMax];
  Solutions = new SQM_solution*[bMax];
  LastChange = new int [bMax];
  LastRunTime = new int [invalid_subset];
}

RefSet::~RefSet () {
  for (int i = 0;i < bNow;i++) {
    delete Solutions[loc[i]];
  }
  delete [] LastRunTime;
  delete [] LastChange;
  delete [] Solutions;
  delete [] Hash;
  delete [] E;
  delete [] loc;
}

void RefSet::Add (SQM_solution &Sol) {
  log_depth++;
  string tag = log_tag("RefSet::Add: ");
 
  logDebug(cout << tag << "Start" << endl);
  int loc0;
  RefSetAdd++;
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
  LastChange[loc0] = NowTime;
  logDebug(cout << tag << "Finish" << endl);
  log_depth--;
}

bool RefSet::Update (SQM_solution &Sol) {
  log_depth++;
  string tag = log_tag("RefSet::Update: ");
  logDebug(cout << tag << "Start" << endl);
  RefSetCall++;
  NewRank = 0;
  E0 = Sol.get_response_time();
  logDebug(cout << tag << "bNow = " << bNow << endl);
  if (bNow == 0) {
    NewRank = 1;
    logDebug(cout << tag << "Calculate Hash (bNow = 0)" << endl);
    Hash0 = Sol.Hash();
    Add (Sol);
    logDebug(cout << tag << "Finish" << endl);
    log_depth--;
    return true;
  }
  else {
    if (E0 >= worst() && bNow == bMax) {
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
      for (int i = bNow-1;i >= 0; i--) {
	if (E[loc[i]] == E0) {
	  DupCheck++;
	  if (Hash[loc[i]] == Hash0) {
	    FullDupCheck++;
	    if (Sol == *Solutions[loc[i]]) {
	      FullDupFound++;
	      logDebug(cout << tag << "Finish" << endl);
	      log_depth--;
	      return false;
	    }
	  }
	}
	else if (E[loc[i]] < E0){
	  NewRank = i+2;
	  Add (Sol);
	  logDebug(cout << tag << "Finish" << endl);
	  log_depth--;
	  return true;
	}
      }
      NewRank = 1; /* This case is when the solution gives the same evaluation 
		      as the best solution */
    }
    Add (Sol);
  }
  logDebug(cout << tag << "Finish" << endl);
  log_depth--;
  return true;
}

void RefSet::SubsetControl () {
  /* Initialization */
  int iLoc;
  logDebug(cout << "RefSet::SubsetControl: Start" << endl);
  for (iLoc = 0;iLoc < bMax;iLoc++)
    LastChange[iLoc] = 0;
  SubsetType = invalid_subset;
  for (int i = 0;i < invalid_subset;i++)
    LastRunTime[i] = 0;
  NowTime = 0;
  StopCondition = 0;
  SubsetType = two_element;
  LocNew = new int [bMax];
  LocOld = new int [bMax];
  while (StopCondition == 0) {

    /*++SubsetType;*/
    NowTime++;

    iNew = 0;
    jOld = 0;
    for (int i = 0;i < bNow;i++) {
      iLoc = loc[i];
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
    logDebug(cout << "Continue SubsetControl with " << 100 * iNew / bNow
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

void RefSet::clean_garbage () {
  logDebug(cout << "Clean garbage: " << garbage.size() << endl);
  while (!garbage.empty()) {
    delete garbage.back();
    garbage.pop_back();
  }
}

void RefSet::algorithm_for_SubsetType1 () {
  int iLoc,jLoc;
  list<SQM_solution*> X;
  logDebug(cout << "Start algorithm for SubsetType1" << endl);
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
	    X.pop_back();
	    /* Optional Check: if LastChange[iLoc] == NowTime, 
	       then can jump to the end of "i loop" to pick up the next i,
	       and generate fewer solutions. */
	    if (LastChange[iLoc] == NowTime)
	      break;
	  }
	}
	X.pop_back();
	clean_garbage();
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
	    X.pop_back();
	    /* Optional Check: if LastChange[iLoc] == NowTime, 
	       then can jump to the end of "i loop" to pick up the next i,
	       and generate fewer solutions. */
	    if (LastChange[iLoc] == NowTime)
	      break;
	  }
	}
	X.pop_back();
	clean_garbage();
      }
    }
}


void RefSet::algorithm_for_SubsetType2 () {
  int loc1,iLoc,jLoc;
  list<SQM_solution*> X;
  logDebug(cout << "Start algorithm for SubsetType2" << endl);
  loc1 = loc[0];
  X.push_back(Solutions[loc1]);
  if (LastChange[loc1] >= LastRunTime[SubsetType]) {
    for (int i = 1;i < bNow;i++) {
      iLoc = loc[i];
      if (LastChange[iLoc] < NowTime) {
	X.push_back(Solutions[iLoc]);
	for (int j = i+1;j < bNow;j++) {
	  jLoc = loc[j];
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
  loc1 = loc[0];
  loc2 = loc[1];
  X.push_back(Solutions[loc1]);
  X.push_back(Solutions[loc2]);
  if (LastChange[loc1] >= LastRunTime[SubsetType] || 
      LastChange[loc2] >= LastRunTime[SubsetType]) {
    for (int i = 2;i < bNow;i++) {
      iLoc = loc[i];
      if (LastChange[iLoc] < NowTime) {
	X.push_back(Solutions[iLoc]);
	for (int j = i+1;j < bNow;j++) {
	  jLoc = loc[j];
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
    iLoc = loc[i];
    X.push_back(Solutions[iLoc]);
    if (LastChange[iLoc] >= LastRunTime[SubsetType]) new_subset = true;
  }
  for (int i = 4;i < bNow;i++) {
    iLoc = loc[i];
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
