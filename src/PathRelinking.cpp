
#include <iostream>
#include "PathRelinking.h"
#include "PerfectMatching.h"
#include "log.h"

using namespace std;

bool incompatible_solutions(SQM_solution*,SQM_solution*);
int* PR_run_perfect_matching(SQM_solution*,SQM_solution*);
int* PR_determine_order_(SQM_solution*,int*,SQM_solution*);
int* PR_determine_order_nf(SQM_solution*,int*,SQM_solution*);
int* PR_determine_order_ff(SQM_solution*,int*,SQM_solution*);

list<SQM_solution*>* Path_Relinking (SQM_solution *X,SQM_solution *Y) {
  logDebug(cout << "START:\t***\tPath_Relinking\t***" << endl);
  if (incompatible_solutions(X,Y)) {
    logError(cerr << "Call Path_Relinking with incompatible solutions" << endl);
    return NULL;
  }

  /* Run perfect matching */
  logDebug(cout << "Run perfect matching" << endl);
  int *pm;
  pm = PR_run_perfect_matching(X,Y);

  /* Determine order of change */
  logDebug(cout << "Determine order of change" << endl);
  int p = X->get_servers();
  int *order = PR_determine_order_nf(X,pm,Y);

  logDebug(cout << "Create solutions" << endl);
  list<SQM_solution*> *Solutions;
  int x,y,loc_x,loc_y;
  SQM_solution *Z;

  Solutions = new list<SQM_solution*>;
  logDebug(cout << "Create list to store solutions" << endl);
  for (int step = 0;step < p - 1;step++) {
    logDebug(cout << "Step " << step+1 << endl);
    x = order[step];
    y = pm[x];
    logDebug(cout << "process server " << x << ", change with srever " << y << endl);
    loc_x = X->get_server_location(x);
    loc_y = Y->get_server_location(y);
    logDebug(cout << "\twith locations " << loc_x << " and " << loc_y << endl);
    if (loc_x != loc_y) {
      Z = X->clone();
      Z->set_server_location(x,loc_y);
      Solutions->push_back(Z);
      X = Z;
    }
  }
  
  delete [] order;
  delete [] pm;
  logDebug(cout << "FINISH:\t***\tPath_Relinking\t***" << endl);
  return Solutions;
}

bool incompatible_solutions(SQM_solution *X,SQM_solution *Y) {
  if (LogDebug) {
    if (X->get_servers() != Y->get_servers())
      cout << "/* No Same number of servers */" << endl;
    if (X->get_instance() != Y->get_instance())
      cout << "/* No Same instance */" << endl;
  }
  return ((X->get_instance() != Y->get_instance()) || /* Same instance */
	  (X->get_servers() != Y->get_servers()) /* Same number of servers */
	  );
}

/* Create distances matrix and call Perfect_Matching procedure */
int* PR_run_perfect_matching(SQM_solution *X,SQM_solution *Y) {
  int p = X->get_servers();
  int q = Y->get_servers();
  SQM_instance *I = X->get_instance();
  point *site = I->site(0);
  double **distances;
  int loc_x,loc_y;
  int *pm;
  distances = new double*[p];
  for (int i = 0;i < p;i++) distances[i] = new double[q];
  for (int i = 0;i < p;i++) {
    loc_x = X->get_server_location(i);
    for (int l = 0;l < q;l++) {
      loc_y = Y->get_server_location(l);
      distances[i][l] = dist(&(site[loc_x]),&(site[loc_y]));
    }
  }
  pm = Perfect_Matching(p,distances);
  for (int i = 0;i < p;i++) delete [] distances[i];
  delete [] distances;
  return pm;
}

int* PR_determine_order_(SQM_solution *X,int *pm,SQM_solution *Y) {
  int p = X->get_servers();
  int *order = new int [p];
  for (int i = 0;i < p;i++) order[i] = i; 
  return order;
}

/* nearest first */
int* PR_determine_order_nf(SQM_solution *X,int *pm,SQM_solution *Y) {
  int *order,loc_x,loc_y;
  int p = X->get_servers();
  double *d;
  SQM_instance *I = X->get_instance();
  point *site = I->site(0);

  order = new int [p];
  d = new double [p];
  for (int i = 0;i < p;i++) {
    loc_x = X->get_server_location(i);
    loc_y = Y->get_server_location(pm[i]);
    d[i] = dist(&(site[loc_x]),&(site[loc_y]));
  }
  sort_dist(p,d,order);
  return order;
}

/* farthest first */
int* PR_determine_order_ff(SQM_solution *X,int *pm,SQM_solution *Y) {
  int *order,loc_x,loc_y;
  int p = X->get_servers();
  double *d;
  SQM_instance *I = X->get_instance();
  point *site = I->site(0);

  order = new int [p];
  d = new double [p];
  for (int i = 0;i < p;i++) {
    loc_x = X->get_server_location(i);
    loc_y = Y->get_server_location(pm[i]);
    d[i] = -dist(&(site[loc_x]),&(site[loc_y]));
  }
  sort_dist(p,d,order);
  return order;
}
