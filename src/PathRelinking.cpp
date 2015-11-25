
#include <iostream>
#include "PathRelinking.h"
#include "PerfectMatching.h"
#include "log.h"

using namespace std;

bool compatible_solutions(SQM_solution*,SQM_solution*);
int* PR_run_perfect_matching(SQM_solution*,SQM_solution*);

list<SQM_solution*>* Path_Relinking (SQM_solution *X,SQM_solution *Y) {
  int p;
  SQM_solution *Z;
  list<SQM_solution*> *Solutions;
  int *pm,x,y;
  int loc_x,loc_y;

  if (!compatible_solutions(X,Y)) {
    logError(cerr << "Call Path_Relinking with incompatible solutions" << endl);
    return NULL;
  }

  /* Run perfect matching */
  pm = PR_run_perfect_matching(X,Y);

  /* Determine order of change */
  int *order;
  order = new int [p];
  /* pending: determine a better way to approach X to Y */
  for (int i = 0;i < p;i++) order[i] = i; 

  Solutions = new list<SQM_solution*>;
  for (int step = 0;step < p - 1;step++) {
    x = order[step];
    y = pm[x];
    loc_x = X->get_server_location(x);
    loc_y = Y->get_server_location(y);
    if (loc_x != loc_y) {
      Z = X->clone();
      Z->set_server_location(x,loc_y);
      Solutions->push_back(Z);
      X = Z;
    }
  }
  
  delete [] order;
  delete [] pm;
  return Solutions;
}

bool compatible_solutions(SQM_solution *X,SQM_solution *Y) {
  if (X->get_servers() != Y->get_servers())
    cout << "/* No Same number of servers */" << endl;

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
