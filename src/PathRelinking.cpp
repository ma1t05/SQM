
#include "PathRelinking.h"

bool compatible_solutions(SQM_Solution*,SQM_Solution*);

list<SQM_solution*>* Path_Relinking (SQM_Solution *X,SQM_Solution *Y) {
  int p;
  SQM_solution *Z,;
  list<SQM_solution*> *Solutions;

  if (!compatible_solutions(X,Y)) {
    logError(cerr << "Call Path_Relinking with incompatible solutions" << endl);
    return NULL;
  }

  /* Run perfect matching */
  int *pm;
  double **distances;
  distances = new double*[p];
  for (int i = 0;i < p;i++) distances[i] = new double[p];
  pm = Perfect_Matching(p,distances);
  
  /* Determine order of change */
  int *order;
  order = new int [p];
  /* pening: determine a better way to approach X to Y */
  for (int i = 0;i < p;i++) order[i] = i; 

  Solutions = new list<SQM_solution*>;
  for (int step = 0;step < p - 1;step++) {
    x = order[step];
    y = pm[x];
    location_x = X->get_server_location(i);
    location_y = Y->get_server_location(pm[i]);
    if (location_x != location_y) {
      Z = X->clone();
      Z->set_server_location(x,location_y);
      Solutions->push_back(Z);
      X = Z;
    }
  }
  
  delete [] order;
  delete [] pm;
  return Solutions;
}

bool compatible_solutions(SQM_Solution*,SQM_Solution*) {
  return ((X->get_instance() != Y->get_instance()) || /* Same instance */
	  (X->get_servers() != Y->get_server()) /* Same number of servers */
	  );
}
