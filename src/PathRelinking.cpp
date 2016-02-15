
#include "PathRelinking.h"

using namespace std;

bool incompatible_solutions(SQM_solution&,SQM_solution&);
double** PR_distances_matrix(SQM_solution&,SQM_solution&);

list<SQM_solution*>* Path_Relinking (SQM_solution &X,SQM_solution &Y) {
  log_depth++;
  string tag = log_tag("Path_Relinking: ");
  logDebug(cout << tag << "Start" << endl);
  if (incompatible_solutions(X,Y)) {
    logError(cerr << tag << "Call with incompatible solutions" << endl);
    log_depth--;
    return NULL;
  }

  /* Run perfect matching */
  logDebug(cout << tag << "Run matching" << endl);
  int *pm;
  pm = matching_function(X,Y);

  /* Determine order of change */
  logDebug(cout << tag << "Determine order of change" << endl);
  int p = X.get_servers();
  int *order = order_function(X,pm,Y);

  logDebug(cout << tag << "Create solutions" << endl);
  list<SQM_solution*> *Solutions;
  int x,y,loc_x,loc_y;
  SQM_solution *Z;

  Solutions = new list<SQM_solution*>;
  logDebug(cout << tag << "Create list to store solutions" << endl);
  Z = &X;
  for (int step = 0;step < p - 1;step++) {
    x = order[step];
    y = pm[x];
    loc_x = Z->get_server_location(x);
    loc_y = Y.get_server_location(y);
    if (loc_x != loc_y) {
      logDebug
	(cout << tag << "Step " << step+1 << " change server " << x 
	 << ", with srever " << y << " with locations " << loc_x << " and " 
	 << loc_y << endl
	 );
      Z = Z->clone();
      Z->set_server_location(x,loc_y);
      Solutions->push_back(Z);
    }
  }
  
  delete [] order;
  delete [] pm;
  logDebug(cout << tag << "Finish" << endl);
  log_depth--;
  return Solutions;
}

bool incompatible_solutions(SQM_solution &X,SQM_solution &Y) {
  if (LogInfo) {
    if (X.get_servers() != Y.get_servers())
      logError(cerr << "Error:No Same number of servers" << endl);
    if (X.get_instance() != Y.get_instance())
      logError(cerr << "Error:No Same instance" << endl);
  }
  return ((X.get_instance() != Y.get_instance()) || /* Same instance */
	  (X.get_servers() != Y.get_servers()) /* Same number of servers */
	  );
}

matching_type& operator++(matching_type &target) {
  target = static_cast<matching_type>(target + 1);
  return target;
}

/* Create distances matrix and call Perfect_Matching procedure */
int* PR_perfect_matching(SQM_solution &X,SQM_solution &Y) {
  int p = X.get_servers();
  int *pm;
  double **distances;
  distances = PR_distances_matrix(X,Y);
  pm = Perfect_Matching(p,distances);
  for (int i = 0;i < p;i++) delete [] distances[i];
  delete [] distances;
  return pm;
}

int* PR_workload_matching(SQM_solution &X,SQM_solution &Y) {
  int p = X.get_servers();
  int *pm;
  int *a,*b;
  double *wl_x,*wl_y;

  pm = new int [p];
  wl_x = X.get_workload();
  wl_y = Y.get_workload();
  a = new int [p];
  b = new int [p];
  sort_dist(p,wl_x,a);
  sort_dist(p,wl_y,b);
  for (int i = 0;i < p;i++) 
    pm[a[i]] = b[i];
  delete [] b;
  delete [] a;
  delete [] wl_y;
  delete [] wl_x;
  return pm;
}

int* PR_random_matching(SQM_solution &X,SQM_solution &Y) {
  int p = X.get_servers();
  int *pm,pos;
  bool *used;

  pm = new int [p];
  used = new bool [p];
  for (int i = 0;i < p;i++) used[i] = false;
  for (int i = 0;i < p;i++) {
    pm[i] = -1;
    do {
      pos = unif(p);
      if (used[pos] == false) {
	pm[i] = pos;
	used[pos] = true;
      }
    } while (pm[i] == -1);
  }
  delete [] used;
  return pm;
}

order_type& operator++(order_type &target) {
  target = static_cast<order_type>(target + 1);
  return target;
}

double PR_perfect_matching_cost(SQM_solution &X,SQM_solution &Y) {
  int p = X.get_servers();
  double cost;
  double **distances;
  distances = PR_distances_matrix(X,Y);
  cost = Perfect_Matching_cost(p,distances);
  for (int i = 0;i < p;i++) delete [] distances[i];
  delete [] distances;
  return cost;
}

/* Create distances matrix  */
double** PR_distances_matrix(SQM_solution &X,SQM_solution &Y) {
  int p = X.get_servers();
  int q = Y.get_servers();
  SQM_instance *I = X.get_instance();
  point *site = I->site(0);
  double **distances;
  int loc_x,loc_y;
  int *pm;
  distances = new double*[p];
  for (int i = 0;i < p;i++) distances[i] = new double[q];
  for (int i = 0;i < p;i++) {
    loc_x = X.get_server_location(i);
    for (int l = 0;l < q;l++) {
      loc_y = Y.get_server_location(l);
      distances[i][l] = dist(&(site[loc_x]),&(site[loc_y]));
    }
  }
  return distances;
}

int* PR_processing_order_random(SQM_solution &X,int *pm,SQM_solution &Y) {
  int p = X.get_servers(),pos;
  int *order = new int [p];
  bool *used = new bool [p];
  for (int i = 0;i < p;i++) used[i] = false;
  for (int i = 0;i < p;i++) {
    order[i] = -1;
    do {
      pos = unif(p);
      if (used[pos] == false) {
	order[i] = pos;
	used[pos] = true;
      }
    } while (order[i] == -1);
  }
  for (int i = 0;i < p;i++) order[i] = i; 
  delete [] used;
  return order;
}

/* nearest first */
int* PR_processing_order_nf(SQM_solution &X,int *pm,SQM_solution &Y) {
  int *order,loc_x,loc_y;
  int p = X.get_servers();
  double *d;
  SQM_instance *I = X.get_instance();
  point *site = I->site(0);

  order = new int [p];
  d = new double [p];
  for (int i = 0;i < p;i++) {
    loc_x = X.get_server_location(i);
    loc_y = Y.get_server_location(pm[i]);
    d[i] = dist(&(site[loc_x]),&(site[loc_y]));
  }
  sort_dist(p,d,order);
  return order;
}

/* farthest first */
int* PR_processing_order_ff(SQM_solution &X,int *pm,SQM_solution &Y) {
  int *order,loc_x,loc_y;
  int p = X.get_servers();
  double *d;
  SQM_instance *I = X.get_instance();
  point *site = I->site(0);

  order = new int [p];
  d = new double [p];
  for (int i = 0;i < p;i++) {
    loc_x = X.get_server_location(i);
    loc_y = Y.get_server_location(pm[i]);
    d[i] = -dist(&(site[loc_x]),&(site[loc_y]));
  }
  sort_dist(p,d,order);
  return order;
}

void SQM_path_relinking(RefSet &EliteSols,list<SQM_solution*> &Solutions) {
  int total_improved_solutions = 0;
  SQM_solution *X,*Y;
  list<SQM_solution*>::iterator it1,it2,end,Z;
  list<SQM_solution*> *pr_sols;
  SQM_solution *Best,*Best_input;
  int N;
  string tag = "SQM_path_relinking: ";

  logDebug(cout << tag << "Start" << endl);
  it1 = Solutions.begin();
  end = Solutions.end();
  while (it1 != end) {
    X = *it1;
    it2 = ++it1;
    while (it2 != end) {
      Y = *it2; it2++;
      pr_sols = Path_Relinking(X,Y);
      if (pr_sols != NULL) {
	for (Z = pr_sols->begin();Z != pr_sols->end();Z++) {
	  Improvement_Method(**Z); /* Improvement method */
	  if (EliteSols.Update(**Z))
	    total_improved_solutions++;
	  else delete *Z;
	}
	delete pr_sols;
      }
    }
  }
  logDebug(cout << tag << "Finish" << endl);
}

void Path_Relinking(RefSet &EliteSols,SQM_solution &X,SQM_solution &Y) {
  list<SQM_solution*>::iterator Z;
  list<SQM_solution*> *pr_sols;
  SQM_solution *Best,*Best_input;
  int N;
  string tag = "RefSet::Path_Relinking: ";

  logDebug(cout << tag << "Start" << endl);
  pr_sols = Path_Relinking(X,Y);
  if (pr_sols != NULL) {
    for (Z = pr_sols->begin();Z != pr_sols->end();Z++) {
      Improvement_Method(**Z); /* Improvement method */
      if (!EliteSols.Update(**Z))
	delete *Z;
    }
    delete pr_sols;
  }
  logDebug(cout << tag << "Finish" << endl);
}

SQM_solution* SQM_best_solution(list<SQM_solution*>* Solutions) {
  SQM_solution *Best;
  list<SQM_solution*>::iterator X;

  Best = NULL;
  for (X = Solutions->begin();X != Solutions->end();X++) 
    if (Best == NULL || *Best > **X)
      Best = *X;

  return Best;
}

SQM_solution* SQM_leave_only_the_best(list<SQM_solution*>* Solutions) {
  SQM_solution *Best;
  list<SQM_solution*>::iterator X;

  Best = NULL;
  for (X = Solutions->begin();X != Solutions->end();X++) 
    if (Best == NULL || *Best > **X) {
      if (Best != NULL) delete Best;
      Best = *X;
    }
    else delete *X;
  delete Solutions;
  return Best;
}

void SQM_delete_sols(list<SQM_solution*>* Solutions) {
  list<SQM_solution*>::iterator X;
  for (X = Solutions->begin();X != Solutions->end();X++) 
    delete *X;
  delete Solutions;
}

double SQM_min_cost_pm(RefSet &Sols,SQM_solution &Sol) {
  double cost,min_cost;
  int bNow = Sols.get_elements();

  min_cost = PR_perfect_matching_cost(*Sols.get_sol(0),Sol);
  for (int i = 1;i < bNow;i++) {
    cost = PR_perfect_matching_cost(*Sols.get_sol(i),Sol);
    if (min_cost > cost) min_cost = cost;
  }

  return min_cost;
}
