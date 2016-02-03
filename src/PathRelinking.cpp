
#include "PathRelinking.h"

using namespace std;

bool incompatible_solutions(SQM_solution*,SQM_solution*);
double** PR_distances_matrix(SQM_solution*,SQM_solution*);

list<SQM_solution*>* Path_Relinking (SQM_solution *X,SQM_solution *Y) {
  log_depth++;
  string tag = log_tag("Path_Relinking: ");
  logDebug(cout << tag << "Start" << endl);
  if (incompatible_solutions(X,Y)) {
    logError(cerr << tag << "Call with incompatible solutions" << endl);
    return NULL;
  }

  /* Run perfect matching */
  logDebug(cout << tag << "Run matching" << endl);
  int *pm;
  pm = matching_function(X,Y);

  /* Determine order of change */
  logDebug(cout << tag << "Determine order of change" << endl);
  int p = X->get_servers();
  int *order = order_function(X,pm,Y);

  logDebug(cout << tag << "Create solutions" << endl);
  list<SQM_solution*> *Solutions;
  int x,y,loc_x,loc_y;
  SQM_solution *Z;

  Solutions = new list<SQM_solution*>;
  logDebug(cout << tag << "Create list to store solutions" << endl);
  for (int step = 0;step < p - 1;step++) {
    x = order[step];
    y = pm[x];
    loc_x = X->get_server_location(x);
    loc_y = Y->get_server_location(y);
    if (loc_x != loc_y) {
      logDebug
	(cout << tag << "Step " << step+1 << " change server " << x 
	 << ", with srever " << y << " with locations " << loc_x << " and " 
	 << loc_y << endl
	 );
      Z = X->clone();
      Z->set_server_location(x,loc_y);
      Solutions->push_back(Z);
      X = Z;
    }
  }
  
  delete [] order;
  delete [] pm;
  logDebug(cout << tag << "Finish" << endl);
  log_depth--;
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

matching_type& operator++(matching_type &target) {
  target = static_cast<matching_type>(target + 1);
  return target;
}

/* Create distances matrix and call Perfect_Matching procedure */
int* PR_run_perfect_matching(SQM_solution *X,SQM_solution *Y) {
  int p = X->get_servers();
  int *pm;
  double **distances;
  distances = PR_distances_matrix(X,Y);
  pm = Perfect_Matching(p,distances);
  for (int i = 0;i < p;i++) delete [] distances[i];
  delete [] distances;
  return pm;
}

int* PR_workload_matching(SQM_solution *X,SQM_solution *Y) {
  int p = X->get_servers();
  int *pm;
  int *a,*b;
  double *wl_x,*wl_y;

  pm = new int [p];
  wl_x = X->get_workload();
  wl_y = Y->get_workload();
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

int* PR_random_matching(SQM_solution *X,SQM_solution *Y) {
  int p = X->get_servers();
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

double PR_perfect_matching_cost(SQM_solution *X,SQM_solution *Y) {
  int p = X->get_servers();
  double cost;
  double **distances;
  distances = PR_distances_matrix(X,Y);
  cost = Perfect_Matching_cost(p,distances);
  for (int i = 0;i < p;i++) delete [] distances[i];
  delete [] distances;
  return cost;
}

/* Create distances matrix  */
double** PR_distances_matrix(SQM_solution *X,SQM_solution *Y) {
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
  return distances;
}

int* PR_processing_order_random(SQM_solution *X,int *pm,SQM_solution *Y) {
  int p = X->get_servers(),pos;
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
int* PR_processing_order_nf(SQM_solution *X,int *pm,SQM_solution *Y) {
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
int* PR_processing_order_ff(SQM_solution *X,int *pm,SQM_solution *Y) {
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

SQM_solution* SQM_path_relinking(list<SQM_solution*>* Solutions) {
  int total_improved_solutions = 0;
  double Tr_x,Tr_y,Best_TR,TR;
  list<SQM_solution*>::iterator X,Y,Z;
  RefSet *EliteSols;
  list<SQM_solution*> *path_relinking_sols,*improved_solutions;
  SQM_solution *Best,*Best_input;
  clock_t beginning,now;
  double best_rt,avg_rt,worst_rt;
  int N;
  string tag = "SQM_path_relinking: ";

  EliteSols = new RefSet(Solutions->size());
  for (X = Solutions->begin();X != Solutions->end();X++)
    EliteSols->Update(**X);

  logLevel = LOG_DEBUG;
  logDebug(cout << tag << "Start" << endl);
  beginning = clock();
  improved_solutions = new list<SQM_solution*>;
  best_rt = 100.0,worst_rt = 0.0,avg_rt = 0.0, N = 0;
  for (X = Solutions->begin();X != Solutions->end();X++) {
    Tr_x = (*X)->get_response_time(); 
    Y = X;
    Y++;
    for (;Y != Solutions->end();Y++) {
      Tr_y = (*Y)->get_response_time();
      Best_TR = (Tr_x > Tr_y ? Tr_y : Tr_x);
      path_relinking_sols = Path_Relinking(*X,*Y);
      if (path_relinking_sols != NULL) {
	for (Z = path_relinking_sols->begin();Z != path_relinking_sols->end();Z++) {
	  SQM_heuristic(*Z);
	  TR = (*Z)->get_response_time();
	  avg_rt += TR; N++;
	  if (best_rt > TR) best_rt = TR;
	  if (worst_rt < TR) worst_rt = TR;
	  if (TR < Best_TR) {
	    total_improved_solutions++;
	    improved_solutions->push_back(*Z);
	  }
	  else delete *Z;
	}
	delete path_relinking_sols;
      }
    }
  }
  now = clock();

  logInfo
    (cout
     << "\t***\tPath Relinking results\t***" << endl
     << "   Best Response time :\t" << best_rt << endl
     << "Average Response time :\t" << avg_rt / N << endl
     << "  Worst Response time :\t" << worst_rt << endl
     << "           time (sec) :\t" << (double)(now - beginning)/CLOCKS_PER_SEC << endl
     << "   Improved solutions :\t" << total_improved_solutions 
     << " (" << N << ")" << endl
     );

  results << N << "," << total_improved_solutions << "," 
	  << best_rt << "," << avg_rt << "," << worst_rt << ",";

  Best_input = SQM_best_solution(Solutions);
  /* Clears Solutions */
  Best = SQM_leave_only_the_best(improved_solutions);
  logDebug(cout << tag << "Improved solutions deleted" << endl);

  logInfo
    (cout << tag << "Diference between best input and best output :\t"
     << 100 * (Best_input->get_response_time() - Best->get_response_time()) /
     Best_input->get_response_time() << " %" << endl
     );
  logDebug(cout << tag << "Finish" << endl);
  logLevel = LOG_INFO;
  return Best;
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

double SQM_min_cost_pm(list<SQM_solution*> *Sols,SQM_solution *Sol) {
  double cost,min_cost;
  list<SQM_solution*>::iterator it = Sols->begin();

  min_cost = PR_perfect_matching_cost(*it,Sol);
  while (it != Sols->end()) {
    cost = PR_perfect_matching_cost(*it,Sol);
    if (min_cost > cost) min_cost = cost;
    it++;
  }

  return min_cost;
}
