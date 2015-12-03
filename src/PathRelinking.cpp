
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
double** PR_distances_matrix(SQM_solution*,SQM_solution*);

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
  int *pm;
  double **distances;
  distances = PR_distances_matrix(X,Y);
  pm = Perfect_Matching(p,distances);
  for (int i = 0;i < p;i++) delete [] distances[i];
  delete [] distances;
  return pm;
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

SQM_solution* SQM_path_relinking(list<SQM_solution*>* Solutions) {
  int total_improved_solutions = 0;
  double Tr_x,Tr_y,Best_TR,TR;
  list<SQM_solution*>::iterator X,Y,Z;
  list<SQM_solution*> *path_relinking_sols,*improved_solutions;
  SQM_solution *Best,*Best_input;
  clock_t beginning,now;
  double best_rt,avg_rt,worst_rt;
  int N;

  beginning = clock();
  improved_solutions = new list<SQM_solution*>;
  best_rt = 100.0,worst_rt = 0.0,avg_rt = 0.0, N = 0;
  for (X = Solutions->begin();X != Solutions->end();X++) {
    Tr_x = (*X)->get_response_time();
    for (Y = X;Y != Solutions->end();Y++) {
      if (Y != X) {
	Tr_y = (*Y)->get_response_time();
	Best_TR = (Tr_x > Tr_y ? Tr_y : Tr_x);
	path_relinking_sols = Path_Relinking(*X,*Y);
	if (path_relinking_sols != NULL) {
	  for (Z = path_relinking_sols->begin();Z != path_relinking_sols->end();Z++) {
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
  }
  now = clock();

  cout << "\t***\tPath Relinking results\t***" << endl;
  cout << "   Best Response time :\t" << best_rt << endl
       << "Average Response time :\t" << avg_rt / N << endl
       << "  Worst Response time :\t" << worst_rt << endl
       << "           time (sec) :\t" << (double)(now - beginning)/CLOCKS_PER_SEC << endl
       << "   Improved solutions :\t" << total_improved_solutions 
       << " (" << N << ")" << endl;


  Best_input = SQM_best_solution(Solutions);
  /* Clears Solutions */
  Best = SQM_leave_only_the_best(improved_solutions);
  logDebug(cout << "Improved solutions deleted" << endl);

  logInfo(cout << "Diference between best input and best output :\t"
	  << 100*(Best_input->get_response_time() - Best->get_response_time()) / Best_input->get_response_time() << " %" << endl);

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
