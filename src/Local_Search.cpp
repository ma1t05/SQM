
#include <iostream>
#include <list>
#include "Local_Search.h"

typedef list<int> Sites;
typedef list<int> Servers;
/* move less workload server near to more workload */
int LS_movement_lm(SQM_solution&);
/* move a adyacent server to the server with moreworkload closer to him */
void LS_movement_mh(SQM_solution&); 
int LS_get_server_with_less_workload(SQM_solution&);
int LS_get_server_with_less_workload(SQM_solution&,Sites*);
int LS_get_server_with_more_workload(SQM_solution&);
Sites* LS_get_adjacent_sites(SQM_solution&,int);
Servers* LS_get_adjacent_servers(SQM_solution&,int);
void LS_print_workloads(SQM_solution&);

void Local_Search (SQM_solution &X) {
  int j;
  double rt;
  do {
    rt = X.get_response_time ();
    j = LS_movement_lm(X);
  } while (X.get_response_time () < rt);
  X.set_server_location(j,X.get_server_past_location(j));
}

int LS_movement_lm(SQM_solution &X) {
  int best_loc = UNASIGNED_LOCATION;
  double best_rt;
  int p = X.get_servers();
  /* Obtain the server with less workload */
  int j = LS_get_server_with_less_workload(X);
  int loc_j  = X.get_server_location(j);
  double v = X.get_server_speed(j);
  double beta = X.get_server_beta(j);
  /* obtain the news workloads */
  int k = LS_get_server_with_more_workload(X);
  /* Put a server near the server with more workload */
  Sites *lst = LS_get_adjacent_sites(X,k);
  for (Sites::iterator it = lst->begin(),end = lst->end(); it != end;it++) {
    X.test_server_location(j,*it);
    if (best_loc == UNASIGNED_LOCATION || X.get_response_time() < best_rt) {
      best_loc = *it;
      best_rt = X.get_response_time();
    }
  }

  if (best_loc == UNASIGNED_LOCATION) {
    cerr << "Warining: No new location" << endl;
    best_loc = loc_j;
  }
  X.test_server_location(j,loc_j); /* The past location of the server */
  X.set_server_location(j,best_loc);
  delete lst;
  return j;
}

void LS_movement_mh(SQM_solution &X) {
  int best_location;
  double best_rt;
  /* Obtain the server with more workload */
  int k = LS_get_server_with_more_workload(X);
  /* Obtain the adyacent servers */
  Servers *servers = LS_get_adjacent_servers(X,k);
  /* Obtain the adyacent server with less workload */
  int j = LS_get_server_with_less_workload(X,servers);
  int loc_j = X.get_server_location(j);
  delete servers;
  /* Put a server near the server with more workload */
  Sites *sites = LS_get_adjacent_sites(X,k);
  for (Sites::iterator it = sites->begin(),end = sites->end(); it != end;it++) {
    X.test_server_location(j,*it);
    if (best_location == UNASIGNED_LOCATION ||
	X.get_response_time() < best_rt) {
      best_location = *it;
      best_rt = X.get_response_time();
    }
  }
  X.test_server_location(j,loc_j); /* The past location of the server */
  X.set_server_location(j,best_location);
  delete sites;
}

int LS_get_server_with_less_workload(SQM_solution &X) {
  int j;
  int p = X.get_servers();
  double *wl;
  wl = X.get_workload();
  j = 0;
  for (int i = 1;i < p;i++)
    if (wl[i] < wl[j])
      j = i;
  delete [] wl;
  return j;
}

int LS_get_server_with_less_workload(SQM_solution &X,Servers *servers) {
  int j;
  int p = X.get_servers();
  double *wl;
  wl = X.get_workload();
  j = servers->front();
  for (Servers::iterator it = servers->begin(), end = servers->end();it != end;it++)
    if (wl[*it] < wl[j])
      j = *it;
  delete [] wl;
  return j;
}

int LS_get_server_with_more_workload(SQM_solution &X) {
  int j;
  int p = X.get_servers();
  double *wl;
  wl = X.get_workload();
  j = 0;
  for (int i = 1;i < p;i++)
    if (wl[i] > wl[j])
      j = i;
  delete [] wl;
  return j;
}

Sites* LS_get_adjacent_sites(SQM_solution &X,int i) {
  SQM_instance *I = X.get_instance();
  int m = I->demand_points();
  int n = I->potential_sites();
  int p = X.get_servers();
  bool *adjacent;
  Sites *lst;
  int loc_i,loc_j;
  int nearest_loc = UNASIGNED_LOCATION;
  int **a = X.preferred_servers();
  int sites,order;

  logDebug(cout << "LS_get_adjacent_sites: Start" << endl);
  lst = new Sites;
  adjacent = new bool [n];
  for (int k = 0;k < n;k++) adjacent[k] = false;

  loc_i = X.get_server_location(i);
  sites = 0;
  order = 0;
  do {
    for (int j = 0;j < m;j++)
      if (a[j][order] == i) {
	double radious = I->distance(loc_i,j);
	for (int loc = 0;loc < n;loc++)
	  if (I->distance(loc,j) <= radious)
	    if (!adjacent[loc]) {
	      adjacent[loc] = true;
	      sites++;
	    }
      }
    order++;
  } while (sites == 0 && order < p);

  for (int loc = 0;loc < n;loc++)
    if (adjacent[loc])
      lst->push_back(loc);
  
  delete [] adjacent;
  logDebug(cout << "LS_get_adjacent_sites: Finish" << endl);
  return lst;
}

/* Return adjacent servers acording delaunay triangulation */
Servers* LS_get_adjacent_servers(SQM_solution &X,int j) {
  Servers *lst;
  
  return lst;
}
