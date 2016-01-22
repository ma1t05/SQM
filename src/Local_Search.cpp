
#include <iostream>
#include <list>
#include "Local_Search.h"

typedef list<int> Sites;
typedef list<int> Servers;
/* move less workload server near to more workload */
void LS_movement_lm(SQM_solution&);
/* move a adyacent server to the server with moreworkload closer to him */
void LS_movement_mh(SQM_solution&); 
int LS_get_server_with_less_workload(SQM_solution&);
int LS_get_server_with_less_workload(SQM_solution&,Sites*);
int LS_get_server_with_more_workload(SQM_solution&);
Sites* LS_get_adjacent_sites(SQM_solution&,int);
Servers* LS_get_adjacent_servers(SQM_solution&,int);
void LS_print_workloads(SQM_solution&);

void Local_Search (SQM_solution &X) {
  int p = X.get_servers();
  double rt;
  do {
    rt = X.get_response_time ();
    LS_movement_lm(X);
  } while (X.get_response_time () < rt);
  X.set_server_location(p-1,X.get_server_past_location(p-1));
}

void LS_movement_lm(SQM_solution &X) {
  int best_loc = UNASIGNED_LOCATION;
  double best_rt;
  int p = X.get_servers();
  /* Obtain the server with less workload */
  int j = LS_get_server_with_less_workload(X);
  int loc_j  = X.get_server_location(j);
  double v = X.get_server_speed(j);
  double beta = X.get_server_beta(j);
  /* Remove them */ 
  X.remove_server(j);
  /* obtain the news workloads */
  int k = LS_get_server_with_more_workload(X);
  /* Put a server near the server with more workload */
  Sites *lst = LS_get_adjacent_sites(X,k);
  X.add_server();
  X.set_speed(v,beta);
  for (Sites::iterator it = lst->begin(),end = lst->end(); it != end;it++) {
    X.test_server_location(p-1,*it);
    if (best_loc == UNASIGNED_LOCATION || X.get_response_time() < best_rt) {
      best_loc = *it;
      best_rt = X.get_response_time();
    }
  }
  X.test_server_location(p-1,loc_j); /* The past location of the server */
  X.set_server_location(p-1,best_loc);
  delete lst;
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
  int n = I->potential_sites();
  int p = X.get_servers();
  Sites *lst;
  lst = new Sites;
  int loc_i,loc_j;
  int nearest_loc = UNASIGNED_LOCATION;
  loc_i = X.get_server_location(i);
  for (int j = 0;j < p;j++)
    if (j != i) {
      loc_j = X.get_server_location(j);
      if (nearest_loc == UNASIGNED_LOCATION || 
	  I->sites_distance(loc_i,loc_j) < I->sites_distance(loc_i,nearest_loc))
	nearest_loc = loc_j;
    }

  double radious = I->sites_distance(loc_i,nearest_loc);
  for (int k = 0;k < n;k++)
    if (I->sites_distance(loc_i,k) < radious + JARVIS_EPSILON)
      lst->push_back(k);
  return lst;
}

/* Return adjacent servers acording delaunay triangulation */
Servers* LS_get_adjacent_servers(SQM_solution &X,int j) {
  Servers *lst;
  
  return lst;
}
