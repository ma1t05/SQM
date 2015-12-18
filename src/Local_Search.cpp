
#include <list>
#include "Local_Search.h"
#include "MST.h"

void LS_movement_lm(SQM_solution *X); /* move less workload server near to more workload */
void LS_movement_mh(SQM_solution *X); /* move a adyacent server to the server with moreworkload closer to him */

void Local_Search (SQM_solution *X) {
  int p = X->get_servers();
  double *wl;
  wl = MST_workload(X);
  /* Obtain the server with less workload */
  /* Removethem */
  /* obtain the news workloads */
  /* Put a server near the server with more workload */
}

void LS_movement_lm(SQM_solution *X) {
  int p = X->get_servers();
  double *wl;
  wl = MST_workload(X);
  /* Obtain the server with less workload */
  /* Removethem */
  /* obtain the news workloads */
  /* Put a server near the server with more workload */
}

void LS_movement_mh(SQM_solution *X) {
  /* Obtain the server with more workload */
  /* Obtain the adyacent servers */
  /*  */
  /* Put a server near the server with more workload */
}

int LS_return_server_with_less_workload(SQM_solution *X) {
  int j;
  int p = X->get_servers();
  double *wl;
  wl = MST_workload(X);
  j = 0;
  for (int i = 1;i < p;i++)
    if (wl[i] < wl[j])
      j = i;
  return j;
}

int LS_return_server_with_more_workload(SQM_solution *X) {
  int j;
  int p = X->get_servers();
  double *wl;
  wl = MST_workload(X);
  j = 0;
  for (int i = 1;i < p;i++)
    if (wl[i] < wl[j])
      j = i;
  return j;
}

list<int> LS_return_adjacent_servers(SQM_solution *X,int j) {
  
}
