
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

