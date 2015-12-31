
#include "PerfectMatching.h"

int* Perfect_Matching(int p,double **C) {

  MinCost<int,double> G(2*p + 2,p * (p+2),NULL);

  /* Add edges from Source Node */
  G.AddNodeExcess(0,p);
  for (int i = 0;i < p;i++)
    G.AddEdge(0,i+1,1,0,0.0);

  /* Add edges between elements */
  for (int i = 0;i < p;i++)
    for (int j = 0;j < p;j++)
      G.AddEdge(i+1,p+j+1,1,0,C[i][j]);
    
  /* Add edges to Sink node */
  G.AddNodeExcess(2*p+1,-p);
  for (int j = 0;j < p;j++)
    G.AddEdge(p+j+1,2*p+1,1,0,0.0);

  G.Solve();

  int *pm,e;
  pm = new int [p];
  for (int i = 0;i < p;i++) {
    for (int j = 0;j < p;j++) {
      e = p * (i+1) + j;
      if (G.GetRCap(e) == 0) {
	pm[i] = j;
	break;
      }
    }
  }

  return pm;
}

double Perfect_Matching_cost(int p,double **C) {
  double cost = 0.0;
  int *pm = Perfect_Matching(p,C);

  for (int i = 0;i < p;i++)
    cost += C[i][pm[i]];

  delete [] pm;
  return cost;
}
