
#include "PerfectMatching.h"

int labeling(int p,double **C,int i,int *pm);

int* Perfect_Matching(int p,double **C) {

  /* Step 1 */
  double min;
  for (int i = 0;i < p;i++) {
    min = C[i][0];
    for (int j = 1;j < p;j++)
      if (C[i][j] < min)
	min = C[i][j];
    if (min > 0)
      for (int j = 1;j < p;j++)
	C[i][j] -= min;
  }

  /* Step 2 */
  for (int j = 1;j < p;j++) {
    min = C[0][j];
    for (int i = 1;i < p;i++)
      if (C[i][j] < min)
	min = C[i][j];
    if (min > 0)
      for (int i = 0;i < p;i++)
	C[i][j] -= min;
  }

  int *pm;
  pm = new int [p];
  for (int i = 0;i < p;i++) pm[i] = -1;
  int lab = labeling(p,C,0,pm);
  if (lab == p)
    return pm; /* pending: transpose pm */
  else {
    
  }

}

int labeling(int p,double **C,int i,int *pm) {
  if (i == p) return 0;
  int min = p-i,depth,op;
  for (int j = 0;j < p;j++)
    if (C[i][j] == 0 && pm[j] == -1) {
      pm[j] = i;
      depth = labeling(p,C,i+1,pm);
      if (depth == 0)
	return 0;
      else if (depth < min) 
	min = depth;
      pm[j] = -1;
    }
  return min;
}
