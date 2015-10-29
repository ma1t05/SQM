
#include "SQM_GRASP.h"

bool GRASP_closest_to_b(SQM_instance *I,int node,int center_a,int center_b);

response_unit* GRASP
(SQM_instance *I,
 int p, // Number of adjusters
 double lambda, // mean rate per unit of time within service calls are generated in Poisson manner
 double Mu_NT, // mean of non-travel time component of the service time
 double v, // Speed
 double alpha // Random factor {1: random, 0: greedy}
 ) {
  int n = I->N,m = I->M;
  int r;
  int element;
  int *rcl;
  double *T_r;
  double beta = 1.5;
  response_unit *X;

  /* Debug cout << endl << endl << "*****Start GRASP*****" << endl << endl << endl; /* */
  if (p < 1) return NULL;
  X = new response_unit[p];
  for (int i = 0;i < p;i++) {
    X[i].v = v;
    X[i].beta = beta;
  }

  // cout << "/* Locate the first server */" << endl;
  X[0].location = unif(n);
  r = 1;
  T_r = new double [n];
  rcl = new int [n];
  while (r < p) 
    {
      // cout << "[" << r << "]/* Evaluate posible locations*/" << "\t";
      for (int i = 0;i < n;i++) {
	X[r].location = i;
	T_r[i] = SQM_response_time(I,r+1,X,lambda,Mu_NT);
      }

      // cout << "/* Sort Restricted Candidates List */" << "\t";
      sort_dist(n,T_r,rcl);
      // cout << "/* Choose random element from the rcl */";
      // cout <<"\r";
      element = unif(ceil(alpha * n));
      X[r++].location = rcl[element];
    }
  // cout << endl;

  delete [] rcl;
  delete [] T_r;
  /* Debug cout << "Finish GRASP" << endl; /* */
  return X;
}

double GRASP_func_NN
(SQM_instance *I,
 int p,
 response_unit *X,
 double lambda,
 double Mu_NT
 ) {
  /* Variables definition */
  int n = I->N,m = I->M;
  int nearest,k;
  double *rho;
  double Obj;
  /* */

  Obj = 0.0;
  for (int j = 0;j < m;j++) {
    /* Obtain the nearest server */
    nearest = X[0].location;
    for (int i = 1;i < p;i++) {
      if (GRASP_closest_to_b(I,j,nearest,X[i].location)) {
	nearest = X[i].location;
	k = i;
      }
    }
    Obj += dist(&(I->V[j]),&(I->W[nearest])) / X[k].v;
  }
  Obj /= MINS_PER_BLOCK * BLOCKS_PER_HORIZON;

  return Obj;
}

bool GRASP_closest_to_b(SQM_instance *I,int node,int center_a,int center_b) {
  return (dist(&(I->V[node]),&(I->W[center_a])) > dist(&(I->V[node]),&(I->W[center_b])));
}