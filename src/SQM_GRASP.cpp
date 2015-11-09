
#include "SQM_GRASP.h"
#include "mp_jarvis.h"
#include "MST.h"
#include "instance-creator.h"
#include "log.h"

bool GRASP_closest_to_b(SQM_instance *I,int node,int center_a,int center_b);
int GRASP_nearest_server(SQM_instance *I,int j,int p,response_unit *X);

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
	T_r[i] = GRASP_func_NN(I,r+1,X,lambda,Mu_NT);
	/* T_r[i] = MST_response_time(I,r+1,X,lambda,Mu_NT);*/
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
  logDebug(cout << "Start GRASP_func_NN" << endl);
  /* Variables definition */
  int n = I->N,m = I->M;
  int nearest,k;
  double *rho;
  double Obj;
  /* */

  Obj = 0.0;
  for (int j = 0;j < m;j++) {
    k = GRASP_nearest_server(I,j,p,X); /* Obtain the nearest server */
    nearest = X[k].location;
    Obj += dist(&(I->V[j]),&(I->W[nearest])) / X[k].v;
  }
  Obj /= MINS_PER_BLOCK * BLOCKS_PER_HORIZON;

  logDebug(cout << "Finish GRASP_func_NN" << endl);
  return Obj;
}

bool GRASP_closest_to_b(SQM_instance *I,int node,int center_a,int center_b) {
  return (dist(&(I->V[node]),&(I->W[center_a])) > dist(&(I->V[node]),&(I->W[center_b])));
}

int GRASP_nearest_server(SQM_instance *I,int j,int p,response_unit *X) {
  int k = 0;
  for (int i = 1;i < p;i++) {
    if (GRASP_closest_to_b(I,j,X[k].location,X[i].location)) {
      k = i;
    }
  }
  return k;
}

double GRASP_func_kNN
(SQM_instance *I,
 int p,
 response_unit *X,
 double lambda,
 double Mu_NT,
 int K
 ) {
  logDebug(cout << "Start GRASP_func_kNN" << endl);
  /* Variable definitions */
  int m = I->M; /* Number of demand points */
  int **a;
  double **Dist;
  double *d;

  a = new int*[m];
  for (int k = 0;k < m;k++)
    a[k] = new int[p];
  Dist = SQM_dist_matrix(I);
  d = new double[p];

  for (int k = 0;k < m;k++) {
    for (int i = 0;i < p;i++)
      d[i] = Dist[k][X[i].location];
    sort_dist(p,d,a[k]);
  }

  /* Pendiente: Calculo de Tao, normalizar lambda * demanda / demanda total */
  for (int i = 0; i < p;i++) {
    rho[i] = 0.0;
    for (int k = 0; k < m;k++)
      if (a[k][0] == i)
	rho[i] += Tao * lambda * I->V[k].demand;
  }

  delete [] d;
  for (int j = 0;j < m;j++) delete [] Dist[j];
  delete [] Dist;
  for (int k = 0;k < m;k++) delete [] a[k];
  delete [] a;
  logDebug(cout << "Finish GRASP_func_kNN" << endl);
}
 
