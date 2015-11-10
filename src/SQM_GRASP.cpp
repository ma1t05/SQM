
#include "SQM_GRASP.h"
#include "mp_jarvis.h"
#include "MST.h"
#include "instance-creator.h"
#include "log.h"

#define min(A,B) ((A)>(B)?(A):(B))

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
	/* T_r[i] = MST_response_time(I,r+1,X,lambda,Mu_NT);*/
	/* T_r[i] = GRASP_func_NN(I,r+1,X,lambda,Mu_NT);*/
	T_r[i] = GRASP_func_kNN(I,r+1,X,lambda,Mu_NT,min(r+1,3));
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
    Obj += I->V[j].demand * dist(&(I->V[j]),&(I->W[nearest])) / X[k].v;
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
  double RT = 0.0; /* Response Time */
  int **a;
  double **Dist;
  double *Lambda;
  double demand;
  logInfo(cout << "Termina definicion de variables" << endl);

  a = new int*[m];
  for (int k = 0;k < m;k++)
    a[k] = new int[p];
  Dist = SQM_dist_matrix(I);

  double *d = new double[p];
  for (int k = 0;k < m;k++) {
    for (int i = 0;i < p;i++)
      d[i] = Dist[k][X[i].location];
    sort_dist(p,d,a[k]);
  }
  delete [] d;

  Lambda = new double[m];
  demand = 0.0;
  for (int k = 0;k < m;k++) demand += I->V[k].demand;
  for (int k = 0;k < m;k++) Lambda[k] = I->V[k].demand * lambda / demand;

  logDebug(cout << "Comienza calculo de rho_i" << endl);
  /* Calculate the first approach for rho */
  double distance;
  double *rho = new double[p];
  for (int i = 0;i < p;i++) rho[i] = 0.0;
  for (int i = 0;i < p;i++) {
    distance = 0.0;
    for (int k = 0;k < m;k++)
      if (a[k][0] == i) {
	rho[i] += Lambda[k] * (1/Mu_NT + (X[i].beta / X[i].v) * Dist[k][X[i].location]/ (MINS_PER_BLOCK * BLOCKS_PER_HORIZON));
	distance += Dist[k][X[i].location];
      }
    RT += rho[i] * distance / (X[i].v * MINS_PER_BLOCK * BLOCKS_PER_HORIZON);
  }
  
  logDebug(cout << "Comienza calculo de rho_i para resto de ordenes" << endl);
  int t = 1;
  double rho_a_ml;
  double *new_rho = new double [p];
  do {
    for (int i = 0;i < p;i++) new_rho[i] = 0.0;
    for (int i = 0;i < p;i++) {
      distance = 0.0;
      for (int k = 0;k < m;k++)
	if (a[k][t] == i) {
	  rho_a_ml = 1;
	  for (int l = 0;l < t;l++) rho_a_ml *= rho[a[k][l]];
	  new_rho[i] += (1 - rho[i]) * rho_a_ml * Lambda[k] * (1/Mu_NT + (X[i].beta / X[i].v) * Dist[k][X[i].location] / (MINS_PER_BLOCK * BLOCKS_PER_HORIZON));
	  distance += Dist[k][X[i].location];
	}
      RT += new_rho[i] * distance / (X[i].v * MINS_PER_BLOCK * BLOCKS_PER_HORIZON);
      rho[i] += new_rho[i];
      if (rho[i] > 1)
	logError(cout << "¡rho_" << i+1 << " > 1! in order " << t << endl);
    }
  } while (++t < K &&  t < p);
  logDebug(cout << "Comienza a liberar memoria" << endl);

  delete [] new_rho;    
  delete [] Lambda;
  for (int j=0;j < m;j++) delete [] Dist[j];
  delete [] Dist;
  for (int k = 0;k < m;k++) delete [] a[k];
  delete [] a;
  logDebug(cout << "Finish GRASP_func_kNN" << endl);
  return RT;
}
