
#include "SQM_heuristic.h"

int comp(const void *a,const void *b);

void sort_dist (int n,double *d,int *c) {
  std::pair<double,int> *x;
  x = new std::pair<double,int>[n];
  for (int i = 0;i < n;i++) {
    x[i].first = d[i];
    x[i].second = i;
  }
  qsort(x,n,sizeof(std::pair<double,int>), comp);
  for (int i = 0;i < n;i++)
    c[i] = x[i].second;
  delete[] x;
}

response_unit* SQM_heuristic
(SQM_instance *I,
 int p, // Number of adjusters
 double lambda, // mean rate per unit of time within service calls are generated in Poisson manner
 double Mu_NT // mean of non-travel time component of the service time
 ) {
  
  /* Variable definitions */
  response_unit *X;
  double *MST,*mst; // mean service time
  double T_r,t_r; // expected response time
  double **Dist; // Matrix of distances
  double **f,**Tao;
  double *d;
  double mu;
  double P_B0; // Pendiente
  double *Lambda;
  double delta_mu;
  int **a;
  int n = I->N,m = I->M;
  point *V = I->V,*W = I->W;
  
  /* Populate matrix of distances */
  Dist = new double*[m];
  for (int j = 0;j < m;j++)
    Dist[j] = new double [n];

  for (int j = 0;j < m;j++) {
    for (int i = 0;i < n;i++) {
      Dist[j][i] = dist(&(V[i]),&(W[j]));
    }
  }

  /* Guess a location */
  X = new response_unit[p];
  for (int i = 0;i < p;i++)
    X[i].location = i;

  T_r = 99999.9;
  MST = new double[p];
  mst = new double[p];
  d = new double[p];
  Lambda = new double[n];
  f = new double*[p];
  for (int i = 0;i < p;i++)
    f[i] = new double[n];
  Tao = new double*[p];
  for (int i = 0;i < p;i++)
    Tao[i] = new double[n];

  a = new int*[n];
  for (int k = 0;k < n;k++)
    a[k] = new int[p];

  for (int k = 0;k < n;k++)
    Lambda[k] = I->V[k].demand * lambda;

  /* SERVICE MEAN TIME CALIBRATION */
  do {
    // Step 0
    for (int i = 0;i < p;i++)
      MST[i] = 1 / Mu_NT;
    do {

      // Step 1: Run the Hypercube Model

      for (int i = 0;i < p;i++) {
	for (int k = 0;k < n;k++) {
	  Tao[i][k] = (Mu_NT + (X[i].beta / X[i].v) * Dist[X[i].location][k]);
	}
      }
      
      for (int k = 0;k < n;k++) {
	for (int i = 0;i < p;i++)
	  d[i] = Dist[X[i].location][k];
	sort_dist(n,d,a[k]);
      }

      f = jarvis_hypercube_approximation(n,p,Lambda,Tao,a);

      // T_R(X)
      t_r = 0.0;
      // the expected travel time component
      for (int i = 0;i < p;i++) {
	for (int k = 0;k < m;k++)
	  t_r += f[i][k] * Dist[X[i].location][k];
      }
      // the mean queue delay component
      mu = 0.0;
      for (int i = 0;i < p;i++)
	mu += 1 / MST[i];
      t_r += P_B0 * mu / pow(mu - lambda,2.0);
      
      // Step 2
      for (int i = 0;i < p;i++) {
	double h = 0.0;
	for (int j = 0;j < n;j++)
	  h += f[i][j];
	mst[i] = 0.0;
	for (int k = 0;k < n;k++)
	  mst[i] = (f[i][k]/h) * (Mu_NT + (X[i].beta / X[i].v) * Dist[X[i].location][k]);
      }
      
      // Step 3
      delta_mu = 0.0;
      for (int i = 0;i < p;i++) {
	if (abs(mst[i] - MST[i]) > delta_mu)
	  delta_mu = abs(mst[i] - MST[i]);
      }

      if (delta_mu > epsilon) {
	for (int i = 0;i < p;i++)
	  MST[i] = mst[i];
      }
      
    } while (delta_mu > epsilon);
  } while (abs(T_r - t_r) > epsilon);
  
  for (int j=0;j < m;j++)
    delete [] Dist[j];
  delete [] Dist;

  return X;
}

void heuristic1
(int p, // number of response units
 network *G, // the transportation network
 double lambda, // mean rate per unit of time within service calls are generated in Poisson manner
 double Mu_NT, // mean of non-travel time component of the service time
 response_unit *X) {

  double *MST,*mst; // mean service time
  double T_r,t_r; // expected response time
  double **f;
  double **Dist; // Matrix of distances
  double mu;
  double P_B0;
  bool exit;
  int n = G->n;

  Dist = new double*[G->n];
  for (int i = 0;i < G->n;i++) {
    Dist[i] = new double[G->n];
  }
  for (int i = 0;i < G->n;i++) {
    for (int j = 0;j < G->n;j++) {

    }
  }
  
  MST = new double[p];
  mst = new double[p];
  f = new double*[p];
  for (int i = 0;i < p;i++)
    f[i] = new double[G->n];
  
  do {
    // Step 0
    for (int i = 0;i < p;i++)
      MST[i] = 1 / Mu_NT;
  
    do {
      // Step 1: Run the Hypercube Model
  
      // T_R(X)
      t_r = 0.0;
      // the expected travel time component
      for (int i = 0;i < p;i++) {
	for (int k = 0;k < G->n;k++)
	  t_r += f[i][k] * Dist[X[i].location][k];
      }
      // the mean queue delay component
      mu = 0.0;
      for (int i = 0;i < p;i++)
	mu += 1 / MST[i];
      t_r += P_B0 * mu / pow(mu - lambda,2.0);
      
      // Step 2
      for (int i = 0;i < p;i++) {
	double h = 0.0;
	for (int j = 0;j < G->n;j++)
	  h += f[i][j];
	mst[i] = 0.0;
	for (int k = 0;k < G->n;k++)
	  mst[i] = (f[i][k]/h) * (Mu_NT + (X[i].beta / X[i].v) * Dist[X[i].location][k]);
      }
  
      // Step 3
      exit = true;
      for (int i = 0;i < p;i++) {
	if (abs(mst[i] - MST[i]) > EPSILON) {
	  exit = false;
	  for (i = 0;i < p;i++)
	    MST[i] = mst[i];
	  break;
	}
      }
    } while (!exit);

    exit = true;
    if (abs(t_r - T_r) > EPSILON) {
      exit = false;
      for (int i = 0;i < p;i++) {

	double T = 0.0,t;
	double h = 0.0;
	double *hi;
	hi = new double[G->n];
	// h_j^i = f_{ij}\sum_{k=1}^{n}{f_{ij}}
	for (int j = 0;j < n;j++)
	  h += f[i][j];
	for (int k = 0;k < n;k++)
	  hi[k] = f[i][k] / h;

	// Solve the 1-median location model with h_j^i
	for (int j = 0;j < n;j++)
	  T += hi[j] * Dist[X[i].location][j];
	for (int k = 0;k < n;k++) {
	  t = 0;
	  for (int j = 0;j < n;j++)
	    t += hi[j] * Dist[X[i].location][j];
	  if (t < T) {
	    T = t;
	    X[i].location = t;
	  }
	}
	delete[] hi;

      }
      T_r = t_r;
    }
  } while (!exit);

  for (int i = 0;i < p;i++)
    delete[] f[i];
  delete[] f;
  delete[] mst;
  delete[] MST;
  for (int i = 0;i < p;i++)
    delete[] Dist[i];
  delete[] Dist;
  
}

int comp(const void *a,const void *b) {
  std::pair<double,int> *x,*y;
  x = (std::pair<double,int>*)a;
  y = (std::pair<double,int>*)b;
  if (x->first > y->first) return 1;
  else if (x->first < y->first) return -1;
  else return 0;
}
