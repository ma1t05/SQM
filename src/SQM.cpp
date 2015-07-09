
#include "SQM.h"
#include "SQM_model.h"

int main() {
  instance *I;
  I = read_points("./../../PMCLAP/Instancias/Q_MCLP_30.txt");
  SQM_model(I,4,mu,1.0,40.0);

  delete[] I->points;
  delete I;
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
  double **dist; // Matrix of distances
  double mu;
  double P_B0;
  bool exit;
  int n = G->n;

  dist = new double*[G->n];
  for (int i = 0;i < G->n;i++) {
    dist[i] = new double[G->n];
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
	  t_r += f[i][k] * dist[X[i].location][k];
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
	  mst[i] = (f[i][k]/h) * (Mu_NT + (X[i].beta / X[i].v) * dist[X[i].location][k]);
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
	  T += hi[j] * dist[X[i].location][j];
	for (int k = 0;k < n;k++) {
	  t = 0;
	  for (int j = 0;j < n;j++)
	    t += hi[j] * dist[X[i].location][j];
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
    delete[] dist[i];
  delete[] dist;
  
}

