
#deinfe epsilon 0.001

struct response_unit {
  int location;
  double v;
  double beta;
};

typedef struct response_unit response_unit;

int *heuristic1
(int p, // number of response units
 network G, // the transportation network
 double lambda, // mean rate per unit of time within service calls are generated in Poisson manner
 double Mu_NT, // mean of non-travel time component of the service time
 response_unit *X) {

  double *MST; // mean service time
  double T_r; // expected response time
  double f;
  bool exit;
  
  MST = new dobule[p];
  f = new double*[p];
  for (int i = 0;i < p;i++)
    f[i] = new dobule[G->n];

  do {
    for (int i = 0;i < p;i++)
      MST[i] = 1 / Mu_NT;
  
    do {
      // Step 1: Run the Hypercube Model
  
      // Step 2
      for (int i = 0;i < p;i++) {
	double h = 0.0;
	for (int j = 0:j < G->n;j++)
	  h += f[i][j];
	mst[i] = 0.0;
	for (int k = 0;k < G->n;k++)
	  mst[i] = (f[i][k]/h) * (Mu_NT + (X[i]->beta / X[i]->v) * dist[i][k]);
      }
  
      // Step 3
      exit = true;
      for (int i = 0;i < p;i++) {
	if (abs(mst[i] - MST[i]) > epsilon) {
	  exit = false;
	  for (i = 0;i < p;i++)
	    MST[i] = mst[i];
	  break;
	}
      }
    } while (!exit);

    exit = true;
    if (abs(t_r - T_r) > epsilon) {
      exit = false;
      for (int i = 0;i < p;i++) {
	// h_j^i = f_{ij}\sum_{k=1}^{n}{f_{ij}}

	// Solve the 1-median location model with h_j^i
      }
      T_r = t_r;
    }
  } while (!eixt);

  for (int i = 0;i < p;i++)
    delete[] f[i];
  delete[] f;
  delete[] MST;
  return X;
}
