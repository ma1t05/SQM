
#include "SQM_heuristic.h"

int unif(int);
int comp(const void*,const void*);
void sort_dist (int,double*,int*);
response_unit* guess_a_location_01(int,int,point*); // Returns the first p
response_unit* guess_a_location_02(int,int,point*); // Returns p random with replace
response_unit* guess_a_location_03(int,int,point*); // Returns p random

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
  /* */
  double *Lambda;
  double delta_mu;
  int **a;
  int m = I->M; /* Number of demand points */
  int n = I->N; /* Number of potencial sites to locate a server*/
  point *V = I->V,*W = I->W;
  /* Variables por ajustar */
  double v = 64.0;
  double beta = 1.5;
 
  /* Guess a location */
  X = guess_a_location_03(p,n,W);
  for (int i = 0;i < p;i++) {
    X[i].v = v;
    X[i].beta = beta;
  }

  /* Populate matrix of distances */
  Dist = new double*[m];
  for (int j = 0;j < m;j++)
    Dist[j] = new double [n];

  for (int j = 0;j < m;j++) {
    for (int i = 0;i < n;i++) {
      Dist[j][i] = dist(&(V[j]),&(W[i]));
    }
  }

  t_r = 99999.9;
  MST = new double[p];
  mst = new double[p];
  d = new double[p];
  Lambda = new double[m];
  f = new double*[p];
  for (int i = 0;i < p;i++)
    f[i] = new double[m];
  Tao = new double*[p];
  for (int i = 0;i < p;i++)
    Tao[i] = new double[m];

  a = new int*[m];
  for (int k = 0;k < m;k++)
    a[k] = new int[p];

  double demand = 0.0;
  for (int k = 0;k < m;k++) demand += I->V[k].demand;
  for (int k = 0;k < m;k++)
    Lambda[k] = I->V[k].demand * lambda / demand;


  for (int i = 0;i < p;i++)
    LogFile << X[i].location << " ";
  LogFile << endl;
  /* SERVICE MEAN TIME CALIBRATION */
  do {

    /* Print Current solution */
    for (int i = 0;i < p;i++)
      cout << X[i].location << " ";
    cout << "\r";

    T_r = t_r;

    /* **Step 0**
       Initialize Mean Service Time */
    for (int i = 0;i < p;i++)
      mst[i] = 1 / Mu_NT;

    do {

      for (int i = 0;i < p;i++)
	MST[i] = mst[i];

      /* **Step 1**:
	 Run the Hypercube Model */
      /* Update matrix of response times */
      for (int i = 0;i < p;i++) {
	for (int k = 0;k < m;k++) {
	  Tao[i][k] = (1 / Mu_NT + (X[i].beta / X[i].v) * Dist[k][X[i].location]/(60*24));
	}
      }
      
      /* Update matrix of pfreferred servers */
      for (int k = 0;k < m;k++) {
	for (int i = 0;i < p;i++)
	  d[i] = Dist[k][X[i].location];
	sort_dist(p,d,a[k]);
	/*
	cout << k << ":";
	for (int i = 0;i < p;i++) cout << " " << a[k][i];
	cout << endl;
	*/
      }

      jarvis_hypercube_approximation(m,p,Lambda,Tao,a,f);

      /* *Expected Response Time* */
      /* + expected travel time component */
      t_r = 0.0;
      for (int i = 0;i < p;i++) {
	for (int k = 0;k < m;k++)
	  t_r += f[i][k] * Dist[k][X[i].location];
      }

      P_B0 = 1.0;
      for (int i = 0;i < p;i++) {
	double rho_i = 0.0;
	for (int k = 0;k < m;k++)
	  rho_i += Lambda[k] * f[i][k] * Tao[i][k];
	P_B0 *= (1 - rho_i);
      }

      /* + mean queue delay component */
      mu = 0.0;
      for (int i = 0;i < p;i++)
	mu += 1 / MST[i];
      t_r += P_B0 * mu / pow(mu - lambda,2.0);
      
      /* Step 2 
	 Update mean service time */
      for (int i = 0;i < p;i++) {
	double h = 0.0;
	for (int k = 0;k < m;k++)
	  h += f[i][k];
	mst[i] = 0.0;
	for (int k = 0;k < m;k++)
	  mst[i] += (f[i][k]/h) * (1 / Mu_NT + (X[i].beta / X[i].v) * Dist[k][X[i].location]/(60*24));
      }
      
      /* Step 3 */
      delta_mu = 0.0;
      for (int i = 0;i < p;i++) {
	if (abs(mst[i] - MST[i]) > delta_mu)
	  delta_mu = abs(mst[i] - MST[i]);
      }

    } while (delta_mu > epsilon);
    
    double *h_i = new double[m];
    for (int i = 0;i < p;i++) {
      // Block A
      //cout << "\t\tBlock A " << i << endl;
      double h = 0.0;
      for (int k = 0;k < m;k++) 
	h += f[i][k];
      if (h == 0) cout << "for i = " << i+1 << " sum over f_ij is 0" << endl;
      for (int k = 0;k < m;k++) 
	h_i[k] = f[i][k]/h;

      // Block B
      /* Solve te 1-median location model with h_i^j */
      int best_location = -1;
      double best_sol,sol;
      for (int j = 0;j < n;j++) {
	sol = 0.0;
	for (int k = 0;k < m;k++) 
	  sol += h_i[k] * Dist[k][j];
	if (best_location == -1 || sol < best_sol) {
	  best_location = j;
	  best_sol = sol;
	}
      }
      X[i].location = best_location;

      /* Print current solution to LogFile */
      for (int i = 0;i < p;i++) {
	LogFile << X[i].location;
	if (X[i].past_location != X[i].location)
	  LogFile << "*";
	LogFile << "\t";
	X[i].past_location = X[i].location;
      }
      LogFile << endl;

    }
    delete [] h_i;

  } while (abs(T_r - t_r) > epsilon);
  
  for (int k = 0;k < n;k++) delete a[k];
  delete [] a;
  for (int i = 0;i < p;i++) delete [] Tao[i];
  delete [] Tao;
  delete [] Lambda;
  delete [] d;
  delete [] mst;
  delete [] MST;
  for (int j=0;j < m;j++) delete [] Dist[j];
  delete [] Dist;
  for (int i = 0;i < p;i++) delete [] f[i];
  delete [] f;

  for (int k = 0;k < n;k++)
    for (int i = 0;i < p;i++)
      if (X[i].location == k) cout << k << " ";
  cout << endl;
  return X;
}

int unif(int a) {
  return floor(double(a) * rand() / RAND_MAX);
}

int comp(const void *a,const void *b) {
  std::pair<double,int> *x,*y;
  x = (std::pair<double,int>*)a;
  y = (std::pair<double,int>*)b;
  if (x->first > y->first) return 1;
  else if (x->first < y->first) return -1;
  else return 0;
}

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

response_unit* guess_a_location_01(int p,int n, point *W){
  response_unit *X;
  X = new response_unit[p];
  for (int i = 0;i < p;i++) {
    X[i].location = i;
    X[i].past_location = i;
  }
  return X;
}

response_unit* guess_a_location_02(int p,int n, point *W){
  response_unit *X;
  for (int i = 0;i < p;i++) {
    X[i].location = unif(n);
  }
  return X;
}

response_unit* guess_a_location_03(int p,int n, point *W){
  response_unit *X;
  bool *option;
  int location,locations = 0;
  X = new response_unit[p];
  option = new bool[n];
  for (int i = 0;i < n;i++) option[i] = false;
  do {
    location = unif(n);
    //cout << "location = " << location << endl;
    if (option[location] == false) {
      locations++;
      option[location] = true;
    }
  } while (locations < p);
  for (int i = 0;i < n;i++) {
    if (option[i]) {
      X[--p].location = i;
      //cout << i << " ";
    }
  }
  //cout << endl;
  delete [] option;
  return X;
}
