
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
  num *MST,*mst; // mean service time
  num T_r,t_r; // expected response time
  double **Dist; // Matrix of distances
  num **f,**Tao;
  double *d;
  double mu;
  double P_B0; // Pendiente
  /* */
  num tmp;
  num delta_mu;
  num *Lambda;
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
  Dist = new  double*[m];
  for (int j = 0;j < m;j++)
    Dist[j] = new double [n];

  for (int j = 0;j < m;j++) {
    for (int i = 0;i < n;i++) {
      Dist[j][i] = dist(&(V[j]),&(W[i]));
    }
  }

  MST = new num[p];
  mst = new num[p];
  d = new double[p];
  Lambda = new num[m];
  Tao = new num*[p];
  for (int i = 0;i < p;i++)
    Tao[i] = new num[m];
  f = new num*[p];
  for (int i = 0;i < p;i++)
    f[i] = new num[m];

  mpf_init(tmp);
  mpf_init(delta_mu);
  for (int i = 0;i < p;i++)
    mpf_init(MST[i]);
  for (int i = 0;i < p;i++)
    mpf_init(mst[i]);
  for (int k = 0;k < m;k++)
    mpf_init(Lambda[k]);
  for (int i = 0;i < p;i++) {
    for (int k = 0;k < m;k++)
      mpf_init(Tao[i][k]);
  }
  for (int i = 0;i < p;i++) {
    for (int k = 0;k < m;k++)
      mpf_init(f[i][k]);
  }
  mpf_init(T_r);
  mpf_init_set_d(t_r,99999.9);

  a = new int*[m];
  for (int k = 0;k < m;k++)
    a[k] = new int[p];

  double demand = 0.0;
  for (int k = 0;k < m;k++) demand += I->V[k].demand;
  for (int k = 0;k < m;k++)
    mpf_set_d(Lambda[k],I->V[k].demand * lambda / demand);


  for (int i = 0;i < p;i++)
    LogFile << X[i].location << " ";
  LogFile << endl;
  /* SERVICE MEAN TIME CALIBRATION */
  do {

    /* Print Current solution */
    for (int i = 0;i < p;i++)
      cout << X[i].location << " ";
    cout << endl;

    mpf_set(T_r,t_r);

    /* **Step 0**
       Initialize Mean Service Time */
    for (int i = 0;i < p;i++)
      mpf_set_d(mst[i],1 / Mu_NT);

    do {

      for (int i = 0;i < p;i++)
	mpf_set(MST[i],mst[i]);

      /* **Step 1**:
	Run the Hypercube Model */
      /* Update matrix of response times */
      for (int i = 0;i < p;i++) {
	for (int k = 0;k < m;k++) {
	  mpf_set_d(Tao[i][k],1 / Mu_NT + (X[i].beta / X[i].v) * Dist[k][X[i].location]/(60*24));
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
      mpf_set_ui(t_r,0);
      for (int i = 0;i < p;i++) {
	for (int k = 0;k < m;k++) {
	  mpf_set_d(tmp,Dist[k][X[i].location]);
	  mpf_mul(tmp,tmp,f[i][k]);
	  mpf_add(t_r,t_r,tmp);
	}
      }

      P_B0 = 1.0;
      for (int i = 0;i < p;i++) {
	double rho_i = 0.0;
	for (int k = 0;k < m;k++) {
	  mpf_mul(tmp,Lambda[k],f[i][k]);
	  mpf_mul(tmp,tmp,Tao[i][k]);
	  rho_i += mpf_get_d(tmp);
	  P_B0 *= (1 - rho_i);
	}
      }
      /* + mean queue delay component */
      mu = 0.0;
      for (int i = 0;i < p;i++)
	mu += 1 / mpf_get_d(MST[i]);
      mpf_set_d(tmp,P_B0 * mu / pow(mu - lambda,2.0));
      mpf_add(t_r,t_r,tmp);
      
      /* Step 2 
	 Update mean service time */
      for (int i = 0;i < p;i++) {
	double h = 0.0;
	for (int k = 0;k < m;k++)
	  h += mpf_get_d(f[i][k]);
	mpf_set_ui(mst[i],0);
	for (int k = 0;k < m;k++) {
	  mpf_set_d(tmp,(mpf_get_d(f[i][k])/h) * (1 / Mu_NT + (X[i].beta / X[i].v) * Dist[k][X[i].location]/(60*24)));
	  mpf_add(mst[i],mst[i],tmp);
	}
      }
      
      /* Step 3 */
      mpf_set_ui(delta_mu,0);
      for (int i = 0;i < p;i++) {
	mpf_sub(tmp,mst[i],MST[i]);
	mpf_abs(tmp,tmp);
	if (mpf_cmp(tmp,delta_mu) > 0)
	  mpf_set(delta_mu,tmp);
      }

      cout << "Delta in mst: " << mpf_get_d(delta_mu) << endl;
    } while (mpf_cmp_d(delta_mu,epsilon) > 0);
    
    double *h_i = new double[m];
    num h;
    mpf_init(h);
    for (int i = 0;i < p;i++) {
      // Block A
      //cout << "\t\tBlock A " << i << endl;
      mpf_set_ui(h,0);
      for (int k = 0;k < m;k++) 
	mpf_add(h,h,f[i][k]);
      if (mpf_cmp_ui(h,0) == 0)
	cout << "for i = " << i+1 << " sum over f_ij is 0" << endl;
      for (int k = 0;k < m;k++) {
	mpf_div(tmp,f[i][k],h);
	h_i[k] = mpf_get_d(tmp);
      }

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
    }

    /* Print current solution to LogFile */
    for (int i = 0;i < p;i++) {
      LogFile << X[i].location;
      if (X[i].past_location != X[i].location)
	LogFile << "*";
      LogFile << "\t";
      X[i].past_location = X[i].location;
    }
    LogFile << endl;

    mpf_clear(h);
    delete [] h_i;

    mpf_sub(tmp,T_r,t_r);
    mpf_abs(tmp,tmp);

  } while (mpf_cmp_d(tmp,epsilon) > 0);
  
  for (int k = 0;k < n;k++) delete a[k];
  delete [] a;

  mpf_clear(T_r);
  mpf_clear(t_r);
  for (int i = 0;i < p;i++) {
    for (int j = 0;j < n;j++)
      mpf_clear(f[i][j]);
  }
  for (int i = 0;i < p;i++) {
    for (int j = 0;j < n;j++)
      mpf_clear(Tao[i][j]);
  }
  for (int k = 0;k < m;k++)
    mpf_clear(Lambda[k]);
  mpf_clear(delta_mu);
  mpf_clear(tmp);

  for (int i = 0;i < p;i++)
    delete [] Tao[i];
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
      X[p].past_location = i;
    }
  }
  //cout << endl;
  delete [] option;
  return X;
}
