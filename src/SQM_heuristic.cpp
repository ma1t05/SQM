
#include "SQM_heuristic.h"

int unif(int);
int comp(const void*,const void*);
void sort_dist (int,double*,int*);
response_unit* guess_a_location_01(int,int,point*); // Returns the first p
response_unit* guess_a_location_02(int,int,point*); // Returns p random with replace
response_unit* guess_a_location_03(int,int,point*); // Returns p random
void SQM_update_mst(mpf_t *mst,int m,int p,double Mu_NT,double **Dist,response_unit* X,mpf_t **f);
void SQM_expected_travel_time(mpf_t,int,int,double**,response_unit*,mpf_t**);
void SQM_mean_queue_delay(mpf_t,int,int,mpf_t*,mpf_t*,mpf_t**,mpf_t**);
/* Populate matrix of distances */
double** SQM_dist_matrix(SQM_instance*);

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
  Dist = SQM_dist_matrix(I);

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

    /* Update matrix of pfreferred servers */
    /* Debug cout << "Step 0.2" << endl; /* */
    for (int k = 0;k < m;k++) {
      for (int i = 0;i < p;i++)
	d[i] = Dist[k][X[i].location];
      sort_dist(p,d,a[k]);
    }

    /* Print Current solution */
    for (int i = 0;i < p;i++)
      cout << X[i].location << " ";
    cout << endl;

    mpf_set(T_r,t_r);

    /* **Step 0**
       Initialize Mean Service Time */
    /* Debug cout << "Step 0" << endl; /* */
    for (int i = 0;i < p;i++)
      mpf_set_d(mst[i],1 / Mu_NT);

    do {

      for (int i = 0;i < p;i++)
	mpf_set(MST[i],mst[i]);

      /* Update matrix of response times */
      /* Debug cout << "Step 0.1" << endl; /* */
      for (int i = 0;i < p;i++) {
	for (int k = 0;k < m;k++) {
	  mpf_set_d(Tao[i][k],(X[i].beta / X[i].v) * Dist[k][X[i].location]);
	  mpf_div_ui(Tao[i][k],Tao[i][k],60*24);
	  mpf_add(Tao[i][k],Tao[i][k],MST[i]);
	}
      }

      /* **Step 1**:
	Run the Hypercube Model */
      /* Debug cout << "Step 1.1" << endl; /* */
      jarvis_hypercube_approximation(m,p,Lambda,Tao,a,f);

      /* Step 2 
	 Update mean service time */
      /* Debug cout << "Step 2" << endl; /* */
      SQM_update_mst(mst,m,p,Mu_NT,Dist,X,f);
      
      /* Step 3 */
      /* Debug cout << "Step 3" << endl; /* */
      mpf_set_ui(delta_mu,0);
      for (int i = 0;i < p;i++) {
	mpf_sub(tmp,mst[i],MST[i]);
	mpf_abs(tmp,tmp);
	if (mpf_cmp(tmp,delta_mu) > 0)
	  mpf_set(delta_mu,tmp);
      }

      cout << "Delta in mst: " << mpf_get_d(delta_mu) << endl;
    } while (mpf_cmp_d(delta_mu,epsilon) > 0);
    
    /* *Expected Response Time* */
    mpf_set_ui(t_r,0);
    /* + expected travel time component */
    SQM_expected_travel_time(t_r,m,p,Dist,X,f);
    /* + mean queue delay component */
    SQM_mean_queue_delay(t_r,m,p,Lambda,mst,Tao,f);

    double *h_i = new double[m];
    num h;
    mpf_init(h);
    for (int i = 0;i < p;i++) {
      // Block A
      //cout << "\t\tBlock A " << i << endl;
      mpf_set_ui(h,0);
      for (int k = 0;k < m;k++) 
	mpf_add(h,h,f[i][k]);
      /*
      if (mpf_cmp_ui(h,0) == 0)
	cout << "for i = " << i+1 << " sum over f_ij is 0" << endl;
      */
      for (int k = 0;k < m;k++) {
	mpf_div(tmp,f[i][k],h);
	h_i[k] = mpf_get_d(tmp);
      }

      // Block B
      /* Solve te 1-median location model with h_i^j */
      X[i].location = Solve_1_median_location_model(m,n,Dist,h_i);
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

double SQM_response_time
(SQM_instance *I,
 int p, // Number of adjusters
 response_unit* X,
 double lambda, // mean rate per unit of time within service calls are generated in Poisson manner
 double Mu_NT // mean of non-travel time component of the service time
 ) {
  
  /* Variable definitions */
  double T_r;
  num *MST,*mst; // mean service time
  num t_r; // expected response time
  double **Dist; // Matrix of distances
  num **f,**Tao;
  double *d;
  double mu;
  /* */
  num tmp;
  num delta_mu;
  num *Lambda;
  int **a;
  int m = I->M; /* Number of demand points */
  int n = I->N; /* Number of potencial sites to locate a server*/
  point *V = I->V,*W = I->W;
 
  /* Populate matrix of distances */
  Dist = SQM_dist_matrix(I);

  /* Populate matrix of pfreferred servers */
  a = new int*[m]; 
  for (int k = 0;k < m;k++)
    a[k] = new int[p];

  for (int k = 0;k < m;k++) {
    for (int i = 0;i < p;i++)
      d[i] = Dist[k][X[i].location];
    sort_dist(p,d,a[k]);
  }

  /* Allocate memory for arrays and matrix */
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

  /* Allocate memory for mpf numbers */
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

  double demand = 0.0;
  for (int k = 0;k < m;k++) demand += I->V[k].demand;
  for (int k = 0;k < m;k++)
    mpf_set_d(Lambda[k],I->V[k].demand * lambda / demand);

  /* SERVICE MEAN TIME CALIBRATION */

  /* **Step 0**
     Initialize Mean Service Time */
  for (int i = 0;i < p;i++)
    mpf_set_d(mst[i],1 / Mu_NT);

  do {
    for (int i = 0;i < p;i++)
      mpf_set(MST[i],mst[i]);

    /* Update matrix of response times */
    for (int i = 0;i < p;i++) {
      for (int k = 0;k < m;k++) {
	//mpf_set_d(Tao[i][k],1 / Mu_NT + (X[i].beta / X[i].v) * Dist[k][X[i].location]/(60*24));
	mpf_set_d(Tao[i][k],(X[i].beta / X[i].v) * Dist[k][X[i].location]/(60*24));
	mpf_add(Tao[i][k],Tao[i][k],MST[i]);
      }
    }
      
    /* **Step 1**:
       Run the Hypercube Model */
    jarvis_hypercube_approximation(m,p,Lambda,Tao,a,f);

      
    /* Step 2 
       Update mean service time */
    SQM_update_mst(mst,m,p,Mu_NT,Dist,X,f);

    /* Step 3 */
    mpf_set_ui(delta_mu,0);
    for (int i = 0;i < p;i++) {
      mpf_sub(tmp,mst[i],MST[i]);
      mpf_abs(tmp,tmp);
      if (mpf_cmp(tmp,delta_mu) > 0)
	mpf_set(delta_mu,tmp);
    }
    
    /* Debug */ cout << "Delta in mst: " << mpf_get_d(delta_mu) << endl;
  } while (mpf_cmp_d(delta_mu,epsilon) > 0);

  /* *Expected Response Time* */
  mpf_init(t_r);
  /* + expected travel time component */
  SQM_expected_travel_time(t_r,m,p,Dist,X,f);
  /* + mean queue delay component */
  SQM_mean_queue_delay(t_r,m,p,Lambda,mst,Tao,f);
  /* Conver mpf to double */
  T_r = mpf_get_d(t_r);

  /* Deallocate memory for mpf numbers */
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
  for (int i = 0;i < p;i++)
    mpf_clear(mst[i]);
  for (int i = 0;i < p;i++)
    mpf_clear(MST[i]);
  mpf_clear(delta_mu);
  mpf_clear(tmp);

  /* Deallocate memory for arrays and matrix */
  for (int i = 0;i < p;i++) delete [] f[i];
  delete [] f;
  for (int i = 0;i < p;i++)
    delete [] Tao[i];
  delete [] Tao;
  delete [] Lambda;
  delete [] d;
  delete [] mst;
  delete [] MST;
  for (int k = 0;k < n;k++) delete a[k];
  delete [] a;
  for (int j=0;j < m;j++) delete [] Dist[j];
  delete [] Dist;
  
  return T_r;
}

void SQM_update_mst(mpf_t *mst,int m,int p,double Mu_NT,double **Dist,response_unit* X,mpf_t **f) {
  mpf_t h,tmp;
  mpf_init(h);
  mpf_init(tmp);
  for (int i = 0;i < p;i++) {

    mpf_set_ui(mst[i],0);
    for (int k = 0;k < m;k++) {
      mpf_set_d(tmp,1 / Mu_NT + (X[i].beta / X[i].v) * Dist[k][X[i].location]/(60*24));
      mpf_mul(tmp,tmp,f[i][k]);
      mpf_add(mst[i],mst[i],tmp);
     }

    mpf_set_ui(h,0);
    for (int k = 0;k < m;k++)
      mpf_add(h,h,f[i][k]);
    mpf_div(mst[i],mst[i],h);
  }
  mpf_clear(tmp);
  mpf_clear(h);
}


/* *Expected Travel Time Component* */
void SQM_expected_travel_time
(mpf_t t_r, /* stores the response time */
 int m, /* number of clients */
 int p, /* number of servers */
 double **Dist,
 response_unit* X,
 mpf_t **f
 ){
  mpf_t tmp;

  mpf_init(tmp);

  for (int i = 0;i < p;i++) {
    for (int k = 0;k < m;k++) {
      mpf_set_d(tmp,Dist[k][X[i].location]/X[i].v);
      mpf_div_ui(tmp,tmp,60*24);
      mpf_mul(tmp,tmp,f[i][k]);
      mpf_add(t_r,t_r,tmp);
    }
  } // m,p,Dist,X,f
}

/* mean queue delay component */
void SQM_mean_queue_delay(mpf_t t_r,int m,int p,mpf_t *Lambda,mpf_t *MST,mpf_t **Tao,mpf_t **f) {
  mpf_t mu;
  mpf_t P_B0,rho_i;
  mpf_t tmp;

  mpf_init(tmp);
  mpf_init_set_ui(mu,0);
  for (int i = 0;i < p;i++) {
    mpf_ui_div(tmp,1,MST[i]);
    mpf_add(mu,mu,tmp);
  }
  mpf_set_ui(tmp,0);
  for (int k = 0;k < m;k++)
    mpf_add(tmp,tmp,Lambda[k]); // lambda
  mpf_sub(tmp,mu,tmp); // mu - lambda
  mpf_pow_ui(tmp,tmp,2);// (mu - lambda)^2
  mpf_div(tmp,mu,tmp); // mu / (mu - lambda)^2

  mpf_init(rho_i);
  mpf_init_set_ui(P_B0,1);
  for (int i = 0;i < p;i++) {
    mpf_set_ui(rho_i,0);
    for (int k = 0;k < m;k++) {
      mpf_mul(tmp,Lambda[k],f[i][k]);
      mpf_mul(tmp,tmp,Tao[i][k]);
      mpf_add(rho_i,rho_i,tmp);
    }
    mpf_ui_sub(tmp,1,rho_i);
    mpf_mul(P_B0,P_B0,tmp);
  }
  mpf_mul(tmp,tmp,P_B0); // P_B0 * mu / (mu - lambda)^2
  mpf_add(t_r,t_r,tmp);

  mpf_clear(P_B0);
  mpf_clear(mu);
  mpf_clear(tmp);
}

double** SQM_dist_matrix(SQM_instance *I) {
  int n = I->N,m = I->M;
  double **Dist;

  Dist = new  double*[m];
  for (int j = 0;j < m;j++)
    Dist[j] = new double [n];

  for (int j = 0;j < m;j++) {
    for (int i = 0;i < n;i++) {
      Dist[j][i] = dist(&(I->V[j]),&(I->W[i]));
    }
  }
  return Dist;
}

int Solve_1_median_location_model(int m,int n,double **Dist,double *h) {
  int best_location = -1;
  double best_sol,sol;
  for (int j = 0;j < n;j++) {
    sol = 0.0;
    for (int k = 0;k < m;k++) 
      sol += h[k] * Dist[k][j];
    if (best_location == -1 || sol < best_sol) {
      best_location = j;
      best_sol = sol;
    }
  }
  return best_location;
}
