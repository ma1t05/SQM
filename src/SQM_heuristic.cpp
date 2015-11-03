
#include "SQM_heuristic.h"
#include "instance-creator.h"
#include "mp_jarvis.h"
#include "MST.h"
#include "gnuplot.h"
#include <sstream>

/* Populate matrix of distances */
int Solve_1_median_location_model(int,int,double**,double*);
void SQM_improve_locations(response_unit*,int,int,int,double**,num**);
void SQM_return_previous_solution(response_unit*,int);

void SQM_heuristic
(SQM_instance *I,
 int p, // Number of adjusters
 double lambda, // mean rate per unit of time within service calls are generated in Poisson manner
 double Mu_NT, // mean of non-travel time component of the service time
 response_unit *X
 ) {
  
  /* Variable definitions */
  int m = I->M; /* Number of demand points */
  int n = I->N; /* Number of potencial sites to locate a server*/
  int it = 0; /* iterator counter */
  num tmp;
  num delta_mu;
  num T_r,t_r; // expected response time
  double *d;
  num *MST,*mst; // mean service time
  double **Dist; // Matrix of distances
  double **F; // Version double of matrix f
  int key = rand(); // key number to naming plots
  char plot_output[32];
  num **f,**Tao;
  /* */
  num *Lambda;
  int **a;
  point *V = I->V,*W = I->W;
 
  /* Debug cout << "Start Berman Heuristic" << "\r"; /* */
  /* Populate matrix of distances */
  Dist = SQM_dist_matrix(I);
  F = new double*[p];
  for (int i = 0;i < p;i++)
    F[i] = new double[m];
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

  do {

    /* Update matrix of pfreferred servers */
    for (int k = 0;k < m;k++) {
      for (int i = 0;i < p;i++)
	d[i] = Dist[k][X[i].location];
      sort_dist(p,d,a[k]);
    }

    /* Print Current solution */
    /*
    for (int i = 0;i < p;i++)
      cout << X[i].location << " ";
    cout << endl;
    */

    mpf_set(T_r,t_r);

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
	  mpf_set_d(Tao[i][k],(X[i].beta / X[i].v) * Dist[k][X[i].location]);
	  mpf_div_ui(Tao[i][k],Tao[i][k],MINS_PER_BLOCK*BLOCKS_PER_HORIZON);
	  mpf_add(Tao[i][k],Tao[i][k],MST[i]);
	}
      }

      /* **Step 1**:
	Run the Hypercube Model */
      jarvis_hypercube_approximation(m,p,Lambda,Tao,a,f);

      /* Plot Iteration */
      for (int i = 0;i < p;i++)
	for (int k = 0;k < m;k++)
	  F[i][k] = mpf_get_d(f[i][k]);
      stringstream Key;
      sprintf(plot_output,"./plots/SQM_Heuristic_%010d",key);
      Key << ++it;
      plot_solution_allocation(I,p,X,F,string(plot_output),Key.str());

      /* Step 2 
	 Update mean service time */
      MST_update_mst(mst,m,p,Mu_NT,Dist,X,f);
      
      /* Step 3 */
      mpf_set_ui(delta_mu,0);
      for (int i = 0;i < p;i++) {
	mpf_sub(tmp,mst[i],MST[i]);
	mpf_abs(tmp,tmp);
	if (mpf_cmp(tmp,delta_mu) > 0)
	  mpf_set(delta_mu,tmp);
      }

      /* Debug cout << "Delta in mst: " << mpf_get_d(delta_mu) << endl; /* */
    } while (mpf_cmp_d(delta_mu,epsilon) > 0);
    
    /* *Expected Response Time* */
    mpf_set_ui(t_r,0);
    /* + expected travel time component */
    MST_expected_travel_time(t_r,m,p,Dist,X,f);
    /* + mean queue delay component */
    MST_mean_queue_delay(t_r,m,p,Lambda,mst,Tao,f);

    /*
    mpf_sub(tmp,T_r,t_r);
    cout << "Berman Heuristic ERT [" << ++it << "]:\t" << mpf_get_d(t_r);
    if (mpf_cmp_ui(tmp,0) < 0)
      cout << "*\t";
    else cout << "\t";
    /* Print current solution 
    for (int i = 0;i < p;i++) {
      cout << X[i].location;
      if (X[i].past_location != X[i].location) {
	cout << "*";
      }
      cout << "\t";
    }
    cout << endl;
    
    /* mpf_abs(tmp,tmp); */ // This is a bad 

    if (mpf_cmp_d(tmp,epsilon) > 0)
      SQM_improve_locations(X,m,n,p,Dist,f);
    else if (mpf_cmp_ui(tmp,0) < 0)
      SQM_return_previous_solution(X,p);

  } while (mpf_cmp_d(tmp,epsilon) > 0);
  //cout << endl;
  
  for (int k = 0;k < n;k++) delete a[k];
  delete [] a;

  mpf_clear(T_r);
  mpf_clear(t_r);
  for (int i = 0;i < p;i++) {
    for (int k = 0;k < k;k++)
      mpf_clear(f[i][k]);
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
  for (int i = 0;i < p;i++) delete [] F[i];
  delete [] F;
  for (int i = 0;i < p;i++) delete [] f[i];
  delete [] f;

  /*
  for (int k = 0;k < n;k++)
    for (int i = 0;i < p;i++)
      if (X[i].location == k) cout << k << " ";
  cout << endl;
  */
  /* Debug cout << "Finish Berman Heuristic" << endl; /* */
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

void SQM_improve_locations(response_unit *X,int m,int n,int p,double **Dist,num **f) {
  int new_location;
  double *h_i;
  num h,tmp;
  mpf_init(h);
  mpf_init(tmp);
  h_i = new double[m];
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
    new_location = Solve_1_median_location_model(m,n,Dist,h_i);
    X[i].past_location = X[i].location;
    X[i].location = new_location;
  }
  mpf_clear(tmp);
  mpf_clear(h);
  delete [] h_i;

  /* Print current solution to LogFile */
  for (int i = 0;i < p;i++) {
    LogFile << X[i].location;
    if (X[i].past_location != X[i].location) {
      LogFile << "*";
    }
    LogFile << "\t";
  }
  LogFile << endl;
}

void SQM_return_previous_solution(response_unit *X,int p) {
  for (int i = 0;i < p;i++)
    X[i].location = X[i].past_location;
}
