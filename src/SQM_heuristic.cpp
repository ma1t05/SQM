
#include "SQM.h"
#include "SQM_heuristic.h"
#include "mp_jarvis.h"
#include "MST.h"
#include "gnuplot.h"
#include "log.h"
#include <sstream>

/* Populate matrix of distances */
int Solve_1_median_location_model(int,int,double**,double*);
void SQM_improve_locations(SQM_solution*,mpf_t**);
void SQM_return_previous_solution(SQM_solution*);

void SQM_heuristic
(SQM_solution *Sol,
 double lambda, // mean rate per unit of time within service calls are generated in Poisson manner
 double Mu_NT // mean of non-travel time component of the service time
 ) {
  
  SQM_instance *I = Sol->get_instance();
  int p = Sol->get_servers(); // Number of adjusters
  /* Variable definitions */
  int m = I->demand_points(); /* Number of demand points */
  int n = I->potential_sites(); /* Number of potencial sites to locate a server*/
  int it = 0; /* iterator counter */
  mpf_t tmp;
  mpf_t delta_mu;
  mpf_t T_r,t_r; // expected response time
  double *d;
  mpf_t *MST,*mst; // mean service time
  double **F; // Version double of matrix f
  int key = rand(); // key number to naming plots
  char plot_output[32];
  mpf_t **f,**Tao;
  /* */
  mpf_t *Lambda;
 
  logDebug(cout << "Start Berman Heuristic" << "\r");
  /* Populate matrix of distances */
  F = new double*[p];
  for (int i = 0;i < p;i++)
    F[i] = new double[m];
  MST = new mpf_t[p];
  mst = new mpf_t[p];
  d = new double[p];
  Lambda = new mpf_t[m];
  Tao = new mpf_t*[p];
  for (int i = 0;i < p;i++)
    Tao[i] = new mpf_t[m];
  f = new mpf_t*[p];
  for (int i = 0;i < p;i++)
    f[i] = new mpf_t[m];

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

  double demand = I->total_demand();
  for (int k = 0;k < m;k++)
    mpf_set_d(Lambda[k],I->demand(k)->demand * lambda / demand);


  for (int i = 0;i < p;i++)
    LogFile << Sol->get_server_location(i) << " ";
  LogFile << endl;

  do {

    /* Update matrix of pfreferred servers */
    Sol->update_preferred_servers();

    /* Print Current solution */
    /*
    for (int i = 0;i < p;i++)
      cout << Sol->get_server_location(i) << " ";
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
	  mpf_set_d(Tao[i][k],Sol->get_server_rate(i) * I->distance(Sol->get_server_location(i),k));
	  mpf_div_ui(Tao[i][k],Tao[i][k],MINS_PER_BLOCK*BLOCKS_PER_HORIZON);
	  mpf_add(Tao[i][k],Tao[i][k],MST[i]);
	}
      }

      /* **Step 1**:
	Run the Hypercube Model */
      jarvis_hypercube_approximation(m,p,Lambda,Tao,Sol->preferred_servers(),f);

      /* Step 2 
	 Update mean service time */
      MST_update_mst(mst,Sol,Mu_NT,f);
      
      /* Step 3 */
      mpf_set_ui(delta_mu,0);
      for (int i = 0;i < p;i++) {
	mpf_sub(tmp,mst[i],MST[i]);
	mpf_abs(tmp,tmp);
	if (mpf_cmp(tmp,delta_mu) > 0)
	  mpf_set(delta_mu,tmp);
      }

      logDebug(cout << "Delta in mst: " << mpf_get_d(delta_mu) << endl);
    } while (mpf_cmp_d(delta_mu,epsilon) > 0);
    
    /* *Expected Response Time* */
    mpf_set_ui(t_r,0);
    /* + expected travel time component */
    MST_expected_travel_time(t_r,Sol,f);
    /* + mean queue delay component */
    MST_mean_queue_delay(t_r,m,p,Lambda,mst,Tao,f);

    mpf_sub(tmp,T_r,t_r);
    logDebug(cout << "Berman Heuristic ERT [" << ++it << "]:\t" << mpf_get_d(t_r) << (mpf_cmp_ui(tmp,0) < 0 ? "*\t" : "\t"));
    /* Print current solution */
    if (LogDebug) {
      for (int i = 0;i < p;i++)
	cout << Sol->get_server_location(i) << (Sol->get_server_past_location(i) != Sol->get_server_location(i) ? "*\t" : "\t");
      cout << endl;
    }
    /* mpf_abs(tmp,tmp); */ // This is a bad 

    if (mpf_cmp_d(tmp,epsilon) > 0) {
      /* Plot Iteration */
      if (LogDebug) {
	for (int i = 0;i < p;i++)
	  for (int k = 0;k < m;k++)
	    F[i][k] = mpf_get_d(f[i][k]);
	stringstream Key;
	sprintf(plot_output,"./plots/SQM_Heuristic_%010d",key);
	Key << it;
	plot_solution_allocation(Sol,F,string(plot_output),Key.str());
      }
      SQM_improve_locations(Sol,f);
    }
    else if (mpf_cmp_ui(tmp,0) < 0)
      SQM_return_previous_solution(Sol);

  } while (mpf_cmp_d(tmp,epsilon) > 0);
  //cout << endl;
  
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
  for (int i = 0;i < p;i++) delete [] F[i];
  delete [] F;
  for (int i = 0;i < p;i++) delete [] f[i];
  delete [] f;

  if (LogDebug) {
    for (int k = 0;k < n;k++)
      for (int i = 0;i < p;i++)
	if (Sol->get_server_location(i) == k) cout << k << " ";
    cout << endl;
    cout << "Finish Berman Heuristic" << endl;
  }
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

void SQM_improve_locations(SQM_solution *X,mpf_t **f) {
  SQM_instance *I = X->get_instance();
  int m = I->demand_points();
  int n = I->potential_sites();
  int p = X->get_servers();
  double **Dist;

  int new_location;
  double *h_i;
  mpf_t h,tmp;
  mpf_init(h);
  mpf_init(tmp);
  h_i = new double[m];
  for (int i = 0;i < p;i++) {
    // Block A
    logDebug(cout << "\t\tBlock A " << i << endl);
    mpf_set_ui(h,0);
    for (int k = 0;k < m;k++) 
      mpf_add(h,h,f[i][k]);
    if (LogInfo && mpf_cmp_ui(h,0) == 0)
      cout << "for i = " << i+1 << " sum over f_ij is 0" << endl;
    for (int k = 0;k < m;k++) {
      mpf_div(tmp,f[i][k],h);
      h_i[k] = mpf_get_d(tmp);
    }

    // Block B
    /* Solve te 1-median location model with h_i^j */
    new_location = Solve_1_median_location_model(m,n,Dist,h_i);
    X->set_server_location(i,new_location);
  }
  mpf_clear(tmp);
  mpf_clear(h);
  delete [] h_i;

  /* Print current solution to LogFile */
  for (int i = 0;i < p;i++) {
    LogFile << X->get_server_location(i);
    if (X->get_server_past_location(i) != X->get_server_location(i)) {
      LogFile << "*";
    }
    LogFile << "\t";
  }
  LogFile << endl;
}

void SQM_return_previous_solution(SQM_solution *X) {
  int p = X->get_servers();
  for (int i = 0;i < p;i++)
    X->set_server_location(i,X->get_server_past_location(i));
}
