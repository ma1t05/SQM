
#include "SQM_heuristic.h"

/* Populate matrix of distances */
int Solve_1_median_location_model(int,int,double**,double*);
void SQM_improve_locations(SQM_solution*,mpf_t**);
void SQM_return_previous_solution(SQM_solution*);

void SQM_heuristic (SQM_solution *Sol) {
  /* Variable definitions */
  mpf_t tmp;
  mpf_t T_r,t_r; // expected response time
  mpf_t *mst; // mean service time
  mpf_t **f,**Tao;
  mpf_t *Lambda;
  /* */
  SQM_instance *I = Sol->get_instance();
  int p = Sol->get_servers(); // Number of adjusters
  int m = I->demand_points(); /* Number of demand points */
  int n = I->potential_sites(); /* Number of potencial sites to locate a server*/
  int it = 0; /* iterator counter */
  int key = rand(); // key number to naming plots
  char plot_output[32];
 
  logDebug(cout << "Start Berman Heuristic" << endl);
  _MST_mpf_init(&mst,p);
  _MST_mpf_init(&Lambda,m);
  _MST_mpf_init(&Tao,p,m);
  _MST_mpf_init(&f,p,m);
  mpf_init(tmp);
  mpf_init(T_r);
  mpf_init_set_d(t_r,99999.9);

  for (int i = 0;i < p;i++)
    LogFile << Sol->get_server_location(i) << "\t";
  LogFile << Sol->get_response_time() << endl;

  do {

    logDebug(cout << "/* Update matrix of pfreferred servers */" << endl);
    Sol->update_preferred_servers();

    /* Print Current solution */
    /*
    for (int i = 0;i < p;i++)
      cout << Sol->get_server_location(i) << " ";
    cout << endl;
    */

    mpf_set(T_r,t_r);
    MST_Calibration(f,mst,Tao,Lambda,Sol);
    
    logDebug(cout << "/* *Expected Response Time* */" << endl);
    mpf_set_ui(t_r,0);
    logDebug(cout << "/* + expected travel time component */" << endl);
    MST_expected_travel_time(t_r,Sol,f);
    logDebug(cout << "/* + mean queue delay component */" << endl);
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
      logDebug(cout << "/* Plot Iteration */" << endl);
      if (LogDebug) {

	double **F; // Version double of matrix f
	F = new double*[p];
	for (int i = 0;i < p;i++) {
	  F[i] = new double[m];
	  for (int k = 0;k < m;k++)
	    F[i][k] = mpf_get_d(f[i][k]);
	}

	stringstream Key;
	sprintf(plot_output,"./plots/SQM_Heuristic_%010d",key);
	Key << it;
	plot_solution_allocation(Sol,F,string(plot_output),Key.str());
	for (int i = 0;i < p;i++) delete [] F[i];
	delete [] F;
      }
      logDebug(cout << "Improve locations" << endl);
      SQM_improve_locations(Sol,f);
    }
    else if (mpf_cmp_ui(tmp,0) < 0)
      SQM_return_previous_solution(Sol);

  } while (mpf_cmp_d(tmp,epsilon) > 0);
  //cout << endl;
  
  mpf_clear(t_r);
  mpf_clear(T_r);
  mpf_clear(tmp);
  _MST_mpf_clear(&f,p,m);
  _MST_mpf_clear(&Tao,p,m);
  _MST_mpf_clear(&Lambda,m);
  _MST_mpf_clear(&mst,p);

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
    /* Here is the error Dist = NULL */
    new_location = Solve_1_median_location_model(m,n,I->get_distances_matrix(),h_i);
    X->set_server_location(i,new_location);
  }
  delete [] h_i;
  mpf_clear(tmp);
  mpf_clear(h);

  /* Print current solution to LogFile */
  for (int i = 0;i < p;i++) {
    LogFile << X->get_server_location(i);
    if (X->get_server_past_location(i) != X->get_server_location(i)) LogFile << "*";
    LogFile << "\t";
  }
  LogFile << X->get_response_time() << endl;
}

void SQM_return_previous_solution(SQM_solution *X) {
  int p = X->get_servers();
  for (int i = 0;i < p;i++)
    X->set_server_location(i,X->get_server_past_location(i));
}
