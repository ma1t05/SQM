
#include "SQM.h"
#include "MST.h"
#include "mp_jarvis.h"
#include "log.h"

void _MST_mpf_init(mpf_t**,int);
void _MST_mpf_init(mpf_t***,int,int);
void _MST_mpf_clear(mpf_t**,int);
void _MST_mpf_clear(mpf_t***,int,int);

double MST_response_time (SQM_solution *X) {
  
  /* Variable definitions */
  double T_r;
  mpf_t t_r; // expected response time
  mpf_t *mst; // mean service time
  mpf_t **f,**Tao;
  double *d;
  double mu;
  /* */
  mpf_t *Lambda;
  SQM_instance *I = X->get_instance();
  int m = I->demand_points(); /* Number of demand points */
  int p = X->get_servers ();
 
  logDebug(cout << "Start: MST_response_time" << endl);
  logDebug(cout << "MST_response_time: Populate matrix of pfreferred servers" << endl);
  X->update_preferred_servers();

  logDebug(cout << "MST_response_time: Allocate memory for mpf numbers" << endl);
  _MST_mpf_init(&mst,p);
  _MST_mpf_init(&Lambda,m);
  _MST_mpf_init(&Tao,p,m);
  _MST_mpf_init(&f,p,m);

  logDebug(cout <<"MST_response_time: Run Service Mean Time Calibration" << endl);
  MST_Calibration(f,mst,Tao,Lambda,X);

  /* *Expected Response Time* */
  logDebug(cout <<"MST_response_time: Calcule Expected Response Time" << endl);
  mpf_init(t_r);
  /* + expected travel time component */
  logDebug(cout <<"MST_response_time: add travel time component" << endl);
  MST_expected_travel_time(t_r,X,f);
  /* + mean queue delay component */
  logDebug(cout <<"MST_response_time: add mean queue delay component" << endl);
  MST_mean_queue_delay(t_r,m,p,Lambda,mst,Tao,f);
  /* Conver mpf to double */
  logDebug(cout <<"MST_response_time: convert mpf_t to double" << endl);
  T_r = mpf_get_d(t_r);

  /* Deallocate memory for mpf numbers */
  mpf_clear(t_r);
  _MST_mpf_clear(&f,p,m);
  _MST_mpf_clear(&Tao,p,m);
  _MST_mpf_clear(&Lambda,m);
  _MST_mpf_clear(&mst,p);

  logDebug(cout << "Finish: MST_response_time" << endl);
  return T_r;
}

void MST_update_mst(mpf_t *mst,SQM_solution *X,mpf_t **f) {
  SQM_instance *I = X->get_instance ();
  int m = I->demand_points();
  int p = X->get_servers();
  double Mu_NT = X->get_non_travel_time();
  mpf_t h,tmp;
  mpf_init(h);
  mpf_init(tmp);
  /* cout << "Start update mst" << endl; */
  for (int i = 0;i < p;i++) {

    mpf_set_ui(mst[i],0);
    for (int k = 0;k < m;k++) {
      mpf_set_d(tmp,1 / Mu_NT + X->get_server_rate(i) * X->distance(i,k)/(MINS_PER_BLOCK*BLOCKS_PER_HORIZON));
      mpf_mul(tmp,tmp,f[i][k]);
      mpf_add(mst[i],mst[i],tmp);
    }

    mpf_set_ui(h,0);
    for (int k = 0;k < m;k++) {
      mpf_add(h,h,f[i][k]);
    }
    /*
    cout << "[" << i << "] h" << X[i].location 
	 << " != 0?" << mpf_cmp_ui(h,0) << "\r";
    */
    if (mpf_cmp_ui(h,0) > 0) /* WARNING: This sum shoudn't be 0 */
      mpf_div(mst[i],mst[i],h);
    else cout << endl;
  }
  mpf_clear(tmp);
  mpf_clear(h);
  /* cout << "End update mst" << endl; */
}


/* *Expected Travel Time Component* */
void MST_expected_travel_time
(mpf_t t_r, /* stores the response time */
 SQM_solution* X,
 mpf_t **f
 ){
  SQM_instance *I = X->get_instance();
  int m = I->demand_points();
  int p = X->get_servers();
  mpf_t tmp;
  mpf_init(tmp);
  for (int i = 0;i < p;i++) {
    for (int k = 0;k < m;k++) {
      mpf_set_d(tmp,X->distance(i,k)/X->get_server_speed(i));
      mpf_div_ui(tmp,tmp,MINS_PER_BLOCK*BLOCKS_PER_HORIZON);
      mpf_mul(tmp,tmp,f[i][k]);
      mpf_add(t_r,t_r,tmp);
    }
  } // m,p,Dist,X,f
  /* Debug cout << "expected travel time: " << MINS_PER_BLOCK * BLOCKS_PER_HORIZON * mpf_get_d(t_r) << endl; /* */
}

/* mean queue delay component */
void MST_mean_queue_delay(mpf_t t_r,int m,int p,mpf_t *Lambda,mpf_t *MST,mpf_t **Tao,mpf_t **f) {
  mpf_t mu;
  mpf_t P_B0,rho_i;
  mpf_t tmp;

  logDebug(cout << "mean queue delay START" << endl);

  mpf_init(tmp);
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

  mpf_mul(tmp,tmp,P_B0); // P_B0 * mu / (mu - lambda)^2
  mpf_add(t_r,t_r,tmp);

  mpf_clear(P_B0);
  mpf_clear(mu);
  mpf_clear(tmp);
  logDebug(cout << "mean queue delay FINISH" << endl);
}

double* MST_workload(SQM_solution *Sol) {
  
  /* Variable definitions */
  mpf_t tmp,rho_i;
  mpf_t *mst; // mean service time
  mpf_t **f,**Tao; 
  mpf_t *Lambda;
  double *rho;
  double mu;
  /* */
  SQM_instance *I = Sol->get_instance();
  int m = I->demand_points(); /* Number of demand points */
  int p = Sol->get_servers ();
 
  logDebug(cout << "Start: MST_workload" << endl);
  logDebug(cout << "MST_workload: Populate matrix of pfreferred servers" << endl);
  Sol->update_preferred_servers();

  rho = new double[p];
  logDebug(cout << "MST_workload: Allocate memory for mpf numbers" << endl);
  _MST_mpf_init(&f,p,m);
  _MST_mpf_init(&mst,p);
  _MST_mpf_init(&Tao,p,m);
  _MST_mpf_init(&Lambda,m);
  mpf_init(tmp);
  mpf_init(rho_i);

  logDebug(cout <<"MST_workload: Run Service Mean Time Calibration" << endl);
  MST_Calibration(f,mst,Tao,Lambda,Sol);

  for (int i = 0;i < p;i++) {
    mpf_set_ui(rho_i,0);
    for (int k = 0;k < m;k++) {
      mpf_mul(tmp,Lambda[k],f[i][k]);
      mpf_mul(tmp,tmp,Tao[i][k]);
      mpf_add(rho_i,rho_i,tmp);
    }
    rho[i] = mpf_get_d(rho_i);
  }

  /* Deallocate memory for mpf numbers */
  mpf_clear(rho_i);
  mpf_clear(tmp);
  _MST_mpf_clear(&Lambda,m);
  _MST_mpf_clear(&Tao,p,m);
  _MST_mpf_clear(&mst,p);
  _MST_mpf_clear(&f,p,m);

  logDebug(cout << "Finish: MST_workload" << endl);
  return rho;
}

void MST_Calibration(mpf_t **f,mpf_t *mst,mpf_t **Tao,mpf_t *Lambda,SQM_solution *X) {
  SQM_instance *I = X->get_instance();
  int m = I->demand_points();
  int p = X->get_servers();
  double lambda = X->get_arrival_rate(); // mean rate per unit of time within service calls are generated in Poisson manner
  double Mu_NT = X->get_non_travel_time();
  mpf_t *MST; // mean service time
  mpf_t delta_mu;
  mpf_t tmp;
  double demand;
  
  logDebug(cout << "Start: MST_Calibration" << endl);
  logDebug(cout << "/* SERVICE MEAN TIME CALIBRATION */" << endl);

  mpf_init(tmp);
  mpf_init(delta_mu);
  _MST_mpf_init(&MST,p);

  demand = I->total_demand();
  for (int k = 0;k < m;k++)
    mpf_set_d(Lambda[k],I->demand(k)->demand * lambda / demand);

  /* **Step 0**
     Initialize Mean Service Time */
  for (int i = 0;i < p;i++)
    mpf_set_d(mst[i],1 / Mu_NT);

  do {
    for (int i = 0;i < p;i++)
      mpf_set(MST[i],mst[i]);

    logDebug(cout << "\t**Step 0** - Update matrix of response times" << endl);
    for (int i = 0;i < p;i++) {
      for (int k = 0;k < m;k++) {
	mpf_set_d(Tao[i][k],X->get_server_rate(i) * X->distance(i,k));
	mpf_div_ui(Tao[i][k],Tao[i][k],MINS_PER_BLOCK*BLOCKS_PER_HORIZON);
	mpf_add(Tao[i][k],Tao[i][k],MST[i]);
      }
    }
      
    logDebug(cout << "\t**Step 1** - Run the Hypercube Model" << endl);
    jarvis_hypercube_approximation(m,p,Lambda,Tao,X->preferred_servers(),f);

    logDebug(cout << "\t**Step 2** - Update mean service time" << endl);
    MST_update_mst(mst,X,f);

    logDebug(cout << "\t**Step 3** - ");
    mpf_set_ui(delta_mu,0);
    for (int i = 0;i < p;i++) {
      mpf_sub(tmp,mst[i],MST[i]);
      mpf_abs(tmp,tmp);
      if (mpf_cmp(tmp,delta_mu) > 0)
	mpf_set(delta_mu,tmp);
    }
    
    logDebug(cout << "Delta in mst: " << mpf_get_d(delta_mu) << endl);
  } while (mpf_cmp_d(delta_mu,epsilon) > 0);

  _MST_mpf_clear(&MST,p);
  mpf_clear(delta_mu);
  mpf_clear(tmp);

}

 void _MST_mpf_init(mpf_t **a,int n) {
   mpf_t *b;

   b = new mpf_t[n];

  for (int i = 0;i < n;i++)
    mpf_init(b[i]);
  *a = b;
}

void _MST_mpf_init(mpf_t ***a,int m,int n) {
  mpf_t **b;

  b = new mpf_t*[m];
  for (int i = 0;i < m;i++)
    b[i] = new mpf_t[n];

  for (int i = 0;i < m;i++)
    for (int k = 0;k < n;k++)
      mpf_init(b[i][k]);
  *a = b;
}

void _MST_mpf_clear(mpf_t **a,int n) {
   mpf_t *b = *a;

  for (int k = 0;k < n;k++)
    mpf_clear(b[k]);

  delete [] b;
  *a = NULL;
}

void _MST_mpf_clear(mpf_t ***a,int m,int n) {
   mpf_t **b = *a;

  for (int i = 0;i < m;i++)
    for (int k = 0;k < n;k++)
      mpf_clear(b[i][k]);

  for (int i = 0;i < m;i++)
    delete [] b[i];
  delete [] b;
  *a = NULL;
}

