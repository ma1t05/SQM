
#include "MST.h"
#include <exception>

double MST_response_time (SQM_instance *Inst,int p,server *Servers,
			  int **preferred_servers) 
{
  /* Variable definitions */
  double T_r;
  mpf_t t_r; // expected response time
  mpf_t *mst; // mean service time
  mpf_t **f,**Tao;
  double *d;
  double mu;
  /* */
  mpf_t *Lambda;
  int m = Inst->demand_points(); /* Number of demand points */
  log_depth++;
  string tag = log_tag("MST_response_time: ");
 
  logDebug(cout << tag << "Start" << endl);

  logDebug(cout << tag << "Allocate memory for mpf numbers" << endl);
  _MST_mpf_init(&mst,p);
  _MST_mpf_init(&Lambda,m);
  _MST_mpf_init(&Tao,p,m);
  _MST_mpf_init(&f,p,m);

  logDebug(cout << tag << "Run Service Mean Time Calibration" << endl);
  MST_Calibration(f,mst,Tao,Lambda,Inst,p,Servers,preferred_servers);

  /* *Expected Response Time* */
  logDebug(cout << tag << "Calcule Expected Response Time" << endl);
  mpf_init(t_r);
  /* + expected travel time component */
  logDebug(cout << tag << "add travel time component" << endl);
  MST_expected_travel_time(t_r,Inst,p,Servers,f);
  /* + mean queue delay component */
  logDebug(cout << tag << "add mean queue delay component" << endl);
  MST_mean_queue_delay(t_r,m,p,Lambda,mst,Tao,f);
  /* Conver mpf to double */
  logDebug(cout << tag << "convert mpf_t to double" << endl);
  T_r = mpf_get_d(t_r);

  /* Deallocate memory for mpf numbers */
  mpf_clear(t_r);
  _MST_mpf_clear(&f,p,m);
  _MST_mpf_clear(&Tao,p,m);
  _MST_mpf_clear(&Lambda,m);
  _MST_mpf_clear(&mst,p);

  logDebug(cout << tag << "Finish" << endl);
  log_depth--;
  return T_r;
}

void MST_update_mst(mpf_t *mst,SQM_instance *Inst,int p,server *Servers,
		    mpf_t **f) 
{
  int m = Inst->demand_points();
  double Mu_NT = Inst->get_service_rate();
  double distance;
  mpf_t h,tmp;
  mpf_init(h);
  mpf_init(tmp);
  log_depth++;
  string tag = log_tag("MST_update_mst: ");
  /* cout << tag << "Start" << endl; */
  for (int i = 0;i < p;i++) {

    mpf_set_ui(mst[i],0);
    for (int k = 0;k < m;k++) {
      distance = Inst->distance(Servers[i].get_location(),k);
      mpf_set_d(tmp,1 / Mu_NT + Servers[i].get_rate() * distance);
      mpf_div_ui(tmp,tmp,MINS_PER_BLOCK * BLOCKS_PER_HORIZON);
      mpf_mul(tmp,tmp,f[i][k]);
      /* cout << " " << mpf_get_d(f[i][k]);*/
      mpf_add(mst[i],mst[i],tmp);
    }
    /*cout << endl;*/

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
    else logError(cout << tag << "\tWARNING! sum over j of f_" << i+1 
		  << "j = 0" << endl);
  }
  mpf_clear(tmp);
  mpf_clear(h);
  /* cout << "End update mst" << endl; */
  log_depth--;
}


/* *Expected Travel Time Component* */
void MST_expected_travel_time
(mpf_t t_r, /* stores the response time */
 SQM_instance *Inst,
 int p,
 server *Servers,
 mpf_t **f
 ){

  int m = Inst->demand_points();
  mpf_t tmp;
  mpf_init(tmp);
  for (int i = 0;i < p;i++) {
    for (int k = 0;k < m;k++) {
      mpf_set_d(tmp,Inst->distance(Servers[i].get_location(),k)/Servers[i].get_speed());
      mpf_div_ui(tmp,tmp,MINS_PER_BLOCK*BLOCKS_PER_HORIZON);
      mpf_mul(tmp,tmp,f[i][k]);
      mpf_add(t_r,t_r,tmp);
    }
  } // m,p,Dist,X,f
  /* Debug cout << "expected travel time: " << MINS_PER_BLOCK * BLOCKS_PER_HORIZON * mpf_get_d(t_r) << endl; /* */
}

/* mean queue delay component */
void MST_mean_queue_delay(mpf_t t_r,int m,int p,mpf_t *Lambda,mpf_t *MST,
			  mpf_t **Tao,mpf_t **f)
{
  mpf_t mu;
  mpf_t P_B0,rho_i;
  mpf_t tmp;
  log_depth++;
  string tag = log_tag("MST_mean_queue_delay: ");

  logDebug(cout << tag << "Start" << endl);
  mpf_init(tmp);
  mpf_init_set_ui(P_B0,1);
  mpf_init(rho_i);
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
  mpf_clear(rho_i);

  mpf_init_set_ui(mu,0);
  for (int i = 0;i < p;i++) {
    if (mpf_cmp_ui(MST[i],0) > 0){
      mpf_ui_div(tmp,1,MST[i]);
      mpf_add(mu,mu,tmp);
    }
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
  logDebug(cout << tag << "Finish" << endl);
  log_depth--;
}

double* MST_workload(SQM_instance *Inst,int p,server *Servers,
		     int **preferred_servers) 
{
  /* Variable definitions */
  mpf_t tmp,rho_i;
  mpf_t *mst; // mean service time
  mpf_t **f,**Tao; 
  mpf_t *Lambda;
  double *rho;
  double mu;
  /* */
  int m = Inst->demand_points(); /* Number of demand points */
  log_depth++;
  string tag = log_tag("MST_workload: ");

  logDebug(cout << tag << "Start" << endl);
  logDebug(cout << tag << "Populate matrix of pfreferred servers" << endl);

  rho = new double[p];
  logDebug(cout << tag << "Allocate memory for mpf numbers" << endl);
  _MST_mpf_init(&f,p,m);
  _MST_mpf_init(&mst,p);
  _MST_mpf_init(&Tao,p,m);
  _MST_mpf_init(&Lambda,m);
  mpf_init(tmp);
  mpf_init(rho_i);

  logDebug(cout << tag << "Run Service Mean Time Calibration" << endl);
  MST_Calibration(f,mst,Tao,Lambda,Inst,p,Servers,preferred_servers);

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

  logDebug(cout << tag << "Finish" << endl);
  log_depth--;
  return rho;
}

void MST_Calibration(mpf_t **f,mpf_t *mst,mpf_t **Tao,mpf_t *Lambda,
		     SQM_instance *Inst,int p,server *Servers,
		     int **preferred_servers) 
{
  int m = Inst->demand_points();
  double lambda = Inst->get_arrival_rate(); // mean rate per unit of time within service calls are generated in Poisson manner
  double Mu_NT = Inst->get_service_rate();
  mpf_t *MST; // mean service time
  mpf_t delta_mu;
  mpf_t tmp;
  double demand;
  double distance;
  log_depth++;
  string tag = log_tag("MST_Calibration: ");
  
  logDebug(cout << tag << "Start" << endl);
  logDebug(cout << tag << "/* SERVICE MEAN TIME CALIBRATION */" << endl);

  mpf_init(tmp);
  mpf_init(delta_mu);
  _MST_mpf_init(&MST,p);

  demand = Inst->total_demand();
  for (int k = 0;k < m;k++)
    mpf_set_d(Lambda[k],Inst->demand(k)->demand * lambda / demand);

  /* **Step 0**
     Initialize Mean Service Time */
  for (int i = 0;i < p;i++) {
    mpf_set_d(mst[i],1 / Mu_NT);
    mpf_div_ui(mst[i],mst[i],MINS_PER_BLOCK*BLOCKS_PER_HORIZON);
  }

  do {
    for (int i = 0;i < p;i++)
      mpf_set(MST[i],mst[i]);

    logDebug(cout << tag << " **Step 0** - Update matrix of response times" << endl);
    for (int k = 0;k < m;k++) {
      for (int i = 0;i < p;i++) {
	distance = Inst->distance(Servers[i].get_location(),k);
	mpf_set_d(Tao[i][k],Servers[i].get_rate() * distance);
	mpf_add(Tao[i][k],Tao[i][k],MST[i]);
	mpf_div_ui(Tao[i][k],Tao[i][k],MINS_PER_BLOCK*BLOCKS_PER_HORIZON);
	/* cout << "\t" << mpf_get_d(Tao[i][k]);*/
      }
      /* cout << endl;*/
    }
      
    logDebug(cout << tag << " **Step 1** - Run the Hypercube Model" << endl);
    jarvis_hypercube_approximation(m,p,Lambda,Tao,preferred_servers,f);

    logDebug(cout << tag << " **Step 2** - Update mean service time" << endl);
    MST_update_mst(mst,Inst,p,Servers,f);

    logDebug(cout << tag << " **Step 3** - ");
    mpf_set_ui(delta_mu,0);
    for (int i = 0;i < p;i++) {
      mpf_sub(tmp,mst[i],MST[i]);
      mpf_abs(tmp,tmp);
      if (mpf_cmp(tmp,delta_mu) > 0)
	mpf_set(delta_mu,tmp);
    }
    
    logDebug(cout << tag << " Delta in mst: " << mpf_get_d(delta_mu) << endl);
  } while (mpf_cmp_d(delta_mu,JARVIS_EPSILON) > 0);

  _MST_mpf_clear(&MST,p);
  mpf_clear(delta_mu);
  mpf_clear(tmp);
  logDebug(cout << tag << " Finish" << endl);
  log_depth--;
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

