
#include "SQM.h"
#include "MST.h"
#include "mp_jarvis.h"
#include "log.h"

double MST_response_time
(SQM_solution *X,
 double lambda, // mean rate per unit of time within service calls are generated in Poisson manner
 double Mu_NT // mean of non-travel time component of the service time
 ) {
  
  /* Variable definitions */
  double T_r;
  mpf_t *MST,*mst; // mean service time
  mpf_t t_r; // expected response time
  mpf_t **f,**Tao;
  double *d;
  double mu;
  /* */
  mpf_t tmp;
  mpf_t delta_mu;
  mpf_t *Lambda;
  int **a;
  SQM_instance *I = X->get_instance();
  int m = I->demand_points(); /* Number of demand points */
  int n = I->potential_sites(); /* Number of potencial sites to locate a server*/ 
  int p = X->get_servers ();
  point *V = I->demand(0),*W = I->site(0);
 
  /* Populate matrix of pfreferred servers */
  X->update_preferred_servers();
  a = X->preferred_servers();

  /* Allocate memory for arrays and matrix */
  MST = new mpf_t[p];
  mst = new mpf_t[p];
  Lambda = new mpf_t[m];
  Tao = new mpf_t*[p];
  for (int i = 0;i < p;i++)
    Tao[i] = new mpf_t[m];
  f = new mpf_t*[p];
  for (int i = 0;i < p;i++)
    f[i] = new mpf_t[m];

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

  double demand = I->total_demand();
  for (int k = 0;k < m;k++)
    mpf_set_d(Lambda[k],I->demand(k)->demand * lambda / demand);

  /* SERVICE MEAN TIME CALIBRATION */
  /* Debug cout << "Solution:"; 
  for (int i = 0;i < p;i++)
    cout << " " << X[i].location;
  cout << endl; /* */

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
	mpf_set_d(Tao[i][k],X->get_server_rate(i) * I->distance(X->get_server_location(i),k));
	mpf_div_ui(Tao[i][k],Tao[i][k],MINS_PER_BLOCK*BLOCKS_PER_HORIZON);
	mpf_add(Tao[i][k],Tao[i][k],MST[i]);
      }
    }
      
    /* **Step 1**
	 Run the Hypercube Model */
    jarvis_hypercube_approximation(m,p,Lambda,Tao,a,f);

    /* 
    mpf_t sum;
    mpf_init(sum);
    cout << "Sum_k[]:";
    for (int i = 0;i < p;i++) {
      mpf_set_ui(sum,0);
      for (int k = 0;k < m;k++)
	mpf_add(sum,sum,f[i][k]);
      cout << " " << mpf_get_d(sum);
    }
    cout << endl;
    if (mpf_cmp_ui(sum,0) == 0) {
      for (int i = 0;i < p;i++) {
	cout << "[" << i+1 << "]:";
	for (int k = 0;k < m;k++)
	  cout << " " << mpf_get_d(Tao[i][k]);
	cout << endl;
      }
    }
    mpf_clear(sum);
    */

    /* Step 2
       Update mean service time */
    MST_update_mst(mst,X,Mu_NT,f);

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
  mpf_init(t_r);
  /* + expected travel time component */
  MST_expected_travel_time(t_r,X,f);
  /* + mean queue delay component */
  MST_mean_queue_delay(t_r,m,p,Lambda,mst,Tao,f);
  /* Conver mpf to double */
  T_r = mpf_get_d(t_r);

  /* Deallocate memory for mpf numbers */
  mpf_clear(t_r);
  for (int i = 0;i < p;i++) {
    for (int k = 0;k < m;k++)
      mpf_clear(f[i][k]);
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
  delete [] mst;
  delete [] MST;
  //cout << "Finish Response time" << endl << endl;
  return T_r;
}

void MST_update_mst(mpf_t *mst,SQM_solution *X,double Mu_NT,mpf_t **f) {
  SQM_instance *I = X->get_instance ();
  int m = I->demand_points();
  int p = X->get_servers();
  mpf_t h,tmp;
  mpf_init(h);
  mpf_init(tmp);
  /* cout << "Start update mst" << endl; */
  for (int i = 0;i < p;i++) {

    mpf_set_ui(mst[i],0);
    for (int k = 0;k < m;k++) {
      mpf_set_d(tmp,1 / Mu_NT + X->get_server_rate(i) * I->distance(X->get_server_location(i),k)/(MINS_PER_BLOCK*BLOCKS_PER_HORIZON));
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
      mpf_set_d(tmp,I->distance(X->get_server_location(i),k)/X->get_server_speed(i));
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
  logDebug(cout << "mean queue delay FINISH" << endl);
}
