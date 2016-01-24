
#include "SQM_GRASP.h"

#define min(A,B) ((A)>(B)?(A):(B))

bool GRASP_closest_to_b(SQM_instance *I,int node,int center_a,int center_b);
int GRASP_nearest_server(SQM_solution *Sol,int j);

SQM_solution* GRASP
(SQM_instance &Inst,
 int p, // Number of adjusters
 double v, // Speed
 double beta, // Beta
 double alpha // Random factor {1: random, 0: greedy}
 ) {
  int n = Inst.potential_sites(),m = Inst.demand_points();
  int r;
  int element;
  int *rcl;
  double *T_r;
  double lambda = Inst.get_arrival_rate(); // mean rate per unit of time within service calls are generated in Poisson manner
  double Mu_NT = Inst.get_service_rate(); // mean of non-travel time component of the service time
  server *serv;
  SQM_solution *Sol;

  logDebug(cout << endl << endl << "*****Start GRASP*****" << endl << endl);
  if (p < 1) return NULL;
  Sol = new SQM_solution(Inst);
  Sol->set_params(lambda,Mu_NT);

  logDebug(cout << "/* Locate the first server */" << endl);
  Sol->add_server();
  Sol->set_speed(v,beta);
  Sol->set_server_location(0,unif(n));
  r = 1;
  T_r = new double [n];
  rcl = new int [n];
  while (r < p) 
    {
      logDebug(cout << "[" << r << "]/* Evaluate posible locations*/" << "\t");
      Sol->add_server();
      Sol->set_speed(v,beta);
      for (int i = 0;i < n;i++) {
	Sol->test_server_location(r,i);
	/* T_r[i] = MST_response_time(Sol);*/
	/* T_r[i] = GRASP_func_NN(Sol);*/
	T_r[i] = GRASP_func_kNN(Sol,min(r+1,GRASP_kNN_param));
      }

      logDebug(cout << "/* Sort Restricted Candidates List */" << "\t");
      sort_dist(n,T_r,rcl);
      logDebug(cout << "/* Choose random element from the rcl */" << "\r");
      element = unif(ceil(alpha * n));
      Sol->set_server_location(r++,rcl[element]);
    }
  logDebug(cout << endl);

  delete [] rcl;
  delete [] T_r;
  logDebug(cout << "Finish GRASP" << endl);
  return Sol;
}

double GRASP_func_pMean (SQM_solution *Sol) {
  logDebug(cout << "Start GRASP_func_pMean" << endl);
  /* Variables definition */
  SQM_instance *I = Sol->get_instance();
  int n = I->demand_points();
  int m = I->potential_sites();
  int p = Sol->get_servers();
  int nearest,k;
  double *rho;
  double Obj;
  /* */

  Obj = 0.0;
  for (int j = 0;j < m;j++) {
    k = GRASP_nearest_server(Sol,p); /* Obtain the nearest server */
    nearest = Sol->get_server_location(k);
    Obj += dist(I->demand(j),I->site(nearest)) / Sol->get_server_speed(k);
  }

  logDebug(cout << "Finish GRASP_func_pMean" << endl);
  return Obj;
}

double GRASP_func_NN (SQM_solution *Sol) {
  logDebug(cout << "Start GRASP_func_NN" << endl);
  /* Variables definition */
  SQM_instance *I = Sol->get_instance();
  int n = I->demand_points();
  int m = I->potential_sites();
  int p = Sol->get_servers();
  int nearest,k;
  double *rho;
  double Obj;
  /* */

  Obj = 0.0;
  for (int j = 0;j < m;j++) {
    k = GRASP_nearest_server(Sol,p); /* Obtain the nearest server */
    nearest = Sol->get_server_location(k);
    Obj += I->demand(j)->demand * dist(I->demand(j),I->site(nearest)) / Sol->get_server_speed(k);
  }
  Obj /= MINS_PER_BLOCK * BLOCKS_PER_HORIZON;

  logDebug(cout << "Finish GRASP_func_NN" << endl);
  return Obj;
}

bool GRASP_closest_to_b(SQM_instance *I,int node,int center_a,int center_b) {
  return (dist(I->demand(node),I->site(center_a)) > dist(I->demand(node),I->site(center_b)));
}

int GRASP_nearest_server(SQM_solution *Sol,int j) {
  SQM_instance *I = Sol->get_instance();
  int p = Sol->get_servers();
  int k = 0;
  for (int i = 1;i < p;i++)
    if (GRASP_closest_to_b(I,j,Sol->get_server_location(k),Sol->get_server_location(i)))
      k = i;
  return k;
}

double GRASP_func_kNN (SQM_solution *Sol,int K) {
  logDebug(cout << "Start GRASP_func_kNN" << endl);
  /* Variable definitions */
  SQM_instance *I = Sol->get_instance();
  int p = Sol->get_servers();
  int m = I->demand_points(); /* Number of demand points */
  double RT = 0.0; /* Response Time */
  double lambda = Sol->get_arrival_rate();
  double Mu_NT = Sol->get_non_travel_time();
  int **a = Sol->preferred_servers();
  double *Lambda;
  double demand;
  logDebug(cout << "Termina definicion de variables" << endl);

  Lambda = new double[m];
  demand = I->total_demand();
  for (int k = 0;k < m;k++) Lambda[k] = I->demand(k)->demand * lambda / demand;

  logDebug(cout << "Comienza calculo de rho_i" << endl);
  /* Calculate the first approach for rho */
  double distance;
  double *rho = new double[p];
  for (int i = 0;i < p;i++) rho[i] = 0.0;
  for (int i = 0;i < p;i++) {
    distance = 0.0;
    for (int k = 0;k < m;k++)
      if (a[k][0] == i) {
	rho[i] += Lambda[k] * (1/Mu_NT + Sol->get_server_rate(i) * Sol->distance(i,k) / (MINS_PER_BLOCK * BLOCKS_PER_HORIZON));
	distance +=  Sol->distance(i,k);
      }
    RT += rho[i] * distance / (Sol->get_server_speed(i) * MINS_PER_BLOCK * BLOCKS_PER_HORIZON);
  }
  
  logDebug(cout << "Comienza calculo de rho_i para resto de ordenes" << endl);
  int t = 1;
  double rho_a_ml;
  double *new_rho = new double [p];
  do {
    for (int i = 0;i < p;i++) new_rho[i] = 0.0;
    for (int i = 0;i < p;i++) {
      distance = 0.0;
      for (int k = 0;k < m;k++)
	if (a[k][t] == i) {
	  rho_a_ml = 1;
	  for (int l = 0;l < t;l++) rho_a_ml *= rho[a[k][l]];
	  new_rho[i] += (1 - rho[i]) * rho_a_ml * Lambda[k] * (1/Mu_NT + Sol->get_server_rate(i) * Sol->distance(i,k) / (MINS_PER_BLOCK * BLOCKS_PER_HORIZON));
	  distance += Sol->distance(i,k);
	}
      RT += new_rho[i] * distance / (Sol->get_server_speed(i) * MINS_PER_BLOCK * BLOCKS_PER_HORIZON);
    }
    
    for (int i = 0;i < p;i++) {
      rho[i] += new_rho[i];
      if (rho[i] > 1)
	logError(cout << "¡rho_" << i+1 << " > 1! in order " << t << endl);
    }
  } while (++t < K &&  t < p);

  logDebug(cout << "Comienza a liberar memoria" << endl);
  delete [] new_rho;    
  delete [] rho;
  delete [] Lambda;
  logDebug(cout << "Finish GRASP_func_kNN" << endl);
  return RT;
}
