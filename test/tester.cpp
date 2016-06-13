
#include "common.h"

int main(int argc,char *argv[]) {
  string filename;
  int M_clients,N_sites;
  int p,l;
  double v;
  stringstream LogName;
  SQM_instance *I;

  read_config_file("SQM.conf");
  filename = "Test";
  M_clients = 50; 
  N_sites = 30;
  p = 6; 
  l = 3;
  Mu_NT = MINS_PER_BLOCK * 3 / 60.0;
  lambda = MINS_PER_BLOCK * 6 / 60.0;
  v = 500.0;

  srand(time(NULL));

  /* Open Log File */
  Log_Simulation.open("Simulation.log",std::ofstream::out);
  I = SQM_load_instance(filename,M_clients,N_sites);
  //Test_exponential(*I,p,v);
  //Simulator(*I,p,v);
  Test_Path_Relinking(*I,p,v);
  delete I;
  /* Log */ Log_Simulation.close();

}

void read_config_file(string configFile) {
  char *envp = NULL;
  Config SQM_conf(configFile,&envp);
  /* Sochastic Queue p-Medina Problem Params */
  MINS_PER_BLOCK = SQM_conf.pInt("MINS_PER_BLOCK");
  BLOCKS_PER_HORIZON = SQM_conf.pInt("BLOCKS_PER_HORIZON");
  BETA = SQM_conf.pInt("BETA");

  /*  Parameter for Jarvis Approximation procedure */
  JARVIS_EPSILON = SQM_conf.pDouble("JARVIS_EPSILON");

  /* Ranges to create random instances */
  MIN_RANGE_X = SQM_conf.pDouble("MIN_RANGE_X");
  MAX_RANGE_X = SQM_conf.pDouble("MAX_RANGE_X");
  MIN_RANGE_Y = SQM_conf.pDouble("MIN_RANGE_Y");
  MAX_RANGE_Y = SQM_conf.pDouble("MAX_RANGE_Y");

  GRASP_kNN_param = SQM_conf.pInt("GRASP_kNN_param");

  /* Cplex params */
  EPSILON = SQM_conf.pDouble("EPSILON");
  TIME_MAX = SQM_conf.pDouble("TIME_MAX");
  
}

void Test_MST(SQM_instance &I,int p,double v) {
  int loc;
  SQM_solution *X;
  double rt;
  cout << "Start: Test_MST" << endl;
  X = new SQM_solution(I,p);
  X->set_speed(v,BETA);
  X->set_params(lambda,Mu_NT);
  for (int i = 0;i < p;i++) {
    rt = X->get_response_time();
    cout << rt << endl;
    loc = X->get_server_location(0);
    X->remove_server(0);
    X->add_server();
    X->set_server_location(p-1,loc);
    X->set_speed(v,BETA);
    cout << X->get_response_time() << "\t" << X->get_response_time() - rt << endl;
  }
  delete X;
}

void Test_exponential(SQM_instance &I,int p,double v) {
  int m,*a;
  double demand;
  double lambda_j;
  double *times;

  m = I.demand_points();
  demand = I.total_demand();
  times = new double [m];
  for (int j = 0;j < m;j++) {
    lambda_j = lambda * I.get_demand(j) / demand;
    times[j] = exponential(lambda_j);
  }
  a = new int [m];
  sort_dist(m,times,a);

  cout << "Site\tRate\t\tTime" << endl;
  for (int j = 0;j < m;j++) {
    lambda_j = lambda * I.get_demand(j) / demand;
    cout << setfill('0') << setw(3) << a[j] << "\t" 
	 << setfill(' ')
	 << setw(10) << lambda_j << "\t" 
	 << setw(10) << times[a[j]] << endl;
  }
  
  delete [] a;
  delete [] times;

}

void Test_Path_Relinking(SQM_instance &I,int p,double v) {
  SQM_solution *X,*Y;
  int *match,*order;
  X = new SQM_solution(I,p);
  X->set_speed(v,BETA);
  Y = X->clone();
  SQM_heuristic(Y);
  match = PR_perfect_matching(X,Y);
  order = PR_processing_order_nf(X,match,Y);
  X->update_preferred_servers ();
  gnuplot_Path_Relinking(*X,match,*Y,"step_00");
  for (int step = 0;step < p;step++) {
    int x = order[step];
    int y = match[x];
    int loc_x = X->get_server_location(x);
    int loc_y = Y->get_server_location(y);
    if (loc_x != loc_y) {
      X->set_server_location(x,loc_y);
      X->update_preferred_servers ();
      char step_file[8];
      sprintf(step_file,"step_%02d",step);
      gnuplot_Path_Relinking(*X,match,*Y,string(step_file));
    }
  }
  delete [] order;
  delete [] match;
  delete Y;
  delete X;
}

void Test_Path_Relinking(SQM_instance &I,int p,double v) {
  SQM_solution *X,*Y;
  int *match,*order;
  X = new SQM_solution(I,p);
  X->set_speed(v,BETA);
  Y = X->clone();
  SQM_heuristic(Y);
  int best_loc = UNASIGNED_LOCATION;
  double best_rt;
  int p = X.get_servers();
  /* Obtain the server with less workload */
  int j = LS_get_server_with_less_workload(X);
  int loc_j  = X.get_server_location(j);
  /* obtain the news workloads */
  int k = LS_get_server_with_more_workload(X);
  /* Put a server near the server with more workload */
  Sites *lst = LS_get_adjacent_sites(X,k);
  for (Sites::iterator it = lst->begin(),end = lst->end(); it != end;it++) {
    X.test_server_location(j,*it);
    if (best_loc == UNASIGNED_LOCATION || X.get_response_time() < best_rt) {
      best_loc = *it;
      best_rt = X.get_response_time();
    }
  }

  if (best_loc == UNASIGNED_LOCATION) {
    cerr << "Warining: No new location" << endl;
    best_loc = loc_j;
  }
  X.test_server_location(j,loc_j); /* The past location of the server */
  X.set_server_location(j,best_loc);
  delete lst;
}
