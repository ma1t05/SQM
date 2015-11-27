
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
/*#include "SQM_model.h"
  #include "Goldberg.h"*/
#include "SQM_heuristic.h"
#include "config.h"
#include "gnuplot.h"
#include "SQM_GRASP.h"
#include "log.h"
#include "PathRelinking.h"
#include <list>
using namespace std;

std::ofstream LogFile;
std::ofstream results;
std::ofstream dat;

void Log_Start_SQMH(int M_clients,int N_sites,int p,double mu,double f);
void Call_SQM_heuristic(SQM_instance* I,int p,double f,double mu);
void Call_SQM_GRASP(SQM_instance *I,int p,double lambda,double Mu_NT,double v);
void Call_SQM_random(SQM_instance *I,int p,double lambda,double Mu_NT,double v);

int main(int argc,char *argv[]) {
  string filename;
  int M_clients,N_sites;
  int p,l;
  double mu,f,v;
  stringstream LogName;
  SQM_instance *I;

  if (argc < 9) {
    filename = "Pba";
    M_clients = 50; N_sites = 30;
    p = 5; l = 3;
    mu = 60.0*24.0/20.0; f = 0.006; v = 40.0;
  }
  else {
    filename = argv[1];
    M_clients = atoi(argv[2]);
    N_sites = atoi(argv[3]);
    p = atoi(argv[4]);
    l = atoi(argv[5]);
    mu = atof(argv[6]);
    f = atof(argv[7]);
    v = atof(argv[8]);
  }
  
  srand(time(NULL));

  /* Open Log File */
  LogName << "SQM_" << M_clients << "_" << N_sites << "_" << p << ".log";
  LogFile.open(LogName.str().c_str(),std::ofstream::app);

  I = SQM_load_instance(filename,M_clients,N_sites);
  // Call_SQM_model(I,p,l,f,mu,v,filename);
  // Call_SQM_GRASP(I,p,f,mu,v);
  Call_SQM_random(I,p,f,mu,v);
  /* Log Log_Start_SQMH(M_clients,N_sites,p,mu,f); /* */
  // Call_SQM_heuristic(I,p,f,mu);
  delete I;
  /* Log */ LogFile.close();
  logInfo(cout << endl << "Saved in LogFile: " << LogName.str() << endl);
}

void read_config_file() {
  fstream config;
  
  config.open("SQM.conf",fstream::in);
  
}

void Log_Start_SQMH(int M_clients,int N_sites,int p,double mu,double f) {
  LogFile << "*** Start SQM Heuristic ***" << endl
	  << "with "
	  << M_clients << " clients, "
	  << N_sites << " sites, "
	  << p << " servers, "
	  << mu << " mean service time, "
	  << f << " arrival rate"
	  << endl;
}

void Call_SQM_heuristic(SQM_instance* I,int p,double f,double mu) {
  SQM_solution *Sol;
  cout << "Calling SQM_Heuristic" << endl;
  Sol = new SQM_solution (I,p);
  SQM_heuristic(Sol,f,mu);
  for (int i = 0;i < p;i++) LogFile << Sol->get_server_location(i) << " ";
  LogFile << endl;
  delete Sol;
}

/* Experiment to determine the necesary GRASP iterations */
void Call_SQM_GRASP(SQM_instance *I,int p,double lambda,double Mu_NT,double v) {
  int m = I->demand_points();
  int n = I->potential_sites();
  double beta = 1.5;
  double T_r1,T_r2,t_r;
  double best_rt = 100.0,worst_rt = 0.0,avg_rt = 0.0;
  int N = 1000;
  SQM_solution *X,*G;
  char GRASP_output[32];

  /* Evaluate GRASP */
  for (double alpha = 0.0;alpha < 0.99;alpha += 0.05) {
    sprintf(GRASP_output,"./plots/GRASP_iterations_%d_%d_%d_%0.2f.dat",m,n,p,alpha);
    dat.open(GRASP_output,std::ofstream::out);
    best_rt = 100.0,worst_rt = 0.0,avg_rt = 0.0;
    for (int r = 0;r < N;r++) {
      G = GRASP(I,p,lambda,Mu_NT,v,alpha); /* */
      T_r1 = G->get_response_time();
      avg_rt += T_r1;
      if (best_rt > T_r1) best_rt = T_r1;
      if (worst_rt < T_r1) worst_rt = T_r1;
      if (r%10 == 0)
	dat << r << " " << best_rt << " " << avg_rt/(r+1) << " " << worst_rt << endl;
      delete G;
    }
    dat.close();
  }

  
}

void Call_SQM_random(SQM_instance *I,int p,double lambda,double Mu_NT,double v) {
  int m = I->demand_points();
  int n = I->potential_sites();
  clock_t beginning,now;
  double beta = 1.5;
  double T_r1,T_r2,t_r,BRT;
  double gap = 0.0,best_rt = 100.0,worst_rt = 0.0,avg_rt = 0.0;
  int N = 100;
  SQM_solution *Best,*Best_RS,*Best_GRASP,*BEST_GRASP;
  SQM_solution *X,*G;
  char GRASP_output[32];
  int *Sol = new int [n];
  list<SQM_solution*> *best_solutions;
  best_solutions = new list<SQM_solution*>;

  cout << "      Total Demand : " << I->total_demand() << endl;

  /* Evaluate GRASP */
  BEST_GRASP = NULL;
  results.open("GRASP_results.csv",std::ofstream::app);
  dat.open("GRASP.dat",std::ofstream::out);
  for (double alpha = 0.0;alpha < 0.99;alpha += 0.05) {
    results << m << "," << n << "," << p << ","
	    << Mu_NT << "," << lambda << "," << alpha << "," << N << ",";
    beginning = clock();
    Best_GRASP = NULL;
    best_rt = 100.0,worst_rt = 0.0,avg_rt = 0.0;
    for (int r = 0;r < N;r++) {
      G = GRASP(I,p,lambda,Mu_NT,v,alpha); /* */
      T_r1 = G->get_response_time();
      /* Log */ Log_Start_SQMH(m,n,p,Mu_NT,lambda); /* */
      SQM_heuristic(G,lambda,Mu_NT);
      T_r2 = G->get_response_time();
      avg_rt += T_r2;
      if (best_rt > T_r2) best_rt = T_r2;
      if (worst_rt < T_r2) worst_rt = T_r2;

      if (Best_GRASP == NULL || T_r2 < t_r) {
	if (Best_GRASP != NULL) delete Best_GRASP;
	Best_GRASP = G;
	t_r = T_r2;
      }
      else delete G;
    }
    now = clock();
    dat << alpha << " " << best_rt << " " << avg_rt/N << " " << worst_rt << endl;
    results << best_rt << "," << avg_rt / N << "," << worst_rt << "," 
	    << (double)(now - beginning)/CLOCKS_PER_SEC << endl;
    cout << "\t *** \t GRASP results alpha = " << alpha << "\t ***" << endl;
    cout << "   Best Response time : " << best_rt << endl
	 << "Average Response time : " << avg_rt / N << endl
	 << "  Worst Response time : " << worst_rt << endl
	 << "           time (sec) : " << (double)(now - beginning)/CLOCKS_PER_SEC << endl;
    if (BEST_GRASP == NULL || t_r < BRT) {
      /*if (BEST_GRASP != NULL) delete BEST_GRASP;*/
      BEST_GRASP = Best_GRASP;
      BRT = t_r;
    }
    /*else delete Best_GRASP;*/
    best_solutions->push_back(Best_GRASP);
  }
  
  /* Plot Best GRASP Solution */
  logDebug(cout << "Plot Best GRASP Solution" << endl);
  for (int i = 0;i < p;i++) cout << BEST_GRASP->get_server_location(i) << " ";
  cout << endl;
  for (int k = 0;k < n;k++) Sol[k] = 0;
  for (int i = 0;i < p;i++) Sol[BEST_GRASP->get_server_location(i)]++;
  logDebug(cout << "Calling plot method" << endl);
  plot_instance_solution(I,Sol,"SQM_Best_GRASP_Sol");

  logInfo(cout << "Start SMQ random" << endl);
  beginning = clock();
  Best_RS = NULL;
  best_rt = 100.0,worst_rt = 0.0,avg_rt = 0.0;
  for (int r = 0;r < N;r++) {
    X = new SQM_solution(I,p);
    X->set_speed(v,beta);
    X->set_params(lambda,Mu_NT);
    T_r1 = X->get_response_time();
    if (LogDebug) {
      for (int k = 0;k < n;k++)
	for (int i = 0;i < p;i++)
	  if (X->get_server_location(i) == k) cout << k << " ";
      cout << endl;
      cout << "Response time : " << T_r1 << endl;
    }
    /* Log */ Log_Start_SQMH(m,n,p,Mu_NT,lambda); /* */
    SQM_heuristic(X,lambda,Mu_NT);
    T_r2 = X->get_response_time();

    avg_rt  += T_r2;
    if (best_rt > T_r2) best_rt = T_r2;
    if (worst_rt < T_r2) worst_rt = T_r2;

    if (Best_RS == NULL || T_r2 < t_r) {
      if (Best_RS != NULL) delete Best_RS;
      Best_RS = X;
      t_r = T_r2;
    }
    else delete X;
  }
  best_solutions->push_back(Best_RS);
  now = clock();
  dat << "1.0 " << best_rt << " " << avg_rt/N << " " << worst_rt << endl;
  dat.close();
  results << m << "," << n << "," << p << ","
	  << Mu_NT << "," << lambda << ",1.0," << N << ","
	  << best_rt << "," << avg_rt / N << "," << worst_rt << "," 
	  << (double)(now - beginning)/CLOCKS_PER_SEC << endl;
  results.close();
  cout << "\t *** \tRandom results\t ***" << endl;
  cout << "   Best Response time : " << best_rt << endl
       << "Average Response time : " << avg_rt / N << endl
       << "  Worst Response time : " << worst_rt << endl
       << "           time (sec) : " << (double)(now - beginning)/CLOCKS_PER_SEC << endl;

  /* Plot Random Best Solution */
  for (int k = 0;k < n;k++) Sol[k] = 0;
  for (int i = 0;i < p;i++) Sol[Best_RS->get_server_location(i)]++;
  plot_instance_solution(I,Sol,"SQM_Best_Random_Sol");

  T_r1 = Best_RS->get_response_time();
  T_r2 = BEST_GRASP->get_response_time();
  Best = (T_r1 < T_r2 ? Best_RS : BEST_GRASP);

  /* Plot Best Solution */
  for (int k = 0;k < n;k++) Sol[k] = 0;
  for (int i = 0;i < p;i++) Sol[Best->get_server_location(i)]++;
  plot_instance_solution(I,Sol,"SQM_Best_Sol");
  sprintf(GRASP_output,"./plots/GRASP_%d_%d_%d_%d",m,n,p,rand());
  gnuplot_GRASP(GRASP_output);

  Best = SQM_path_relinking(best_solutions);
  /* Plot Best Solution */
  for (int k = 0;k < n;k++) Sol[k] = 0;
  for (int i = 0;i < p;i++) Sol[Best->get_server_location(i)]++;
  plot_instance_solution(I,Sol,"SQM_Best_Sol_PR");
  sprintf(GRASP_output,"./plots/GRASP_%d_%d_%d_%d",m,n,p,rand());
  gnuplot_GRASP(GRASP_output);

  delete [] Sol;
  delete Best;
}
