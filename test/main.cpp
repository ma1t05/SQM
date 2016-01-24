
#include <ctime>
#include <iostream>
#include <unistd.h> /* getopt */
#include <getopt.h> /* getopt_long */

using namespace std;
#include "config.h"
#include "Goldberg.h"
#include "SQM_model.h"
#include "SQM_GRASP.h"
#include "PathRelinking.h"
#include "Local_Search.h"

/* Global Variables */
string Instance_Name;
int M_clients,N_sites;
int p,l;
double v;
/* Global variables read as aguments */
double lambda;
double Mu_NT;

void print_usage ();
void process_command_line(int,char**);
void read_config_file(string configFile);
void Log_Start_SQMH(int M_clients,int N_sites,int p,double mu,double f);
typedef list<SQM_solution*> Solutions;

/* Test Functions */
void Test_SQM_model(SQM_instance&,int,double);
void Test_SQM_heuristic(SQM_instance&,int,double);
void Test_SQM_GRASP(SQM_instance&,int,double);
void Test_SQM_random(SQM_instance&,int,double);
void Test_SQM_Path_Relinking(SQM_instance&,int,double);
void Test_SQM_Local_Search(SQM_instance&,int,double);
/*  */
void (*Test_Function)
(SQM_instance&, /* The instance */
 int,           /* The number of servers */
 double         /* The speed */
 );

/* log.h extern variables */
std::ofstream LogFile;
std::ofstream results;
std::ofstream dat;
/* Simulation.h extern variables */
std::ofstream Log_Simulation;

/* PathRelinking extern variables */
int* (*matching_function)(SQM_solution*,SQM_solution*); /* function for match */
int* (*order_function)(SQM_solution*,int*,SQM_solution*); /* function for proccess */

/* Global variables read from config */
int MINS_PER_BLOCK;
int BLOCKS_PER_HORIZON;
double BETA;
double JARVIS_EPSILON;
double MIN_RANGE_X;
double MIN_RANGE_Y;
double MAX_RANGE_X;
double MAX_RANGE_Y;
int GRASP_kNN_param;
double EPSILON;

double TIME_MAX;

int main(int argc,char **argv) {
  stringstream LogName;
  SQM_instance *I;

  process_command_line(argc,argv);
  read_config_file("SQM.conf");
  srand(time(NULL));

  /* Open Log File */
  LogName << "SQM_" << M_clients << "_" << N_sites << "_" << p << ".log";
  LogFile.open(LogName.str().c_str(),std::ofstream::app);

  I = SQM_load_instance(Instance_Name,M_clients,N_sites);
  Test_Function(*I,p,v);
  delete I;

  /* Log */ LogFile.close();
  logInfo(cout << endl << "Saved in LogFile: " << LogName.str() << endl);

  return 0;
}

void print_usage () {
  cout << "Usage: SQM [options] <command>" << endl;

  /* Print Options */
  cout << "Options:" << endl
       << "  -f FILE, --prefix=FILE      " 
       << "Use FILE as file prefix, for read/write instance." << endl
       << "  -m M, --demand=M            " 
       << "Sets the number of Demand Points to M." << endl
       << "  -n N, --sites=N             " 
       << "Sets the number of Potential Sites to N." << endl
       << "  -p P, --servers=P           " 
       << "Sets the number of servers to P." << endl
       << "  -k K, --consider=K          " 
       << "For commands that apply, only considered the K" << endl
       << "                              closest servers." << endl
       << "  -l LAMBDA, --arrival-rate=LAMBDA" << endl 
       << "                              " 
       << "Sets the General Arrival Rate to LAMBDA." << endl
       << "  -m MU, --service-rate=MU    " 
       << "Sets the on scene Service Rate to MU." << endl
       << "  -s SPEED, --speed=SPEED     " 
       << "Sets the speed of transfer to SPEED." << endl
       << endl;

  /* Print available commands */
  cout << "The available commands are:" << endl
       << "  model             " << endl
       << "  heuristic         " << endl
       << "  GRASP             " << endl
       << "  random            " << endl
       << "  Path_Relinking    " << endl
       << "  Local_Search      " << endl;

}

void process_command_line(int argc,char **argv) {
  int c;
  static int verbose_flag;

  /* Default Values */
  Instance_Name = "Test";
  M_clients = 50; 
  N_sites = 50;
  p = 10; 
  l = 10;
  Mu_NT = 3;
  lambda = 6;
  v = 500.0;

  verbose_flag = LOG_INFO;
  while (1) 
    {
      int long_index = 0;
      static struct option long_options[] =
	{
	  {"verbose"     ,no_argument      ,&verbose_flag , LOG_DEBUG},
	  {"brief"       ,no_argument      ,&verbose_flag , LOG_INFO },
	  {"superbrief"  ,no_argument      ,&verbose_flag , LOG_QUIET},
	  {"help"        ,no_argument      ,0 , 'h'},
	  {"prefix"      ,required_argument,0 , 'f'},
	  {"demand"      ,required_argument,0 , 'M'},
	  {"sites"       ,required_argument,0 , 'N'},
	  {"servers"     ,required_argument,0 , 'p'},
	  {"aux"         ,required_argument,0 , 'k'},
	  {"arrival-rate",required_argument,0 , 'l'},
	  {"service-rate",required_argument,0 , 'm'},
	  {"spped"       ,required_argument,0 , 's'},
	  {0             ,0                ,0 , 0  }
	};

      c = getopt_long (argc,argv,"hf:M:N:p:k:l:m:s:",long_options,&long_index);
      if (c == -1) break;
      switch (c) 
	{
	case 0: /* If this option set a flag, do nothing else now. */
	  if (long_options[long_index].flag != 0)
	    break;
	  cout << "option " << long_options[long_index].name;
	  if (optarg)
	    cout << " with arg " << optarg;
	  cout << endl;
	  break;
	case 'f': Instance_Name = string(optarg); /* file-prefix */
	  break;
	case 'M': M_clients = atoi(optarg);
	  break;
	case 'N': N_sites = atoi(optarg);
	  break;
	case 'p': p = atoi(optarg);
	  break;
	case 'k': l = atoi(optarg);
	  break;
	case 'l': lambda = atof(optarg);
	  break;
	case 'm': Mu_NT = atof(optarg);
	  break;
	case 's': v = atof(optarg);
	  break;
	case 'h': print_usage ();
	  exit (0);
	case '?':
	  /* getopt_long already printed an error message */
	  cout << "Option " << optarg << endl;
	  break;
	default: print_usage ();
	  abort ();
	}
    }
  
  switch (verbose_flag) 
    {
    case LOG_QUIET: logLevel = LOG_QUIET;
      break;
    case LOG_INFO: logLevel = LOG_INFO;
      break;
    case LOG_DEBUG: logLevel = LOG_DEBUG;
      break;
    default: logLevel = LOG_INFO;
      break;
    }

  /* Print any remaining command line arguments (not options). */
  if (optind < argc)
    {
      while (optind < argc) {
	string method = string(argv[optind]);
	if (method == "model")
	  Test_Function = Test_SQM_model;
	else if (method == "heuristic")
	  Test_Function = Test_SQM_heuristic;
	else if (method == "GRASP")
	  Test_Function = Test_SQM_GRASP;
	else if (method == "random")
	  Test_Function = Test_SQM_random;
	else if (method == "Path_Relinking")
	  Test_Function = Test_SQM_Path_Relinking;
	else if (method == "Local_Search")
	  Test_Function = Test_SQM_Local_Search;
	else
	  Test_Function = Test_SQM_heuristic;
	optind++;
      }
      cout << endl;
    }
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

void Test_SQM_model(SQM_instance &Inst,int p,double v) {
  int *Sol;
  char sub[16];

  results.open("results_SQM_model.csv",std::ofstream::app);
  results << Inst.demand_points()
	  << "," << Inst.potential_sites()
	  << "," << p 
	  << "," << l
	  << "," << Inst.get_service_rate()
	  << "," << Inst.get_arrival_rate()
	  << "," << v
	  << "," << Instance_Name << "_demand.ins"
	  << "," << Instance_Name << "_facility.ins";
  
  Sol = SQM_model(Inst,p,l,v);
  sprintf(sub,"_%02d_%02d",p,l);
  plot_instance_solution(Inst,Sol,Instance_Name+sub);
  delete[] Sol;

  if (l == p) {
    results << Inst.demand_points()
	    << "," << Inst.potential_sites()
	    << "," << p 
	    << "," << l
	    << "," << Inst.get_service_rate()
	    << "," << Inst.get_arrival_rate()
	    << "," << v
	    << "," << Instance_Name << "_demand.ins" 
	    << "," << Instance_Name << "_facility.ins";
    Goldberg(Inst,p,Mu_NT,lambda);
  }
  results.close();
}

void Test_SQM_heuristic(SQM_instance &Inst,int p,double v) {
  SQM_solution *Sol;
  logInfo(cout << "Calling SQM_Heuristic" << endl);
  Sol = new SQM_solution (Inst,p);
  Sol->set_speed(v,BETA);
  Sol->set_params(lambda,Mu_NT);
  cout << Sol->get_response_time() << endl;
  for (int i = 0;i < p;i++) 
    cout << Sol->get_server_location(i) << "\t";
  cout << endl;
  SQM_heuristic(Sol);
  cout << Sol->get_response_time() << endl;
  for (int i = 0;i < p;i++) 
    cout << Sol->get_server_location(i) << "\t";
  cout << endl;
  if (LogInfo) {
    for (int i = 0;i < p;i++) LogFile << Sol->get_server_location(i) << "\t";
    LogFile << endl;
  }
  delete Sol;
}

/* Experiment to determine the necesary GRASP iterations */
void Test_SQM_GRASP(SQM_instance &Inst,int p,double v) {
  int m = Inst.demand_points();
  int n = Inst.potential_sites();
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
      G = GRASP(Inst,p,v,BETA,alpha); /* */
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

void Test_SQM_random(SQM_instance &Inst,int p,double v) {
  int m = Inst.demand_points();
  int n = Inst.potential_sites();
  clock_t beginning,now;
  double T_r1,T_r2,t_r,BRT;
  double gap = 0.0,best_rt = 100.0,worst_rt = 0.0,avg_rt = 0.0,avg;
  int N = 1000;
  SQM_solution *Best,*Best_RS,*Best_GRASP,*BEST_GRASP;
  SQM_solution *X,*G;
  char GRASP_output[32];
  //int *Sol = new int [n];
  Solutions *best_solutions;
  best_solutions = new Solutions;

  cout << "      Total Demand : " << Inst.total_demand() << endl;

  /* Evaluate GRASP */
  BEST_GRASP = NULL;
  results.open("GRASP_results.csv",std::ofstream::app);
  dat.open("GRASP.dat",std::ofstream::out);
  for (double alpha = 0.0;alpha < 0.99;alpha += 0.05) {
    results << m << "," << n << "," << p << ","
	    << Mu_NT << "," << lambda << "," << alpha << "," << N << ",";
    beginning = clock();
    Best_GRASP = NULL;
    best_rt = 100.0,worst_rt = 0.0,avg_rt = 0.0,avg = 0.0;
    for (int r = 0;r < N;r++) {
      G = GRASP(Inst,p,v,BETA,alpha); /* */
      T_r1 = G->get_response_time();
      avg += T_r1;
      /* Log */ Log_Start_SQMH(m,n,p,Mu_NT,lambda); /* */
      SQM_heuristic(G);
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
    dat << alpha << " " << best_rt << " " << avg_rt/N << " " << avg/N << " " << worst_rt << endl;
    results << best_rt << "," << avg_rt / N << "," << worst_rt << "," 
	    << (double)(now - beginning)/CLOCKS_PER_SEC << endl;
    cout << "\t *** \t GRASP results alpha = " << alpha << "\t ***" << endl;
    cout << "   Best Response time : " << best_rt << endl
	 << "Average Response time : " << avg_rt / N << endl
	 << "  Worst Response time : " << worst_rt << endl
	 << "           time (sec) : " << (double)(now - beginning)/CLOCKS_PER_SEC << endl;
    if (BEST_GRASP == NULL || t_r < BRT) {
      if (BEST_GRASP != NULL) delete BEST_GRASP;
      BEST_GRASP = Best_GRASP;
      BRT = t_r;
    }
    else delete Best_GRASP;
  }
  
  /* Plot Best GRASP Solution *//*
  logDebug(cout << "Plot Best GRASP Solution" << endl);
  for (int i = 0;i < p;i++) cout << BEST_GRASP->get_server_location(i) << " ";
  cout << endl;
  for (int k = 0;k < n;k++) Sol[k] = 0;
  for (int i = 0;i < p;i++) Sol[BEST_GRASP->get_server_location(i)]++;
  logDebug(cout << "Calling plot method" << endl);
  plot_instance_solution(I,Sol,"SQM_Best_GRASP_Sol");*/

  logInfo(cout << "Start SQM random" << endl);
  beginning = clock();
  Best_RS = NULL;
  best_rt = 100.0,worst_rt = 0.0,avg_rt = 0.0,avg = 0.0;
  for (int r = 0;r < N;r++) {
    X = new SQM_solution(Inst,p);
    X->set_speed(v,BETA);
    X->set_params(lambda,Mu_NT);
    T_r1 = X->get_response_time();
    avg += T_r1;
    if (LogDebug) {
      for (int k = 0;k < n;k++)
	for (int i = 0;i < p;i++)
	  if (X->get_server_location(i) == k) cout << k << " ";
      cout << endl;
      cout << "Response time : " << T_r1 << endl;
    }
    /* Log */ Log_Start_SQMH(m,n,p,Mu_NT,lambda); /* */
    SQM_heuristic(X);
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
  dat << "1.0 " << best_rt << " " << avg_rt/N << " " << avg/N << " " << worst_rt << endl;
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

  /* Plot Random Best Solution *//*
  for (int k = 0;k < n;k++) Sol[k] = 0;
  for (int i = 0;i < p;i++) Sol[Best_RS->get_server_location(i)]++;
  plot_instance_solution(I,Sol,"SQM_Best_Random_Sol");*/

  T_r1 = Best_RS->get_response_time();
  T_r2 = BEST_GRASP->get_response_time();
  Best = (T_r1 < T_r2 ? Best_RS : BEST_GRASP);

  /* Plot Best Solution *//*
  for (int k = 0;k < n;k++) Sol[k] = 0;
  for (int i = 0;i < p;i++) Sol[Best->get_server_location(i)]++;
  plot_instance_solution(I,Sol,"SQM_Best_Sol");*/
  sprintf(GRASP_output,"./plots/GRASP_%d_%d_%d_%d",m,n,p,rand());
  gnuplot_GRASP(GRASP_output);

  //delete [] Sol;
  delete Best_RS;
  delete BEST_GRASP;
}

void Test_SQM_Path_Relinking(SQM_instance &Inst,int p,double v) {
  int N = 400,num_elite = 10;
  double rt,worst_rt;
  clock_t beginning,now;
  SQM_solution *X,*Best;
  Solutions *elite_sols,*various_sols;

  worst_rt = 100.0;
  beginning = clock();
  elite_sols = new Solutions;
  for (int r = 0;r < N;r++) {
    X = new SQM_solution(Inst,p);
    X->set_speed(v,BETA);
    SQM_heuristic(X);
    rt = X->get_response_time();
    if (rt < worst_rt) {
      if (elite_sols->size() == num_elite) {
	delete elite_sols->back();
	elite_sols->pop_back();
	worst_rt = elite_sols->back()->get_response_time();
      }

      if (elite_sols->size() > 0) {
	for (Solutions::iterator it = elite_sols->begin(),end = elite_sols->end();it != end;it++)
	  if (**it > *X) {
	    elite_sols->insert(it,X);
	    break;
	  }
      }
      else {
	elite_sols->push_back(X);
	worst_rt = X->get_response_time();
      }
    }
    else delete X;
  }
  cout << endl;

  now = clock();
  cout << "After " << (double)(now - beginning)/CLOCKS_PER_SEC << " seconds" << endl;
  cout << "The best " << num_elite << " response times:" << endl;
  for (Solutions::iterator it = elite_sols->begin();it != elite_sols->end();it++)
    cout << (*it)->get_response_time() << endl;
  Best = elite_sols->front();

  various_sols = new Solutions;
  for (int r = 0;r < 5*num_elite;r++) {
    X = new SQM_solution(Inst,p);
    X->set_speed(v,BETA);
    X->pm_cost = SQM_min_cost_pm(elite_sols,X);
    if (various_sols->size() > 0) {
      Solutions::iterator it;
      for (it = various_sols->begin();it != various_sols->end();it++)
	if ((*it)->pm_cost < X->pm_cost) {
	  various_sols->insert(it,X);
	  break;
	}
      if (it == various_sols->end())
	various_sols->push_back(X);

      if (various_sols->size() > num_elite) {
	delete various_sols->back();
	various_sols->pop_back();
      }
    }
    else various_sols->push_back(X);
  }
  
  cout << "The various " << num_elite << " response times:" << endl;
  for (Solutions::iterator it = various_sols->begin(), end = various_sols->end();it != end;it++) {
    cout << (*it)->get_response_time() << endl;
    elite_sols->push_back(*it);
  }
  delete various_sols;

  int Case = 1;
  double time,gap;

  results.open("results_PathRelinking.csv",std::ofstream::app);
  for (matching_type matching = perfect_matching;matching != invalid_matching;++matching) {
    for (order_type order = nearest_first;order != invalid_order;++order) {
      results << Best;

      switch (matching) 
	{
	case perfect_matching:
	  matching_function = PR_run_perfect_matching;
	  results << "perfect,";
	  break;
	case workload_matching:
	  matching_function = PR_workload_matching;
	  results << "workload,";
	  break;
	case random_matching:
	  matching_function = PR_random_matching;
	  results << "random,";
	  break;
	}

      switch (order)
	{
	case nearest_first:
	  order_function = PR_processing_order_nf;
	  results << "nf,";
	  break;
	case farthest_first:
	  order_function = PR_processing_order_ff;
	  results << "ff,";
	  break;
	case random_order:
	  order_function = PR_processing_order_random;
	  results << "random,";
	  break;
	}

      beginning = clock();
      X = SQM_path_relinking(elite_sols);
      //SQM_heuristic(X);
      gap = 100*(Best->get_response_time() - X->get_response_time()) / Best->get_response_time();
      now = clock();
      time = (double)(now - beginning)/CLOCKS_PER_SEC;
      results << time << "," << gap << "," << X->get_response_time() << endl;
      delete X;
    }
  }
  results.close();

  SQM_delete_sols(elite_sols);
}

void Test_SQM_Local_Search(SQM_instance &Inst,int p,double v) {
  int N = 1;
  SQM_solution *X,*Y;
  double rt,h_rt,ls_rt;

  cout << "Start: Local_Search" << endl;
  for (int r = 0;r < N;r++) {
    X = new SQM_solution(Inst,p);
    X->set_speed(v,BETA);
    rt = X->get_response_time();
    Y = X->clone();
    cout << "\t        Method: " << "Response Time\t" << "% Improvement" << endl
	 << "\t response time: " << rt << endl;

    SQM_heuristic(Y);
    h_rt = Y->get_response_time();
    Local_Search(*Y);
    ls_rt = Y->get_response_time();
    cout << "\t     heuristic: " << h_rt << "\t" << 100.0*(rt-h_rt)/rt << endl
	 << "\t +local search: " << ls_rt << "\t+" << 100.0*(h_rt-ls_rt)/rt << endl;

    Local_Search(*X);
    ls_rt = X->get_response_time();
    SQM_heuristic(X);	
    h_rt = X->get_response_time();
    cout << "\t  local search: " << ls_rt << "\t" << 100.0*(rt-ls_rt)/rt << endl 
	 << "\t    +heuristic: " << h_rt << "\t+" << 100.0*(ls_rt-h_rt)/rt << endl;
    
    cout << endl;
    delete X;
    delete Y;
  }
}
