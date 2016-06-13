
#include "common.h"

/* Global Variables */
string Instance_Name;
int M_clients,N_sites;
int p,l;
double v;
/* Global variables read as aguments */
double lambda;
double Mu_NT;

long SQM_solution::calls_to_grt = 0;
clock_t SQM_solution::processing_time = 0;

/* log.h extern variables */
std::ofstream LogFile;
std::ofstream results;
std::ofstream dat;
int log_depth;
/* Simulation.h extern variables */
std::ofstream Log_Simulation;

/* RefSet extern variables */
void (*Improvement_Method)(SQM_solution&);
SolList* (*Combine_Solutions)(SQM_solution&,SQM_solution&);
/* PathRelinking extern variables */
int* (*matching_function)(SQM_solution&,SQM_solution&); /* function for match */
int* (*order_function)(SQM_solution&,int*,SQM_solution&); /* function for proccess */

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
       << "  Local_Search      " << endl
       << endl;

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
	  {"consider"    ,required_argument,0 , 'k'},
	  {"arrival-rate",required_argument,0 , 'l'},
	  {"service-rate",required_argument,0 , 'm'},
	  {"speed"       ,required_argument,0 , 's'},
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
	else if (method == "multi_start")
	  Test_Function = Test_SQM_multistart;
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
	else if (method == "Tune_PR")
	  Test_Function = Tune_Path_Relinking;
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

void Test_SQM_multistart(SQM_instance &Inst,int p,double v) {
  int N = 5000,step = 100;
  SQM_solution *Sol;
  RefSet Top(10);

  for (int i = step;i <= N;i+=step) {
    for (int j = 0;j < step;j++) {
      Sol = new SQM_solution (Inst,p);
      Sol->set_speed(v,BETA);
      if (!Top.TryAdd(*Sol,Sol->get_response_time())) delete Sol;
    }
    Top.clean_garbage();
    cout << "The best solution at iteration " << i << " is: " << Top.best()
	 << endl;
  }
  cout << "  RefSetCall:" << Top.get_Calls() << endl
       << "   RefSetAdd:" << Top.get_Adds() << endl
       << "    DupCheck:" << Top.get_Checks() << endl
       << "FullDupCheck:" << Top.get_FullCheck() << endl
       << "FullDupFound:" << Top.get_DupFound() << endl
       << "    Best Sol:" << Top.best() << endl
       << "   Worst Sol:" << Top.worst() << endl
       << endl;
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
  SQM_heuristic(*Sol);
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
      SQM_heuristic(*G);
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
    SQM_heuristic(*X);
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
  int N = 500,num_elite = 10;
  double rt,worst_rt;
  clock_t beginning,now;
  double seconds;
  SQM_solution *X;
  SubsetControl *SC;
  RefSet *elite_sols;
  double best_multistart_rt,best_rt,improvement;
  matching_type match;
  order_type order;

  /* Determine combination method */
  match = perfect_matching;
  order = nearest_first;
  Combine_Solutions = Path_Relinking;
  Improvement_Method = No_Improvement;
  /* Determine methods to use in Path_Relinking */
  /* PR_{perfect|random|workload}_matching */
  matching_function = PR_perfect_matching; 
  /* PR_processing_order_{random|nf|ff} */
  order_function = PR_processing_order_nf;

  /* Generate Pool */
  beginning = clock();
  SolList pool;
  for (int r = 0;r < N;r++) {
    X = new SQM_solution(Inst,p);
    X->set_speed(v,BETA);
    Improvement_Method(*X);
    X->get_response_time();
    pool.push_back(X);
  }
  now = clock();
  seconds = (double)(now - beginning)/CLOCKS_PER_SEC;
  cout << "Generate " << N << " solutions in " << seconds << " secs." << endl;

  beginning = clock();
  pool.sort(compare_SQMSols);
  now = clock();
  seconds = (double)(now - beginning)/CLOCKS_PER_SEC;
  cout << "Sort solutions in " << seconds << " secs." << endl;

  best_multistart_rt = pool.front()->get_response_time();
  cout << "Best multistart rt: " << best_multistart_rt << endl;

  beginning = clock();
  SC = new TwoTier_SC (num_elite,num_elite,pool);
  now = clock();
  seconds = (double)(now - beginning)/CLOCKS_PER_SEC;

  elite_sols = SC->get_RefSet ();
  best_rt = elite_sols->best();
  improvement = 100*(best_multistart_rt - best_rt) / best_multistart_rt;
  cout << "After " << seconds << " seconds" << endl
       << "the best response time is: " << best_rt << endl
       << " with an % improvement of: " << improvement << endl;

  /* Write Results */
  results.open("results_PathRelinking.csv",std::ofstream::app);
  results
    /* Instance */
    << (*elite_sols)[0]
    /* Path Relinking params */
    << "Path_Relinking,"
    << match << ","
    << order << ","
    /* Data Results */
    << best_multistart_rt << ","
    << best_rt << ","
    << improvement << ","
    << seconds << ","
    << SQM_solution::get_calls_to_grt() << ","
    << SQM_solution::get_processing_time()
    << endl;
  results.close();

  delete SC;
}

void Test_SQM_Local_Search(SQM_instance &Inst,int p,double v) {
  int iLoc;
  int top_sols = 10;
  int N = 500, step = 100;
  SQM_solution *Best,*Sol,*X,*Y;
  double rt,bh_rt,bh_ls_rt,ls_rt,ls_bh_rt;
  clock_t beginning,now;
  clock_t start,end;
  clock_t bh_clocks,bhls_clocks,ls_clocks,lsbh_clocks;
  double ms_seconds,ls_seconds;
  RefSet Top(top_sols);

  logInfo(cout << "Test Local_Search: Start" << endl);

  beginning = clock();
  for (int r = step;r <= N;r+=step) {
    for (int i = 0;i < step;i++) {
      Sol = new SQM_solution(Inst,p);
      Sol->set_speed(v,BETA);
      iLoc = Top.TryAdd(*Sol,Sol->get_response_time());
      if (iLoc < 0) delete Sol;
    }
    logDebug(cout << "Clean garbage" << endl);
    Top.clean_garbage();
  }
  now = clock();
  ms_seconds = (double)(now - beginning)/CLOCKS_PER_SEC;
  logInfo(cout << "Ends multistart in " << ms_seconds << " seconds" << endl);

  beginning = clock();
  rt = bh_rt = bh_ls_rt = ls_rt = ls_bh_rt = 0;
  bh_clocks = bhls_clocks = ls_clocks = lsbh_clocks = 0;
  for (int i = 0;i < top_sols;i++) {
    iLoc = Top.location(i);
    Sol = Top[iLoc];
    X = Sol->clone();
    Y = Sol->clone();

    rt += Sol->get_response_time();

    logDebug(cout << "Step " <<i+1<< " (a)Apply\n\t Berman Heuristic" << endl);
    start = clock();
    SQM_heuristic(*Y);
    end = clock();
    bh_clocks += end - start;
    bh_rt += Y->get_response_time();
    logDebug(cout << "\t Local Search" << endl);
    start = clock();
    Local_Search(*Y);
    end = clock();
    bhls_clocks += end - start;
    bh_ls_rt += Y->get_response_time();

    logDebug(cout << "Step " << i+1 << " (b)Apply\n\t Local Search" << endl);
    start = clock();
    Local_Search(*X);
    end = clock();
    ls_clocks += end - start;
    ls_rt += X->get_response_time();
    logDebug(cout << "\t Berman Heuristic" << endl);
    start = clock();
    SQM_heuristic(*X);	
    end = clock();
    lsbh_clocks += end - start;
    ls_bh_rt += X->get_response_time();

    delete X;
    delete Y;
  }
  now = clock();
  ls_seconds = (double)(now - beginning)/CLOCKS_PER_SEC;
    
  results.open("results_LocalSearch.csv",std::ofstream::app);
  results 
    /* Instance */
    << Sol
    /* Best multistart */
    << rt/top_sols << ","
    /* Berman Heuristic + Local Search */
    << bh_rt/top_sols << ","
    << bh_ls_rt/top_sols << ","
    /* Local Search + Berman Heuristic */
    << ls_rt/top_sols << ","
    << ls_bh_rt/top_sols << ","
    /* Improvement */
    << 100.0*(rt-bh_ls_rt)/rt << ","     /* BH+LS */
    << 100.0*(rt-bh_rt)/rt << ","        /* BH */
    << 100.0*(bh_rt-bh_ls_rt)/rt  << "," /* +LS */
    << 100.0*(rt-ls_bh_rt)/rt << ","     /* LS+BH */
    << 100.0*(rt-ls_rt)/rt << ","        /* LS */
    << 100.0*(ls_rt-ls_bh_rt)/rt << ","  /* +BH */
    /* Processing Times */
    /* Multi Start */
    << ms_seconds << ","
    /* Berman Heuristic & Local Search */
    << ls_seconds << ","
    << (double)bh_clocks/CLOCKS_PER_SEC << ","    /* Berman Heuristic */
    << (double)bhls_clocks/CLOCKS_PER_SEC << ","  /* + Local Search */
    << (double)ls_clocks/CLOCKS_PER_SEC << ","    /* Local Search */
    << (double)lsbh_clocks/CLOCKS_PER_SEC         /* + Berman Heuristic */
    << endl;
  results.close();

  logInfo
    (cout 
     << "Ends Local Search after " << ls_seconds << " seconds" << endl
     << "       Method: " << "Response Time\t" << "% Improvement" << endl
     << "response time: " << rt/top_sols << endl
     << "    heuristic: " << bh_rt/top_sols << "\t" << 100.0*(rt-bh_rt)/rt
     << endl
     << "+local search: " << bh_ls_rt/top_sols << "\t+" << 100.0*(bh_rt-bh_ls_rt)/rt
     << endl
     << " local search: " << ls_rt/top_sols << "\t" << 100.0*(rt-ls_rt)/rt
     << endl
     << "   +heuristic: " << ls_bh_rt/top_sols << "\t+" << 100.0*(ls_rt-ls_bh_rt)/rt
     << endl
     );
  logInfo(cout << "Test Local_Search: Finish" << endl);
}

void Tune_Path_Relinking(SQM_instance &Inst,int p,double v) {
  int N = 500,num_elite = 10;
  double rt,worst_rt;
  clock_t beginning,now;
  double secs;
  SQM_solution *X;
  SubsetControl *SC;
  RefSet *elite_sols;
  double best_multistart_rt,best_rt,improvement;
  matching_type match;
  order_type order;

  logInfo(cout << "++ Start Path-relinking ++" << endl);

  Combine_Solutions = Path_Relinking;
  Improvement_Method = No_Improvement;

  /* Determine combination method */
  match = perfect_matching;
  order = nearest_first;
  do {
    set_match_method(match);
    set_order_method(order);

    logInfo(cout << "match: " << match << endl
	    << "order: " << order << endl);
    /* Generate Pool */
    beginning = clock();
    SolList pool;
    for (int r = 0;r < N;r++) {
      X = new SQM_solution(Inst,p);
      X->set_speed(v,BETA);
      Improvement_Method(*X);
      X->get_response_time();
      pool.push_back(X);
    }
    now = clock();
    secs = (double)(now - beginning)/CLOCKS_PER_SEC;
    logInfo
      (cout << "Generate " << N << " solutions in " << secs << " secs." << endl);

    beginning = clock();
    pool.sort(compare_SQMSols);
    now = clock();
    secs = (double)(now - beginning)/CLOCKS_PER_SEC;
    logInfo(cout << "Sort solutions in " << secs << " secs." << endl);

    best_multistart_rt = pool.front()->get_response_time();
    logInfo(cout << "Best multistart rt: " << best_multistart_rt << endl);

    beginning = clock();
    SC = new TwoTier_SC (num_elite,num_elite,pool);
    now = clock();
    secs = (double)(now - beginning)/CLOCKS_PER_SEC;

    elite_sols = SC->get_RefSet ();
    best_rt = elite_sols->best();
    improvement = 100*(best_multistart_rt - best_rt) / best_multistart_rt;
    logInfo
      (cout
       << "After " << secs << " seconds" << endl
       << "the best response time is: " << best_rt << endl
       << " with an % improvement of: " << improvement << endl);

    /* Write Results */
    results.open("results_PathRelinking.csv",std::ofstream::app);
    results
      /* Instance */
      << (*elite_sols)[0]
      /* Path Relinking params */
      << "Path_Relinking,"
      << match << ","
      << order << ","
      /* Data Results */
      << best_multistart_rt << ","
      << best_rt << ","
      << improvement << ","
      << secs << ","
      << SQM_solution::get_calls_to_grt() << ","
      << SQM_solution::get_processing_time()
      << endl;
    results.close();

    delete SC;
  } while ((++order != invalid_order) ||
	   (!(order = nearest_first) && (++match != invalid_matching)));
}
