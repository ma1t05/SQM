
#include "SQM.h"
#include "SQM_model.h"
#include "Goldberg.h"
#include "SQM_heuristic.h"
#include "config.h"
#include "gnuplot.h"
#include "SQM_GRASP.h"
#include "MST.h"

std::ofstream LogFile;
std::ofstream results;

bool file_exists (const string&);
SQM_instance* Load_instance(string filename,int M_clients,int N_sites);
void Call_SQM_heuristic(SQM_instance* I,int p,double f,double mu);
void Call_SQM_model(SQM_instance* I,int p,int l,double f,double mu,double v,string filename);
void Log_Start_SQMH(int M_clients,int N_sites,int p,double mu,double f);
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

  I = Load_instance(filename,M_clients,N_sites);
  // Call_SQM_model(I,p,l,f,mu,v,filename);
  Call_SQM_random(I,p,f,mu,v);
  /* Log Log_Start_SQMH(M_clients,N_sites,p,mu,f); /* */
  // Call_SQM_heuristic(I,p,f,mu);
  IC_delete_instance(I);
  /* Log */ LogFile.close();
  cout << endl << "Saved in LogFile: " << LogName.str() << endl;
}

void read_config_file() {
  fstream config;
  
  config.open("SQM.conf",fstream::in);
  
}

bool file_exists (const string& name) {
  if (FILE *file = fopen(name.c_str(), "r")) {
    fclose(file);
    return true;
  }
  return false;
}

SQM_instance* Load_instance(string filename,int M_clients,int N_sites) {
  SQM_instance *I;
  string Path = "../git/PMCLAP/Instancias/Q_MCLP_";
  /*I = read_points(demad_file.c_str());*/
  if ((M_clients == N_sites) && (filename == "Q_MCLP")) {
    switch (M_clients) {
    case 30 : 
      filename = Path+"30.txt";
      break;
    case 324 : 
      filename = Path+"324.txt";
      break;
    case 818 : 
      filename = Path+"818.txt";
      break;
    case 3283 : 
      filename = Path+"3283.txt";
      break;
    defautl:
      return NULL;
    }
    cout << "Read file: " << filename << endl;
    I = IC_load_instance(filename);
    return I;
  }
  if (file_exists(filename+"_demand.ins") &&
      file_exists(filename+"_facility.ins")) {
    I = IC_read_instance (filename+"_demand.ins",filename+"_facility.ins");
  }
  else {
    I = IC_create_instance(M_clients,N_sites);
    IC_write_instance(I,filename+"_demand.ins",filename+"_facility.ins");
  }
  return I;
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
  response_unit *X;
  cout << "Calling SQM_Heuristic" << endl;
  X = guess_a_location_03(p,I->N,I->W);
  SQM_heuristic(I,p,f,mu,X);
  for (int i = 0;i < p;i++) LogFile << X[i].location << " ";
  LogFile << endl;
  delete [] X;
}

void Call_SQM_model(SQM_instance* I,int p,int l,double f,double mu,double v,string filename) {
  int *Sol;

  results.open("results.csv",std::ofstream::app);
  results << I->M
	  << "," << I->N
	  << "," << p 
	  << "," << l
	  << "," << mu
	  << "," << f 
	  << "," << v;

  results << "," << filename << "_demand.ins," << filename << "_facility.ins";
  
  Sol = SQM_model(I,p,l,mu,f,v);
  char sub[16];
  sprintf(sub,"_%02d_%02d",p,l);
  plot_instance_solution(I,Sol,filename+sub);
  delete[] Sol;

  if (l == p) {
    results << I->M
	    << "," << I->N
	    << "," << p 
	    << "," << l
	    << "," << mu
	    << "," << f 
	    << "," << v;
    results << "," << filename << "_demand.ins," << filename << "_facility.ins";
    Goldberg(I,p,mu,f);
  }
  results.close();
}

void Call_SQM_random(SQM_instance *I,int p,double lambda,double Mu_NT,double v) {
  double beta = 1.5;
  double T_r1,T_r2,t_r;
  double gap = 0.0,best_gap = -1.0,worst_gap = 1.0,avg_gap = 0.0;
  int N = 10;
  response_unit *Best,*Best_RS,*Best_GRASP;
  response_unit *X,*G;

  // cout << "Start SMQ random" << endl;
  Best_RS = NULL;
  for (int r = 0;r < N;r++) {
    X = guess_a_location_03(p,I->N,I->W);
    for (int i = 0;i < p;i++) {
      X[i].v = v;
      X[i].beta = beta;
    }
    T_r1 = MST_response_time(I,p,X,lambda,Mu_NT);
    /*
    for (int k = 0;k < I->N;k++)
      for (int i = 0;i < p;i++)
	if (X[i].location == k) cout << k << " ";
    cout << endl;
    */
    // cout << "Response time : " << T_r << endl;
    /* Log */ Log_Start_SQMH(I->M,I->N,p,Mu_NT,lambda); /* */
    SQM_heuristic(I,p,lambda,Mu_NT,X);
    T_r2 = MST_response_time(I,p,X,lambda,Mu_NT);

    gap = (T_r1 - T_r2) / T_r1;
    avg_gap  += gap;
    if (best_gap < gap) best_gap = gap;
    if (worst_gap > gap) worst_gap = gap;

    if (Best_RS == NULL || T_r2 < t_r) {
      if (Best_RS != NULL) delete [] Best_RS;
      Best_RS = X;
      t_r = T_r2;
    }
    else delete [] X;
  }

  cout << "\t *** \tRandom results\t ***" << endl;
  cout << "Best Response time : " << t_r << endl;
  cout << "          Best gap : " << 100 * best_gap << endl
       << "       Average gap : " << 100 * avg_gap / N << endl
       << "         Worst gap : " << 100 * worst_gap << endl;

  /* Plot Random Best Solution */
  int *Sol = new int [I->N];
  for (int k = 0;k < I->N;k++) Sol[k] = 0;
  for (int i = 0;i < p;i++) Sol[Best_RS[i].location]++;
  plot_instance_solution(I,Sol,"SQM_Best_Random_Sol");

  /* Evaluate GRASP */
  Best_GRASP = NULL;
  best_gap = -1.0,worst_gap = 1.0,avg_gap = 0.0;
  for (int r = 0;r < N;r++) {
    G = GRASP(I,p,lambda,Mu_NT,v,0.25); /* */
    T_r1 = MST_response_time(I,p,G,lambda,Mu_NT);
    /* Log */ Log_Start_SQMH(I->M,I->N,p,Mu_NT,lambda); /* */
    SQM_heuristic(I,p,lambda,Mu_NT,G);
    T_r2 = MST_response_time(I,p,G,lambda,Mu_NT);

    gap = (T_r1 - T_r2) / T_r1;
    avg_gap += gap;
    if (best_gap < gap) best_gap = gap;
    if (worst_gap > gap) worst_gap = gap;

    if (Best_GRASP == NULL || T_r2 < t_r) {
      if (Best_GRASP != NULL) delete [] Best_GRASP;
      Best_GRASP = G;
      t_r = T_r2;
    }
    else delete [] G;
  }

  cout << "\t *** \t GRASP results\t ***" << endl;
  cout << "Best Response time : " << t_r << endl;
  cout << "          Best gap : " << 100 * best_gap << endl
       << "       Average gap : " << 100 * avg_gap / N << endl
       << "         Worst gap : " << 100 * worst_gap << endl;

  /* Plot Best GRASP Solution */
  for (int k = 0;k < I->N;k++) Sol[k] = 0;
  for (int i = 0;i < p;i++) Sol[Best_GRASP[i].location]++;
  plot_instance_solution(I,Sol,"SQM_Best_GRASP_Sol");

  T_r1 = MST_response_time(I,p,Best_RS,lambda,Mu_NT);
  T_r2 = MST_response_time(I,p,Best_GRASP,lambda,Mu_NT);
  Best = (T_r1 < T_r2 ? Best_RS : Best_GRASP);

  /* Plot Best Solution */
  for (int k = 0;k < I->N;k++) Sol[k] = 0;
  for (int i = 0;i < p;i++) Sol[Best[i].location]++;
  plot_instance_solution(I,Sol,"SQM_Best_Sol");

  /*
  double demand;
  demand = 0.0;
  for (int k = 0;k < I->M;k++) demand += (I->V)[k].demand;
  cout << "      Total Demand : " << demand << endl;
  */

  delete [] Sol;
  delete [] Best_RS;
  delete [] Best_GRASP;
}
