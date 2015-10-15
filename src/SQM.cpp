
#include "SQM.h"
#include "SQM_model.h"
#include "Goldberg.h"
#include "SQM_heuristic.h"
#include "config.h"

std::ofstream LogFile;
std::ofstream results;

bool file_exists (const string&);
SQM_instance* Load_instance(string filename,int M_clients,int N_sites);
void Call_SQM_heuristic(SQM_instance* I,int p,double f,double mu);
void Call_SQM_model(SQM_instance* I,int p,int l,double f,double mu,double v,string filename);
void Log_Start_SQMH(int M_clients,int N_sites,int p,double mu,double f);
void Call_SQM_random(SQM_instance *I,int p,double lambda,double Mu_NT);

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
  /* Log Log_Start_SQMH(M_clients,N_sites,p,mu,f); /* */
  // Call_SQM_heuristic(I,p,f,mu);
  Call_SQM_random(I,p,f,mu);
  IC_delete_instance(I);
  /* Log */ LogFile.close();
  cout << LogName.str() << endl;

  delete[] I->V;
  delete[] I->W;
  delete I;
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
  /*I = read_points(demad_file.c_str());*/
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
  X = SQM_heuristic(I,p,f,mu);
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
  IC_plot_instance(I,Sol,filename+sub);
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

void Call_SQM_random(SQM_instance *I,int p,double lambda,double Mu_NT) {
  double T_r;
  cout << "Start SMQ random" << endl;
  response_unit *X = guess_a_location_03(p,I->N,I->W);
  T_r = SQM_response_time(I,p,X,lambda,Mu_NT);
  cout << "Response time : " << T_r << endl;
  delete [] X;
}
