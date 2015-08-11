
#include "SQM.h"
#include "SQM_model.h"
#include "config.h"

std::ofstream LogFile;

int main(int argc,char *argv[]) {
  string filename;
  stringstream LogName;
  SQM_instance *I;
  int p;
  double mu;
  double f;
  double v;
  char *env;
  env = NULL;
  
  Config config("SQM.conf",&env);
  char *penv;
  cout << "Termina de leer archivo de configuraciones" << endl;
  int M_clients = config.pInt("M"); 
  penv = getenv("clients");
  if (penv != NULL) M_clients = atoi(penv);
  int N_sites = config.pInt("N");
  penv = getenv("facilities");
  if (penv != NULL) N_sites = atoi(penv);

  if (argc < 6) {
    filename = "./../PMCLAP/Instancias/Q_MCLP_30.txt";
    p = 5;
    mu = 60.0*24.0/20.0;
    f = 0.016;
    v = 40.0;
  }
  else {
    filename = argv[1];
    p = atoi(argv[2]);
    mu = atof(argv[3]);
    f = atof(argv[4]);
    v = atof(argv[5]);
  }
  
  LogName << "SQM_" << M_clients << "_" << N_sites << "_" << p << ".log";
  LogFile.open(LogName.str().c_str());

  /*I = read_points(demad_file.c_str());*/
  /* I = IC_read_instance(demand_file,facility_file); */
  I = IC_create_instance(M_clients,N_sites);
  IC_write_instance(I,filename+"_demand.ins",filename+"_facility.ins");
  SQM_model(I,p,3,mu,f,v);
  LogFile.close();

  delete[] I->V;
  delete[] I->W;
  delete I;
}

void read_config_file() {
  fstream config;
  
  config.open("SQM.conf",fstream::in);
  
}
