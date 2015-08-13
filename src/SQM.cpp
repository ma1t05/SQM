
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
  int M_clients;
  int N_sites;

  if (argc < 8) {
    filename = "Pba";
    M_clients = 50;
    N_sites = 30;
    p = 5;
    mu = 60.0*24.0/20.0;
    f = 0.006;
    v = 40.0;
  }
  else {
    filename = argv[1];
    M_clients = atoi(argv[2]);
    N_sites = atoi(argv[3]);
    p = atoi(argv[4]);
    mu = atof(argv[5]);
    f = atof(argv[6]);
    v = atof(argv[7]);
  }
  
  srand(time(NULL));
  LogName << "SQM_" << M_clients << "_" << N_sites << "_" << p << ".log";
  LogFile.open(LogName.str().c_str());

  /*I = read_points(demad_file.c_str());*/
  /* I = IC_read_instance(demand_file,facility_file); */
  I = IC_create_instance(M_clients,N_sites);
  IC_write_instance(I,filename+"_demand.ins",filename+"_facility.ins");
  /*IC_plot_instance(filename,filename+"_demand.ins",filename+"_facility.ins");*/
  /*SQM_model(I,p,3,mu,f,v);*/
  Goldberg(I,p,mu,f);
  LogFile.close();
  cout << LogName.str() << endl;

  delete[] I->V;
  delete[] I->W;
  delete I;
}

void read_config_file() {
  fstream config;
  
  config.open("SQM.conf",fstream::in);
  
}
