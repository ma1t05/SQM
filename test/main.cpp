
#include <iostream>
#include "SQM.h"
using namespace std;

/* log.h extern variables */
std::ofstream LogFile;
std::ofstream results;
std::ofstream dat;

/* PathRelinking extern variables */
int* (*matching_function)(SQM_solution*,SQM_solution*); /* function for match */
int* (*order_function)(SQM_solution*,int*,SQM_solution*); /* function for proccess */

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
  //Call_SQM_random(I,p,f,mu,v);
  //Call_SQM_Path_Relinking(I,p,f,mu,v);
  Call_SQM_Local_Search(I,p,f,mu,v);
  /* Log  Log_Start_SQMH(M_clients,N_sites,p,mu,f); /* */
  //Call_SQM_heuristic(I,p,f,mu,v);
  delete I;
  /* Log */ LogFile.close();
  logInfo(cout << endl << "Saved in LogFile: " << LogName.str() << endl);
}
