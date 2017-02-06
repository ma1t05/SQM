
#include "common.h"

command *Test_Function;

int main(int argc,char **argv) {
  stringstream LogName;
  SQM_instance *I;

  process_command_line(argc,argv);
  read_config_file("SQM.conf");
  srand(time(NULL));

  /* Open Log File */
  LogName << "SQM_" << M_clients << "_" << N_sites << "_" << p << ".log";
  LogFile.open(LogName.str().c_str(),std::ofstream::app);
  log_depth = 0;
  logInfo(cout << "Creates LogFile" << endl);

  I = SQM_load_instance(Instance_Name,M_clients,N_sites);
  logInfo(cout << "Call Test Function" << endl);
  Test_Function(*I,p,v);
  delete I;

  /* Log */ LogFile.close();
  logInfo(cout << endl << "Saved in LogFile: " << LogName.str() << endl);

  return 0;
}
