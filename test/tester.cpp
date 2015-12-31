
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

void Test_MST(SQM_instance *I,int p,double lambda,double Mu_NT,double v);
void Test_exponential(SQM_instance *I,int p,double lambda,double Mu_NT,double v);

int main(int argc,char *argv[]) {
  string filename;
  int M_clients,N_sites;
  int p,l;
  double mu,f,v;
  stringstream LogName;
  SQM_instance *I;

  filename = "Test";
  M_clients = 50; N_sites = 50;
  p = 10; l = 10;
  mu = 60.0*24.0/20.0; f = 40; v = 500.0;

  srand(time(NULL));

  /* Open Log File */
  LogFile.open("Test_SQM.log",std::ofstream::app);
  I = SQM_load_instance(filename,M_clients,N_sites);
  Test_MST(I,p,f,mu,v);
  delete I;
  /* Log */ LogFile.close();
  logInfo(cout << endl << "Saved in LogFile: " << LogName.str() << endl);
}

void Test_MST(SQM_instance *I,int p,double lambda,double Mu_NT,double v) {
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

void Test_exponential(SQM_instance *I,int p,double lambda,double Mu_NT,double v) {
  int m,*a;
  double demand;
  double lambda_j;
  double *times;

  m = I->demand_points();
  demand = I->total_demand();
  /*
  times = new double [m];
  for (int j = 0;j < m;j++) {
    lambda_j = lambda * I->get_demand(j) / demand;
    times[j] = exponential(lambda_j);
  }
  a = new int [m];
  sort_dist(m,times,a);

  /*
  cout << "Site\tRate\t\tTime" << endl;
  for (int j = 0;j < m;j++) {
    lambda_j = lambda * I->get_demand(j) / demand;
    cout << setfill('0') << setw(3) << a[j] << "\t" 
	 << setfill(' ')
	 << setw(10) << lambda_j << "\t" 
	 << setw(10) << times[a[j]] << endl;
  }
  
  delete [] a;
  delete [] times;
  */

  double current_time,etime;
  int it = 0;
  current_time = 0;
  lambda_j = lambda * I->get_demand(0) / demand;
  cout << "Lambda_0 = " << lambda_j << endl;
  while (current_time < 100) {
    etime = exponential(lambda_j);
    current_time += etime;
    cout << ++it << "\t" << current_time << endl;
  }
}
