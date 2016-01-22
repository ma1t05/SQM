
#include <ctime>
#include <iostream>

using namespace std;
#include "config.h"
#include "Goldberg.h"
#include "SQM_model.h"
#include "SQM_GRASP.h"
#include "PathRelinking.h"
#include "Local_Search.h"
#include "Simulation.h"

/* log.h extern variables */
std::ofstream LogFile;
std::ofstream results;
std::ofstream dat;
/* Simulation.h exter variable */
std::ofstream Log_Simulation;

/* PathRelinking extern variables */
int* (*matching_function)(SQM_solution*,SQM_solution*); /* function for match */
int* (*order_function)(SQM_solution*,int*,SQM_solution*); /* function for proccess */

/* Global extern variables read from config */
int MINS_PER_BLOCK;
int BLOCKS_PER_HORIZON;
/* SQM_Solution */
double BETA;
/* SQM_Instance */
double MIN_RANGE_X;
double MIN_RANGE_Y;
double MAX_RANGE_X;
double MAX_RANGE_Y;
/* mp_jarvis */
double JARVIS_EPSILON;
/* SQM_GRASP */
int GRASP_kNN_param;
/* cplex */
double EPSILON;
double TIME_MAX;
/* Global exter variables to read from arguments */
double lambda;
double Mu_NT;

void read_config_file(string configFile);
void Test_MST(SQM_instance *I,int p,double lambda,double Mu_NT,double v);
void Test_exponential(SQM_instance *I,int p,double lambda,double Mu_NT,double v);

int main(int argc,char *argv[]) {
  string filename;
  int M_clients,N_sites;
  int p,l;
  double v;
  stringstream LogName;
  SQM_instance *I;

  read_config_file("SQM.conf");
  filename = "Test";
  M_clients = 50; N_sites = 50;
  p = 10; l = 10;
  Mu_NT = MINS_PER_BLOCK * 3 / 60.0;
  lambda = MINS_PER_BLOCK * 6 / 60.0;
  v = 500.0;

  srand(time(NULL));

  /* Open Log File */
  Log_Simulation.open("Simulation.log",std::ofstream::out);
  I = SQM_load_instance(filename,M_clients,N_sites);
  //Test_exponential(I,p,lambda,Mu_NT,v);
  Simulator(*I,p,v);
  delete I;
  /* Log */ Log_Simulation.close();

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

void Test_MST(SQM_instance *I,int p,double lambda,double Mu_NT,double v) {
  int loc;
  SQM_solution *X;
  double rt;
  cout << "Start: Test_MST" << endl;
  X = new SQM_solution(*I,p);
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
  times = new double [m];
  for (int j = 0;j < m;j++) {
    lambda_j = lambda * I->get_demand(j) / demand;
    times[j] = exponential(lambda_j);
  }
  a = new int [m];
  sort_dist(m,times,a);

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

}
