
#include <ctime>
#include <iostream>

using namespace std;
#include "config.h"
#include "Goldberg.h"
#include "SQM_model.h"
#include "SQM_GRASP.h"
#include "PathRelinking.h"
#include "Local_Search.h"

/* log.h extern variables */
std::ofstream LogFile;
std::ofstream results;
std::ofstream dat;

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

struct call {
  double at_time;
  int node;
};

typedef struct call call;

struct release {
  double at_time;
  int server;
};

typedef struct release release;

struct status {
  SQM_solution *Sol;
  bool *busy;
  double current_time;
  list<call*> queue;
};

typedef struct status status;

void read_config_file(string configFile);
void Test_MST(SQM_instance *I,int p,double lambda,double Mu_NT,double v);
void Test_exponential(SQM_instance *I,int p,double lambda,double Mu_NT,double v);
void Simulator(SQM_instance *I,int p,double lambda,double Mu_NT,double v);
list<call*>* Generate_calls(SQM_instance *I,double lambda);

void Simulator_release_server(status &state,release *rel);
void Simulator_address_call(status &state,call *event);

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
  mu = 60.0*24.0/20.0; f = 4; v = 500.0;

  srand(time(NULL));

  /* Open Log File */
  LogFile.open("Test_SQM.log",std::ofstream::app);
  I = SQM_load_instance(filename,M_clients,N_sites);
  Test_exponential(I,p,f,mu,v);
  delete I;
  /* Log */ LogFile.close();
  logInfo(cout << endl << "Saved in LogFile: " << LogName.str() << endl);
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

void Simulator(SQM_instance *I,int p,double lambda,double Mu_NT,double v) {
  SQM_solution *X = new SQM_solution(I,p);
  status state;
  list<release*> releases;
  list<call*> *calls = Generate_calls(I,lambda);
  list<release*>::iterator rel;

  X->set_speed(v,BETA);
  X->set_params(lambda,Mu_NT);
  SQM_heuristic(X);
  state.Sol = X;
  state.busy = new bool[p];
  state.current_time = 0;

  rel = releases.begin();
  for (list<call*>::iterator it = calls->begin();it != calls->end();it++) {
    //cout << (*it)->node << "\t" << (*it)->at_time << endl;
    while (rel != releases.end() && (*rel)->at_time < (*it)->at_time) {
      Simulator_release_server(state,*rel);
      rel++;
    }
    Simulator_address_call(state,*it);
    delete *it;
  }
  delete calls;
  
}

void Simulator_release_server(status &state,release *rel) {
  state.busy[rel->server] = false;
  state.current_time = rel->at_time;
  // Revisar si queue esta vacio
  // Si no lo esta seleccionar llamada a atender
     // llamar Simulator_address_call
}

void Simulator_address_call(status &state,call *event) {
  // Buscar servidor mas cercano
  /* Si hay servidor
     marcarlo como ocupado
     agregar release a lista
     borrar evento
  */
  /* Si no hay servido
     enviar evento a la cola
  */
}

list<call*>* Generate_calls(SQM_instance *I,double lambda) {
  int m;
  int events;
  double Simulation_Time = 10;
  double current_time,etime;
  double demand,lambda_j;
  list<call*> *calls;
  list<call*>::iterator it;
  call *event;

  logDebug(cout << "Start Generate_calls" << endl);
  calls = new list<call*>;
  m = I->demand_points();
  demand = I->total_demand();
  for (int j = 0;j < m;j++) {
    lambda_j = lambda * I->get_demand(j) / demand;
    current_time = exponential(lambda_j);
    it = calls->begin();
    events = 0;
    while (current_time < Simulation_Time) {
      /* Create event */
      logDebug(cout << "Create event" << endl);
      event = new call;
      event->at_time = current_time;
      event->node = j;
      /* Insert event */
      logDebug(cout << "Insert event with lambda_j = " << lambda_j << endl);
      while (it != calls->end() && current_time > (*it)->at_time) it++;
      if (it != calls->end())
	calls->insert(it,event);
      else calls->push_back(event);
      /* Determine time of next event */
      etime = exponential(lambda_j);
      current_time += etime;
      events++;
      logDebug(cout << "Determine time of next event " << current_time << endl);
    }
    logDebug(cout << "Create " << events << " for demand point " << j << endl);
  }
  
  return calls;
}
