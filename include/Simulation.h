
#ifndef _SIMULATION_H
#define _SIMULATION_H

#include <list>
#include <fstream>
#include "SQM_heuristic.h"

enum type_event {CALL,RELEASE,QUEUING};

struct event {
  type_event type;
  int node;
  double at_time;
};

typedef struct event event;

struct status {
  SQM_solution *Sol;
  bool *busy;
  double current_time;
  double *busy_time;
  list<event*> *events;
  list<event*> queue;
};

typedef struct status status;

void Simulator(SQM_instance *I,int p,double lambda,double Mu_NT,double v);
list<event*>* Generate_calls(SQM_instance *I,double lambda);

extern std::ofstream Log_Simulation;

#endif

/* eof */
