
#ifndef _SIMULATION_H
#define _SIMULATION_H

#include <list>
#include <iostream>
#include <fstream>
#include "SQM_heuristic.h"

#define Start_Time 7
#define Simulation_Time 10

enum event_type {CALL,RELEASE};

class event {
 protected:
  event_type type;
  double at_time;
  int node;
 public:
  event (event_type et,double t,int i) : type(et), at_time(t) ,node(i) {};
  event_type get_type ();
  double get_time ();
  int get_node ();
};

std::ostream& operator<<(std::ostream &os,event &e);

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
