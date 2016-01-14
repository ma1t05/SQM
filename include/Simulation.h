
#ifndef _SIMULATION_H
#define _SIMULATION_H

#include <list>
#include <iostream>
#include <fstream>
#include "SQM_heuristic.h"

#define Start_Time 7
#define Simulation_Time 10

extern std::ofstream Log_Simulation;
typedef list<event*> events_list;

class event {
public:
  virtual double get_time () = 0;
  virtual void process (status&) = 0;
};

class call : public event {
private:
  double at_time;
  int demand_point;
  bool queued;
  double queued_at_time;
public:
  call (double t,int i) : at_time(t), demand_point(i) {};
  double get_time () const;
  void process (status&);
  bool is_queued () const;
  double get_waiting_time () const;
  void queue_call (event_list&);
};
std::ostream& operator<<(std::ostream &os,call &c);

class release : public event {
private:
  int server;
  int demand_point;
  double arrival_time;
  double service_time;
  double return_time;
public:
  release (call&);
  double get_time () const;
  double get_arrival_time () const;
  void process (status&);
};
std::ostream& operator<<(std::ostream &os,release &r);

void Simulator(SQM_instance *I,int p,double lambda,double Mu_NT,double v);
list<event*>* Generate_calls(SQM_instance *I,double lambda);

struct status {
  SQM_solution *Sol;
  bool *busy;
  double current_time;
  double *busy_time;
  list<event*> *events;
  list<event*> queue;
  int calls_sent_to_queue;
  double waiting_time;
  int total_calls;
  double response_time;
  double service_time;
};

typedef struct status status;


#endif

/* eof */
