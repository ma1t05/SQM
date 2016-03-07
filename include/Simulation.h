
#ifndef _SIMULATION_H
#define _SIMULATION_H

#include <list>
#include <iostream>
#include <fstream>
#include "SQM_Solution.h"

#define Start_Time 7
#define Simulation_Time 10

#define NOT_QUEUED -1
extern std::ofstream Log_Simulation;
typedef struct status status;

class event {
public:
  virtual void process (status&) = 0;
  virtual int get_point () const = 0;
  virtual double get_time () const = 0;
};
typedef list<event*> list_events;

class call : public event {
private:
  double at_time;
  int demand_point;
  bool queued;
  double queued_at_time;
  void queue_call (status&);
  void set_time (double t);
public:
  call (double t,int i);
  void process (status&);
  void dequeue_call (status&,int);
  double get_time () const;
  int get_point () const;
  bool is_queued () const;
  double get_waiting_time () const;
};
std::ostream& operator<<(std::ostream &os,call &c);
typedef list<call*> queued_calls;

class release : public event {
private:
  int server;
  int demand_point;
  double travel_time;
  double on_scene_st;
  double follow_up_travel_time;
  double service_time;
  double return_time;
public:
  release (status&,int,int);
  void process (status&);
  int get_point () const;
  double get_time () const;
  double get_travel_time () const;
};
std::ostream& operator<<(std::ostream &os,release &r);

void Simulator(SQM_instance &Inst,int p,double v);
list_events* Generate_calls(SQM_instance &Inst,double lambda);

struct status {
  SQM_solution *Sol;
  bool *busy;
  double current_time;
  double *busy_time;
  list_events *events;
  queued_calls queue;
  int total_calls;
  int calls_sent_to_queue;
  double waiting_time;
  double arrival_time;
  double service_time;
};

extern void (*Improvement_Method)(SQM_solution&);

#endif

/* eof */
