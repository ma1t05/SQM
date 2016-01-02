
#ifndef _SQM_SOLUTION_H
#define _SQM_SOLUTION_H 1

#include <iostream>
#include "SQM_Instance.h"

#define UNASIGNED_LOCATION -1

extern double BETA;

class server {
 private:
  int location;
  int past_location;
  double v;
  double beta;
 public:
  server();
  server(int i);
  ~server();
  void set_speed (double v,double beta);
  void set_location (int i);
  void test_location (int i);
  int get_location () const;
  int get_past_location () const;
  double get_speed () const;
  double get_beta () const;
  double get_rate () const;
};

class SQM_solution {
 private:
  SQM_instance *Inst;
  int p;
  server *Servers;
  int **a;
  double response_time;
  double lambda;
  double Mu_NT;
 public:
  SQM_solution (SQM_instance *I);
  SQM_solution (SQM_instance *I,int p);
  SQM_solution (SQM_solution *Sol);
  ~SQM_solution ();
  SQM_solution* clone();
  double pm_cost; /* cost of perfect matching */
  void set_params(double lambda,double Mu_NT);
  void set_speed (double v,double beta);
  void set_server_location (int i,int j);
  void test_server_location (int i,int j);
  void add_server ();
  void remove_server (int);
  void update_preferred_servers ();
  void delete_preferred_servers ();
  int get_servers () const;
  int get_server_location (int i) const;
  int get_server_past_location (int i) const;
  SQM_instance* get_instance () const;
  double get_server_speed (int i) const;
  double get_server_beta (int i) const;
  double get_server_rate (int i) const;
  double get_arrival_rate () const;
  double get_non_travel_time () const;
  double get_response_time ();
  double distance (int i,int k) const;
  int ** preferred_servers () const;
  bool operator<(SQM_solution&);
  bool operator>(SQM_solution&);
};

ostream& operator<<(ostream&,SQM_solution*);

#endif

/* eof */
