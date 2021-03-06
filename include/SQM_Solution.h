
#ifndef _SQM_SOLUTION_H
#define _SQM_SOLUTION_H 1

#include <iostream>
#include "MST.h"


class SQM_solution {
 private:
  SQM_instance *Inst;
  int p;
  server *Servers;
  int **a;
  double response_time;
  double lambda;
  double Mu_NT;
  static long calls_to_grt;
  static clock_t processing_time;
 public:
  SQM_solution (SQM_instance&);
  SQM_solution (SQM_instance&,int);
  SQM_solution (SQM_solution&);
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
  double Hash() const;
  double get_response_time ();
  double* get_workload ();
  server* get_Servers();
  double distance (int i,int k) const;
  int ** preferred_servers () const;
  bool operator<(SQM_solution&);
  bool operator>(SQM_solution&);
  bool operator==(SQM_solution&);
  bool operator!=(SQM_solution&);
  static long get_calls_to_grt () {return calls_to_grt;};
  static double get_processing_time () {
    return double(processing_time)/CLOCKS_PER_SEC;
  };
};

ostream& operator<<(ostream&,SQM_solution*);

#endif

/* eof */
