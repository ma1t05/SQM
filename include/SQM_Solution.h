
#ifndef _SQM_SOLUTION_H
#define _SQM_SOLUTION_H 1

#include "SQM_Instance.h"

/* Note response_unit was changed to server */
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
  int get_location ();
  int get_past_location ();
  double get_speed ();
  double get_beta ();
  double get_rate ();
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
  void set_params(double lambda,double Mu_NT);
  void set_speed (double v,double beta);
  void set_server_location (int i,int j);
  void test_server_location (int i,int j);
  void add_server ();
  void update_preferred_servers ();
  int get_servers ();
  int get_server_location (int i);
  int get_server_past_location (int i);
  SQM_instance* get_instance ();
  double get_server_speed (int i);
  double get_server_beta (int i);
  double get_server_rate (int i);
  double get_response_time ();
  int ** preferred_servers ();
  bool operator<(SQM_solution&);
  bool operator>(SQM_solution&);
};


#endif

/* eof */
