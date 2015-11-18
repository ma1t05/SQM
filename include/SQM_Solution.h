
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
  void set_speed(double v,double beta);
  int get_location();
  int get_past_location();
  void set_location(int i);
};

class SQM_solution {
 private:
  SQM_instance *Inst;
  int p;
  server *Servers;
  int **a;
 public:
  SQM_solution (SQM_instance *I);
  SQM_solution (SQM_instance *I,int p);
  ~SQM_solution ();
  SQM_instance* get_instance ();
  int get_servers ();
  void set_speed (double v,double beta);
  void set_server_location (int i,int j);
  int get_server_location (int i);
  double get_server_speed(int i);
  void add_server ();
};

#endif

/* eof */
