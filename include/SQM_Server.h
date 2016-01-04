
#ifndef _SQM_SERVER_H
#define _SQM_SERVER_H 1

#include "SQM_Instance.h"
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

#endif

/* eof */
