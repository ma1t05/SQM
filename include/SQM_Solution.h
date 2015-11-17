
#ifndef _SQM_SOLUTION_H
#define _SQM_SOLUTION_H 1

#include "SQM_Solution"

/* Note response_unit was changed to server */
class server {
  int location;
  int past_location;
  double v;
  double beta;
}

class SQM_solution {
 private:
  SQM_instance *I;
  int p;
  server *Servers;
  int **a;
 public:
  SQM_solution(SQM_instance *I,int p);
  ~SQM_solution();
  void set_server_location(int i,int j);
  int get_server_location(int i);
}

#endif

/* eof */
