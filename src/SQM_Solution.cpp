
#include "SQM_Solution.h"

SQM_solution::SQM_solution (SQM_instance *I) {
  Inst = I;
  Servers = new server[p];
  p = 0;
}

SQM_solution::SQM_solution (SQM_instance *I,int servers) {
  bool *option;
  int location,locations = 0;
  Inst = I;
  p = servers;
  Servers = new server[p];
  option = new bool[n];
  for (int i = 0;i < n;i++) option[i] = false;
  do {
    location = unif(n);
    if (option[location] == false) {
      locations++;
      option[location] = true;
    }
  } while (locations < p);
  for (int i = 0;i < n;i++) {
    if (option[i]) {
      Server[--p].location = i;
      Server[p].past_location = i;
    }
  }
  delete [] option;
}

SQM_solution::~SQM_solution () {
  delete [] X;
  for (int k = 0;k < m;k++) delete [] a[k];
  delete [] a;
}

void SQM_solution::set_server_location(int i,int j) {
  if ((i > 0) && (i < p)) {
    Server[i].past_location = Server[i].location;
    Server[i].location = j;
  }
}

int SQM_solution::get_server_location(int i) {
  if ((i > 0) && (i < p))
    return Server[i].location;
  return -1;
}

void SQM_solution:set_speed(double v,double beta) {
  for (int i = 0;i < p;i++) {
    Servers[i].v = v;
    Servers[i].beta = beta;
  }
}
