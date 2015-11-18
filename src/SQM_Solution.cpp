x
#include "SQM_Solution.h"

server::server () {
  location = 0;
  past_location = 0;
  v = 1.0;
  beta = 2.0;
}

server::server (int i) {
  location = i;
  past_location = i;
  v = 1.0;
  beta = 2.0;
}

void server::set_speed(double speed,double b) {
  v = speed;
  beta = b;
}

double server::get_speed () {
  return beta / v;
}

int server::get_location () {
  return location;
}

int server::get_past_location () {
  return past_location;
}

void server::set_location (int i) {
  past_location = location;
  location = i;
}

SQM_solution::SQM_solution (SQM_instance *I) {
  Inst = I;
  Servers = NULL;
  p = 0;
}

SQM_solution::SQM_solution (SQM_instance *I,int servers) {
  bool *option;
  int location,locations = 0;
  int l = 0;

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
    if (option[i]) 
      Servers[l++].set_location(i);
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
  for (int i = 0;i < p;i++)
    Servers[i].set_speed(v,beta);
}

double SM_solution::get_server_speed(int i) {
  return Servers[i].get_speed();
}

SQM_instance* SQM_solution::get_instance () {
  return Inst;
}

int SQM_solution::get_servers () {
  return p;
}
