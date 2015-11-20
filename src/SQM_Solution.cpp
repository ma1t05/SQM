
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

server::~server () {
  
}

void server::set_speed(double speed,double b) {
  v = speed;
  beta = b;
}

double server::get_speed () {
  return v;
}

double server::get_rate () {
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
  p = 0;
  Servers = NULL;
  a = NULL;
}

SQM_solution::SQM_solution (SQM_instance *I,int servers) {
  int n = I->potential_sites();
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
  if (Servers != NULL)
    delete [] Servers;
  if (a != NULL) {
    int M = Inst->demand_points();
    for (int k = 0;k < M;k++) delete [] a[k];
    delete [] a;
  }
}

void SQM_solution::set_server_location (int i,int j) {
  if ((i > 0) && (i < p))
    Servers[i].set_location(j);
}

int SQM_solution::get_server_location (int i) {
  if ((i > 0) && (i < p))
    return Servers[i].get_location();
  return -1;
}

int SQM_solution::get_server_past_location (int i) {
  if ((i > 0) && (i < p))
    return Servers[i].get_past_location();
  return -1;
}

void SQM_solution::set_speed (double v,double beta) {
  for (int i = 0;i < p;i++)
    Servers[i].set_speed(v,beta);
}

double SQM_solution::get_server_speed (int i) {
  return Servers[i].get_speed();
}

double SQM_solution::get_server_rate (int i) {
  return Servers[i].get_rate();
}

SQM_instance* SQM_solution::get_instance () {
  return Inst;
}

int SQM_solution::get_servers () {
  return p;
}

void SQM_solution::update_preferred_servers () {
  if (p == 0) return;
  double *d;
  int m = Inst->demand_points();
  if (a == NULL) {
    a = new int*[m];
    for (int k = 0;k < m;k++)
      a[k] = new int[p];
  }
  d = new double[p];
  for (int k = 0;k < m;k++) {
    for (int i = 0;i < p;i++)
      d[i] = Inst->distance(Servers[i].get_location(),k);
    sort_dist(p,d,a[k]);
  }
  delete [] d;
}

int** SQM_solution::preferred_servers () {
  return a;
}

// response_unit* guess_a_location_01(int p,int n, point *W){
//   response_unit *X;
//   X = new response_unit[p];
//   for (int i = 0;i < p;i++) {
//     Sol->get_server_location(i) = i;
//     X[i].past_location = i;
//   }
//   return X;
// }

// response_unit* guess_a_location_02(int p,int n, point *W){
//   response_unit *X;
//   X = new response_unit[p];
//   for (int i = 0;i < p;i++) {
//     Sol->get_server_location(i) = unif(n);
//   }
//   return X;
// }

// response_unit* guess_a_location_03(int p,int n, point *W){
//   response_unit *X;
//   bool *option;
//   int location,locations = 0;
//   X = new response_unit[p];
//   option = new bool[n];
//   for (int i = 0;i < n;i++) option[i] = false;
//   do {
//     location = unif(n);
//     logDebug(cout << "location = " << location << endl);
//     if (option[location] == false) {
//       locations++;
//       option[location] = true;
//     }
//   } while (locations < p);
//   for (int i = 0;i < n;i++) {
//     if (option[i]) {
//       X[--p].location = i;
//       X[p].past_location = i;
//     }
//   }
//   logDebug(cout << endl);
//   delete [] option;
//   return X;
// }

