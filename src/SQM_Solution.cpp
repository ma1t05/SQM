
#include "SQM_Solution.h"
#include "MST.h"
#include "random.h"

#define UNASIGNED_LOCATION -1

server::server () {
  location = 0;
  past_location = UNASIGNED_LOCATION;
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

double server::get_beta () {
  return beta;
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

void server::test_location (int i) {
  if (past_location == UNASIGNED_LOCATION) past_location = i;
  location = i;
}

SQM_solution::SQM_solution (SQM_instance *I) {
  Inst = I;
  p = 0;
  Servers = NULL;
  a = NULL;
  response_time = -1;
}

SQM_solution::SQM_solution (SQM_instance *I,int servers) {
  int n = I->potential_sites();
  bool *option;
  int location,locations = 0;
  int l = 0;

  Inst = I;
  p = servers;
  Servers = new server[p];
  a = NULL;
  response_time = -1;

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

SQM_solution::SQM_solution (SQM_solution *Sol) {

  Inst = Sol->get_instance();
  p = Sol->get_servers();
  Servers = new server[p];
  a = NULL;
  response_time = -1;

  for (int i = 0;i < p;i++) {
    Servers[i].set_location(Sol->get_server_past_location(i));
    Servers[i].set_location(Sol->get_server_location(i));
    Servers[i].set_speed(Sol->get_server_speed(i),Sol->get_server_beta(i));
  }
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

SQM_solution* SQM_solution::clone() {
  SQM_solution *clon;
  clon = new SQM_solution(this);
  clon->set_params(lambda,Mu_NT);
  return clon;
}

void SQM_solution::set_params(double Lambda,double Mu) {
  lambda = Lambda;
  Mu_NT = Mu;
}

void SQM_solution::set_server_location (int i,int j) {
  if ((i >= 0) && (i < p)) {
    Servers[i].set_location(j);
    response_time = -1;
  }
}

void SQM_solution::test_server_location (int i,int j) {
  if ((i >= 0) && (i < p)) {
    Servers[i].test_location(j);
    response_time = -1;
  }
}

int SQM_solution::get_server_location (int i) {
  if ((i >= 0) && (i < p))
    return Servers[i].get_location();
  return -1;
}

int SQM_solution::get_server_past_location (int i) {
  if ((i >= 0) && (i < p))
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

double SQM_solution::get_server_beta (int i) {
  return Servers[i].get_beta();
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

void SQM_solution::add_server () {
  server *aux;
  int k = p;
  aux = new server [++p];
  for (int i = 0;i < k;i++)
    aux[i].set_location(Servers[i].get_location());
  if (Servers != NULL) delete [] Servers;
  Servers = aux;
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

double SQM_solution::get_arrival_rate () {
  return lambda;
}

double SQM_solution::get_non_travel_time () {
  return Mu_NT;
}

double SQM_solution::get_response_time() {
  if (response_time == -1)
    response_time = MST_response_time(this);
  return response_time;
}

int** SQM_solution::preferred_servers () {
  return a;
}

bool SQM_solution::operator>(SQM_solution& X) {
  return get_response_time() > X.get_response_time();
}

bool SQM_solution::operator<(SQM_solution& X) {
  return X > (*this);
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

