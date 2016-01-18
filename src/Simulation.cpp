
#include <iomanip>
#include <cmath>
#include "Simulation.h"
#include "log.h"

void print_time(double current_time);
void insert_event(list_events&,event&);
void Simulator_release_server(status &state,event &release);
void Simulator_attend_call(status &state,event &call);
bool server_is_free(status &state,int server);
call* get_nearest_event_from_queue(status &state,int server);

void print_time(std::ostream &os,double current_time) {
  int hour = floor(current_time);
  int minute = floor(60*(current_time - hour));
  os << "[" << setfill('0') << setw(2) << hour
     << ":" << setfill('0') << setw(2) << minute
     << "]";
}

void insert_event(list_events &events,event &x) {
  list_events::iterator it = events.begin(),end = events.end();
  while (it != end && x.get_time() > (*it)->get_time()) it++;
  events.insert(it,&x);
}

/* Class 'call' methods */
call::call (double t,int i)  : at_time(t), demand_point(i) {
  queued = false;
  queued_at_time = NOT_QUEUED;
}

void call::process (status &state) {
  int p = state.Sol->get_servers();
  int **pref_srvrs = state.Sol->preferred_servers();
  int *ps = pref_srvrs[demand_point];
  double service_time;

  Log_Simulation << *this << endl;
  state.current_time = at_time;
  state.total_calls++;

  /* Lookup for nearest free server */
  for (int i = 0;i < p;i++) {
    if (server_is_free(state,ps[i])) {
      new release (state,ps[i],demand_point);
      delete this;
      return;
    }
  }

  queue_call(state);
}

int call::get_point () const {
  return demand_point;
}

double call::get_time () const{
  return at_time;
}

void call::queue_call (status &state) {
  queued = true;
  queued_at_time = at_time;

  state.queue.push_back(this);
  state.calls_sent_to_queue++;

  Log_Simulation << "\tCall send to queue:";
  for (queued_calls::iterator it = state.queue.begin(),end = state.queue.end();
       it != end;
       it++)
    Log_Simulation << " " << (*it)->get_point();
  Log_Simulation << endl;
}

void call::dequeue_call (status &state,int server) {
  set_time(state.current_time);
  state.waiting_time += get_waiting_time();
  Log_Simulation << *this << endl;
  new release (state,server,demand_point);
  delete this;
}

void call::set_time (double t) {
  at_time = t;
}

bool call::is_queued () const {
  return queued;
}

double call::get_waiting_time () const {
  return at_time - queued_at_time;
}

std::ostream& operator<<(std::ostream &os,call &incident) {
  print_time(os,incident.get_time() + Start_Time);

  if (incident.is_queued()) {
    os << " Queued call since ";
    print_time(os,incident.get_waiting_time());
  }
  else 
    os << " Incoming call";
  os << " from demand point " << incident.get_point();
    
  return os;
}

/* Class 'release' methods */
release::release(status &state,int id_server,int id_demand) {
  double distance = state.Sol->distance(id_server,id_demand);
  server = id_server;
  demand_point = id_demand;
  travel_time = distance / state.Sol->get_server_speed(server);
  on_scene_st = exponential(state.Sol->get_non_travel_time());
  follow_up_travel_time =
    + distance * state.Sol->get_server_rate(server) 
    - travel_time;
  service_time = travel_time + on_scene_st + follow_up_travel_time;
  return_time = state.current_time + service_time;
  insert_event(*(state.events),*this);
  Log_Simulation << "\tCall hosted by server " << id_server << endl;

  /* data */
  state.busy[server] = true;
  state.busy_time[server] += service_time;
  state.arrival_time += travel_time;
  state.service_time += service_time;
}

void release::process (status &state) {
  state.busy[server] = false;
  state.current_time = get_time();

  Log_Simulation << *this << endl;
  int p = state.Sol->get_servers();
  Log_Simulation << "\tBusy servers:";
  for (int i = 0;i < p;i++)
    if (state.busy[i]) Log_Simulation << " " << i;
  Log_Simulation << endl;

  if (!state.queue.empty()) {
    call *queued_call = NULL;
    queued_call = get_nearest_event_from_queue(state,server);
    queued_call->dequeue_call(state,server);
  }
  delete this;
}

int release::get_point () const {
  return server;
}

double release::get_time () const {
  return return_time;
}

double release::get_travel_time () const {
  return travel_time;
}

std::ostream& operator<<(std::ostream &os,release &incident) {
  print_time(os,incident.get_time() + Start_Time);
  os << " Server " << incident.get_point() << " released";
  return os;
}

void Simulator(SQM_instance *I,int p,double lambda,double Mu_NT,double v) {
  int N = 500;
  SQM_solution *X;
  status state;
  event *incident;

  X = new SQM_solution(I,p);
  X->set_speed(v,BETA);
  X->set_params(lambda,Mu_NT);
  logLevel = LOG_DEBUG;
  SQM_heuristic(X);
  logLevel = LOG_INFO;

  state.Sol = X;
  state.busy = new bool[p];
  state.busy_time = new double[p];
  state.total_calls = 0;
  state.calls_sent_to_queue = 0;
  state.waiting_time = 0.0;
  state.arrival_time = 0.0;
  state.service_time = 0.0;

  for (int i = 0; i < p;i++) state.busy[i] = false;
  for (int i = 0; i < p;i++) state.busy_time[i] = 0.0;

  for (int i = 0;i < N;i++) {
    state.events = Generate_calls(I,lambda);
    state.current_time = 0;

    while (!state.events->empty()) {
      incident = state.events->front();
      state.events->pop_front();
      incident->process(state);
    }
    delete state.events;
  }

  double *wl = X->get_workload();
  cout << "Workload\tBusy times" << endl;
  for (int i = 0; i < p;i++)
    cout << 10 * wl[i] << "\t" << state.busy_time[i]/ (N * Simulation_Time) << endl;
  cout << endl;
  delete [] wl;
  cout << "       Attended calls: " 
       << (double) state.total_calls / N << endl
       << "  Calls send to queue: " 
       << (double) state.calls_sent_to_queue / N << endl
       << "\t*** Average times ***" << endl
       << "service\twaiting\ttravel\tarrival" << endl;
  print_time(cout,state.service_time / state.total_calls); /* service */
  cout << "\t";
  print_time(cout,state.waiting_time / state.calls_sent_to_queue); /* waiting */
  cout << "\t";
  print_time(cout,state.arrival_time / state.total_calls); /* travel */
  cout << "\t";
  print_time(cout,(state.arrival_time + state.waiting_time) / state.total_calls); /* arrival */
  cout << endl;

  delete [] state.busy;
  delete [] state.busy_time;
  delete X;
}

list_events* Generate_calls(SQM_instance *I,double lambda) {
  int m;
  int events;
  double current_time,etime;
  double demand,lambda_j;
  list_events *calls;
  list_events::iterator it;
  call *incoming_call;

  logDebug(cout << "Start Generate_calls" << endl);
  calls = new list_events;
  m = I->demand_points();
  demand = I->total_demand();
  for (int j = 0;j < m;j++) {
    lambda_j = lambda * I->get_demand(j) / demand;
    current_time = exponential(lambda_j);
    it = calls->begin();
    events = 0;
    while (current_time < Simulation_Time) {
      /* Create Call */
      logDebug(cout << "Create event" << endl);
      incoming_call = new call(current_time,j);
      /* Insert Call */
      logDebug(cout << "Insert event with lambda_j = " << lambda_j << endl);
      while (it != calls->end() && current_time > (*it)->get_time()) it++;
      calls->insert(it,incoming_call);
      /* Determine time of next Call */
      etime = exponential(lambda_j);
      current_time += etime;
      events++;
      logDebug(cout << "Determine time of next event " << current_time << endl);
    }
    logDebug(cout << "Create " << events << " for demand point " << j << endl);
  }
  
  return calls;
}

bool server_is_free(status &state,int server) {
  return !state.busy[server];
}

call* get_nearest_event_from_queue(status &state,int server) {
  double distance;
  call *queued;
  queued_calls::iterator it,aux;

  Log_Simulation << "\tQueue isn't empty! [";
  for (it = state.queue.begin();it != state.queue.end();it++)
    Log_Simulation << " " << (*it)->get_point();
  Log_Simulation << " ]" << endl;

  queued = NULL;
  for (it = state.queue.begin();it != state.queue.end();it++) {
    if (queued == NULL ||
	state.Sol->distance(server,(*it)->get_point()) < distance) {
      queued = *it;
      distance = state.Sol->distance(server,queued->get_point());
      aux = it;
    }
  }
  state.queue.erase(aux);
  return queued;
}

event* get_first_event_from_queue(status &state,int server) {
  event *queued;
  queued = state.queue.front();
  state.queue.pop_front();
  return queued;
}
