
#include <iomanip>
#include <cmath>
#include "Simulation.h"
#include "log.h"

void Simulator_release_server(status &state,event &release);
void Simulator_attend_call(status &state,event &call);
bool server_is_free(status &state,int server);
event* get_nearest_event_from_queue(status &state,int server);

event_type event::get_type () {
  return type;
}

double event::get_time () {
  return at_time;
}

int event::get_node () {
  return node;
}

std::ostream& operator<<(std::ostream &os,event &incident) {
  double at_time = incident.get_time();
  os << "[" << setfill('0') << setw(2) << floor(at_time) + Start_Time
     << ":" << setfill('0') << setw(2) << floor(60 * (at_time - floor(at_time))) << "] ";
  switch (incident.get_type()) {
  case CALL: 
    os << "Incoming call from demand point " << incident.get_node();
    break;
  case RELEASE:
    os << "Server " << incident.get_node() << " released";
    break;
  default:
    os << "Unkow " << incident.get_node() << " event";
    break;
  }
  return os;
}

void Simulator(SQM_instance *I,int p,double lambda,double Mu_NT,double v) {
  int N = 500;
  SQM_solution *X = new SQM_solution(I,p);
  status state;
  event *incident;
  list<event*> releases;

  X->set_speed(v,BETA);
  X->set_params(lambda,Mu_NT);
  SQM_heuristic(X);
  state.Sol = X;
  state.busy = new bool[p];
  state.busy_time = new double[p];
  state.calls_sent_to_queue = 0;
  state.waiting_time = 0.0;

  for (int i = 0; i < p;i++) state.busy[i] = false;
  for (int i = 0; i < p;i++) state.busy_time[i] = 0.0;

  for (int i = 0;i < N;i++) {
    state.events = Generate_calls(I,lambda);
    state.current_time = 0;

    while (!state.events->empty()) {
      incident = state.events->front();
      state.events->pop_front();
      switch (incident->get_type()) {
      case CALL:
	Simulator_attend_call(state,*incident);
	break;
      case RELEASE:
	Simulator_release_server(state,*incident);
	break;
      default:
	cout << "Unknow Type!" << endl;
	break;
      }
    }
    delete state.events;
  }

  double *wl = X->get_workload();
  cout << "Workload\tBusy times" << endl;
  for (int i = 0; i < p;i++)
    cout << 10 * wl[i] << "\t" << state.busy_time[i]/N << endl;
  cout << endl;
  cout << " Calls send to queue: " << (double) state.calls_sent_to_queue / N << endl
       << "Average waiting time: " << state.waiting_time / state.calls_sent_to_queue << endl;
  delete [] state.busy;
  delete [] state.busy_time;
  delete X;
}

list<event*>* Generate_calls(SQM_instance *I,double lambda) {
  int m;
  int events;
  double current_time,etime;
  double demand,lambda_j;
  list<event*> *calls;
  list<event*>::iterator it;
  event *call;

  logDebug(cout << "Start Generate_calls" << endl);
  calls = new list<event*>;
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
      call = new event(CALL,current_time,j);
      /* Insert Call */
      logDebug(cout << "Insert event with lambda_j = " << lambda_j << endl);
      while (it != calls->end() && current_time > (*it)->get_time()) it++;
      /*
	if (it != calls->end())
	calls->insert(it,call);
	else calls->push_back(call);
      */
      calls->insert(it,call);
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

void Simulator_release_server(status &state,event &release) {
  int server = release.get_node();
  Log_Simulation << release << endl;
  state.busy[server] = false;
  state.current_time = release.get_time();
  delete &release;

  int p = state.Sol->get_servers();
  Log_Simulation << "\tBusy servers:";
  for (int i = 0;i < p;i++)
    if (state.busy[i]) Log_Simulation << " " << i;
  Log_Simulation << endl;

  if (!state.queue.empty()) {
    event *queued = NULL;

    Log_Simulation << "\tQueue isn't empty! [";
    queued = get_nearest_event_from_queue(state,server);
    /*
      queued = state.queue.front();
      state.queue.pop_front();
    */
    Log_Simulation << "]" << endl;
    state.waiting_time += state.current_time - queued->get_time();
    Simulator_attend_call(state,*queued);
  }
}

void Simulator_attend_call(status &state,event &call) {
  int p = state.Sol->get_servers();
  int **pref_srvrs = state.Sol->preferred_servers();
  int *ps = pref_srvrs[call.get_node()];
  double service_time;

  Log_Simulation << call << endl;
  if (call.get_type() == CALL)
    state.current_time = call.get_time();

  // Lookup for nearest free server
  for (int i = 0;i < p;i++) {
    if (server_is_free(state,ps[i])) {
      delete &call;
      Log_Simulation << "\tCall hosted by server " << ps[i] << endl;
      
      state.busy[ps[i]] = true;
      service_time = exponential(state.Sol->get_non_travel_time())
	+ state.Sol->distance(i,call.get_node()) * state.Sol->get_server_rate(i);
	
      event *release = new event(RELEASE,state.current_time + service_time,ps[i]);
      state.busy_time[ps[i]] += service_time; /* data */

      list<event*>::iterator it = state.events->begin();
      while (it != state.events->end() && release->get_time() > (*it)->get_time()) it++;
      state.events->insert(it,release);

      return;
    }
  }

  state.queue.push_back(&call);
  state.calls_sent_to_queue++;
  Log_Simulation << "\tCall send to queue:";
  for (list<event*>::iterator it = state.queue.begin();it != state.queue.end();it++){
    Log_Simulation << " " << (*it)->get_node();
  }
  Log_Simulation << endl;

}

bool server_is_free(status &state,int server) {
  return !state.busy[server];
}

event* get_nearest_event_from_queue(status &state,int server) {
  event *queued;
  list<event*>::iterator it,aux;

  queued = NULL;
  for (it = state.queue.begin();it != state.queue.end();it++) {
    Log_Simulation << " " << (*it)->get_node();
    if (queued == NULL || state.Sol->distance(server,(*it)->get_node()) < state.Sol->distance(server,queued->get_node())) {
      queued = *it;
      aux = it;
    }
  }
  state.queue.erase(aux);
  return queued;
}
