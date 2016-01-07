
#include "Simulation.h"
#include "log.h"

void Simulator_release_server(status &state,event *rel);
void Simulator_attend_call(status &state,event *event);
bool server_is_free(status &state,int server);

void Simulator(SQM_instance *I,int p,double lambda,double Mu_NT,double v) {
  int N = 100;
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

  for (int i = 0; i < p;i++) state.busy[i] = false;
  for (int i = 0; i < p;i++) state.busy_time[i] = 0.0;

  for (int i = 0;i < N;i++) {
    state.events = Generate_calls(I,lambda);
    state.current_time = 0;

    while (!state.events->empty()) {
      incident = state.events->front();
      state.events->pop_front();
      switch (incident->type) {
      case CALL:
	Simulator_attend_call(state,incident);
	break;
      case RELEASE:
	Simulator_release_server(state,incident);
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

  delete [] state.busy;
  delete [] state.busy_time;
  delete X;
}

list<event*>* Generate_calls(SQM_instance *I,double lambda) {
  int m;
  int events;
  double Simulation_Time = 10;
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
      call = new event;
      call->type = CALL;
      call->at_time = current_time;
      call->node = j;
      /* Insert Call */
      logDebug(cout << "Insert event with lambda_j = " << lambda_j << endl);
      while (it != calls->end() && current_time > (*it)->at_time) it++;
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

void Simulator_release_server(status &state,event *rel) {
  Log_Simulation << "["<< rel->at_time << "] "
		 << "Server " << rel->node << " released" << endl;
  state.busy[rel->node] = false;
  state.current_time = rel->at_time;
  delete rel;
  int p = state.Sol->get_servers();
  Log_Simulation << "\tBusy servers:";
  for (int i = 0;i < p;i++)
    if (state.busy[i]) Log_Simulation << " " << i;
  Log_Simulation << endl;

  if (!state.queue.empty()) {
    Log_Simulation << "\tQueue isn't empty!" << endl;
    event *queued;
    // Determinar evento a atender
    queued = state.queue.front();
    state.queue.pop_front();
    // llamar Simulator_attend_call
    Simulator_attend_call(state,queued);
  }
}

void Simulator_attend_call(status &state,event *call) {
  int p = state.Sol->get_servers();
  int **pref_srvrs = state.Sol->preferred_servers();
  int *ps = pref_srvrs[call->node];

  if (call->type == CALL) {
    Log_Simulation << "[" << call->at_time << "] "
		   << "Incoming call from demand point " << call->node << endl;
    state.current_time = call->at_time;
  }
  else {
    Log_Simulation << "[" << state.current_time << "] "
		   << "Queued call from demand point " << call->node << endl;
  }

  // Lookup for nearest free server
  for (int i = 0;i < p;i++) {
    if (server_is_free(state,ps[i])) {
      delete call;
      Log_Simulation << "\tCall hosted by server " << ps[i] << endl;
      
      state.busy[ps[i]] = true;
      event *release = new event;
      release->type = RELEASE;
      release->node = ps[i];
      release->at_time = state.current_time 
	+ state.Sol->distance(i,call->node) * state.Sol->get_server_rate(i)
	+ exponential(state.Sol->get_non_travel_time());
      state.busy_time[ps[i]] += release->at_time - state.current_time; /* data */

      list<event*>::iterator it = state.events->begin();
      while (it != state.events->end() && release->at_time > (*it)->at_time) it++;
      state.events->insert(it,release);

      return;
    }
  }

  call->type = QUEUING;
  state.queue.push_back(call);
  Log_Simulation << "\tCall send to queue:";
  for (list<event*>::iterator it = state.queue.begin();it != state.queue.end();it++){
    Log_Simulation << " " << (*it)->node;
  }
  Log_Simulation << endl;

}

bool server_is_free(status &state,int server) {
  return !state.busy[server];
}
