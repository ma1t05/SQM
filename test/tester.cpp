
#include "common.h"

command *Test_Function;

void Plot_Local_Search(SQM_instance &I,int p,double v);

int main(int argc,char *argv[]) {
  string filename;
  int M_clients,N_sites;
  int p,l;
  double v;
  stringstream LogName;
  SQM_instance *I;

  read_config_file("SQM.conf");
  filename = "Test";
  M_clients = 50; 
  N_sites = 30;
  p = 6; 
  l = 3;
  Mu_NT = MINS_PER_BLOCK * 3 / 60.0;
  lambda = MINS_PER_BLOCK * 6 / 60.0;
  v = 500.0;

  srand(time(NULL));

  /* Open Log File */
  I = SQM_load_instance(filename,M_clients,N_sites);
  //Test_exponential(*I,p,v);
  //Simulator(*I,p,v);
  Plot_Local_Search(*I,p,v);
  delete I;

}

void Test_MST(SQM_instance &I,int p,double v) {
  int loc;
  SQM_solution *X;
  double rt;
  cout << "Start: Test_MST" << endl;
  X = new SQM_solution(I,p);
  X->set_speed(v,BETA);
  X->set_params(lambda,Mu_NT);
  for (int i = 0;i < p;i++) {
    rt = X->get_response_time();
    cout << rt << endl;
    loc = X->get_server_location(0);
    X->remove_server(0);
    X->add_server();
    X->set_server_location(p-1,loc);
    X->set_speed(v,BETA);
    cout << X->get_response_time() << "\t" << X->get_response_time() - rt << endl;
  }
  delete X;
}

void Test_exponential(SQM_instance &I,int p,double v) {
  int m,*a;
  double demand;
  double lambda_j;
  double *times;

  m = I.demand_points();
  demand = I.total_demand();
  times = new double [m];
  for (int j = 0;j < m;j++) {
    lambda_j = lambda * I.get_demand(j) / demand;
    times[j] = exponential(lambda_j);
  }
  a = new int [m];
  sort_dist(m,times,a);

  cout << "Site\tRate\t\tTime" << endl;
  for (int j = 0;j < m;j++) {
    lambda_j = lambda * I.get_demand(j) / demand;
    cout << setfill('0') << setw(3) << a[j] << "\t" 
	 << setfill(' ')
	 << setw(10) << lambda_j << "\t" 
	 << setw(10) << times[a[j]] << endl;
  }
  
  delete [] a;
  delete [] times;

}

void Test_Path_Relinking(SQM_instance &I,int p,double v) {
  SQM_solution *X,*Y;
  int *match,*order;
  X = new SQM_solution(I,p);
  X->set_speed(v,BETA);
  Y = X->clone();
  SQM_heuristic(*Y);
  match = PR_perfect_matching(*X,*Y);
  order = PR_processing_order_nf(*X,match,*Y);
  X->update_preferred_servers ();
  gnuplot_Path_Relinking(*X,match,*Y,"step_00");
  for (int step = 0;step < p;step++) {
    int x = order[step];
    int y = match[x];
    int loc_x = X->get_server_location(x);
    int loc_y = Y->get_server_location(y);
    if (loc_x != loc_y) {
      X->set_server_location(x,loc_y);
      X->update_preferred_servers ();
      char step_file[8];
      sprintf(step_file,"step_%02d",step);
      gnuplot_Path_Relinking(*X,match,*Y,string(step_file));
    }
  }
  delete [] order;
  delete [] match;
  delete Y;
  delete X;
}

void Plot_Local_Search(SQM_instance &I,int p,double v) {
  SQM_solution *X,*Y;
  int *match,*order;
  X = new SQM_solution(I,p);
  X->set_speed(v,BETA);
  SQM_heuristic(*X);
  int m = I.demand_points();
  int n = I.potential_sites();
  Y = X->clone();
  int best_loc = UNASIGNED_LOCATION;
  double best_rt;
  /* Obtain the server with less workload */
  int j = LS_get_server_with_less_workload(*X);
  int loc_j  = X->get_server_location(j);
  /* obtain the news workloads */
  int k = LS_get_server_with_more_workload(*X);
  /* Put a server near the server with more workload */
  Sites *lst = LS_get_adjacent_sites(*X,k);

  if (best_loc == UNASIGNED_LOCATION) {
    cerr << "Warining: No new location" << endl;
    best_loc = loc_j;
  }
  X->test_server_location(j,loc_j); /* The past location of the server */
  X->set_server_location(j,best_loc);

  FILE *gnuPipe;
  fstream pointsfile;
  char demand_output[32],facility_output[32],
    centers_output[32],backup_output[32];
  sprintf(demand_output,"Tmp_demand.dat");
  gnuplot_write_points_file(demand_output,I.demand(0),m);
  sprintf(facility_output,"Tmp_facility.dat");
  gnuplot_write_points_file(facility_output,I.site(0),n);

  sprintf(centers_output,"Tmp_centers.dat");
  pointsfile.open(centers_output,fstream::out);
  for (int l = 0;l < n;l++) {
    int r = 0;
    for (int i = 0;i < p;i++) 
      if (X->get_server_location(i) == l)
	r++;
    if (r > 0)
      pointsfile << I.site(l)->x << " "
		  << I.site(l)->y << " "
		  << r << endl;
  }
  pointsfile.close();

  sprintf(backup_output,"Tmp_backup.dat");
  pointsfile.open(backup_output,fstream::out);
  for (Sites::iterator it = lst->begin(),end = lst->end(); it != end;it++) {
    int r = *it;
    pointsfile << I.site(r)->x << " "
	       << I.site(r)->y << " "
	       << r << endl;
  }
  pointsfile.close();

  gnuPipe = popen("gnuplot","w");
  gnuplot_sets(gnuPipe);
  fprintf(gnuPipe,"set output './plots/Local_Search.jpeg'\n");
  gnuplot_unsets(gnuPipe);

  fprintf(gnuPipe,"plot ");
  fprintf(gnuPipe,"'%s'",facility_output);
  /*fprintf(gnuPipe," using 1:2 w p lt 2 pt 10 title 'Facility'");*/
  fprintf(gnuPipe," using 1:2:(15) with circles lc rgb 'gray'");
  fprintf(gnuPipe," title 'Facility'");
  fprintf(gnuPipe,", '%s'",demand_output);
  /*fprintf(gnuPipe," using 1:2 w p lt 2 pt 7 lc rgb 'blue' title 'Demand'");*/
  fprintf(gnuPipe," using 1:2:(10) with circles lc rgb 'yellow'");
  fprintf(gnuPipe," title 'Demand'");
  fprintf(gnuPipe,", '%s'",centers_output);
  fprintf(gnuPipe," using 1:2:(8) with circles lc rgb 'dark-grey'");
  fprintf(gnuPipe," title 'Opened'");
  fprintf(gnuPipe,", '%s'",backup_output);
  fprintf(gnuPipe," using 1:2:(5) with circles lc rgb 'dark-green'");
  fprintf(gnuPipe," title 'Option'");
  fprintf(gnuPipe,"\n");
  fprintf(gnuPipe,"set object circle at %d,%d",I.site(loc_j)->x,I.site(loc_j)->y);
  fprintf(gnuPipe," size scr 0.1 fc rgb 'navy'\n");
  pclose(gnuPipe);
  delete lst;
  remove(facility_output);
  remove(demand_output);
  remove(centers_output);
  remove(backup_output);
}
