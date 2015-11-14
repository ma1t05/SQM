 
#include <fstream>
#include <cstdlib>
#include "instance-creator.h"

float unif(float,float);

SQM_instance(string);
SQM_instance(string,string);

SQM_instance::SQM_instance (int m/* demand points */,int n/* sites */) {
  /* Create random demand points */
  M = m;
  V = new point[m];
  for (int i = 0;i < m;i++) {
    V[i].x = unif(MIN_X,MAX_X);
    V[i].y = unif(MIN_Y,MAX_Y);
    V[i].demand = unif(32,1024);
  }

  /* Create random potential facility sites */
  N = n;
  W = new point[n];
  for (int i = 0;i < n;i++) {
    W[i].x = unif(MIN_X,MAX_X);
    W[i].y = unif(MIN_Y,MAX_Y);
    W[i].demand = 0.0;
  }
  set_distances();
  set_preferred_servers();
}

SQM_instance::SQM_instance (string nodes) {
  double S;
  fstream demandfile;

  demandfile.open(Demand_nodes.c_str(),fstream::in);
  demandfile >> M >> S;

  /* Read demand points */
  V = new point[m];
  for (int i = 0;i < m;i++)
    demandfile >> V[i].x >> V[i].y >> V[i].demand;
  demandfile.close();
  
  /* Copy demand points to potential facility points */
  N = M;
  W = new point[N];
  for (int j = 0;j < N;j++) {
    W[j].x = V[j].x;
    W[j].y = V[j].y;
  }

  set_distances();
  set_preferred_servers();
}

SQM_instance::SQM_instance (string Demand_nodes,string facility_nodes) {
  fstream demandfile,facilityfile;

  /* Read demand points */
  demandfile.open(Demand_nodes.c_str(),fstream::in);
  demandfile >> M;
  V = new point[M];
  for (int i = 0;i < M;i++) {
    demandfile >> V[i].x >> V[i].y >> V[i].demand;
  }
  demandfile.close();
  
  /* Read potenctial location points */
  facilityfile.open(facility_nodes.c_str(),fstream::in);
  facilityfile >> N;
  W = new point[N];
  for (int j = 0;j < N;j++) {
    facilityfile >> W[j].x >> W[j].y;
  }
  facilityfile.close();

  set_distances();
  set_preferred_servers();
}

SQM_instance::~SQM_instance () {
  delete [] V;
  delete [] W;
}

void IC_write_instance(SQM_instance* I,string demand_output,string facility_output) {
  fstream demandfile,facilityfile;
  
  demandfile.open(demand_output.c_str(),fstream::out);
  demandfile << I->M << endl;
  for (int i = 0;i < I->M;i++) {
    demandfile << V[i].x << " " 
	       << V[i].y << " "
	       << V[i].demand 
	       << endl;
  }
  demandfile.close();

  facilityfile.open(facility_output.c_str(),fstream::out);
  facilityfile << I->N << endl;
  for (int i = 0;i < I->N;i++) {
    facilityfile << W[i].x << " " 
		 << W[i].y
		 << endl;
  }
  facilityfile.close();
}

float unif(float a,float b) {
  return a + (b - a) * rand() / RAND_MAX;
}

void IC_delete_instance(SQM_instance* I) {
}

double** SQM_dist_matrix(SQM_instance *I) {
  int n = I->N,m = I->M;
  double **Dist;

  Dist = new  double*[m];
  for (int j = 0;j < m;j++)
    Dist[j] = new double [n];

  for (int j = 0;j < m;j++) {
    for (int i = 0;i < n;i++) {
      Dist[j][i] = dist(&(I->V[j]),&(I->W[i]));
    }
  }
  return Dist;
}
