
#include "instance-creator.h"

float unif(float,float);

SQM_instance* IC_create_instance(int m,int n) {
  SQM_instance *I;
  I = new SQM_instance;
  I->M = m;
  I->N = n;
  I->V = new point[m];
  I->W = new point[n];
  for (int i = 0;i < m;i++) {
    (I->V)[i].x = unif(MIN_X,MAX_X);
    (I->V)[i].y = unif(MIN_Y,MAX_Y);
    (I->V)[i].demand = unif(32,1024);
  }
  for (int i = 0;i < n;i++) {
    (I->W)[i].x = unif(MIN_X,MAX_X);
    (I->W)[i].y = unif(MIN_Y,MAX_Y);
    (I->W)[i].demand = 0.0;
  }
  return I;
}

SQM_instance* IC_read_instance (string Demand_nodes,string facility_nodes) {
  int n,m;
  fstream demandfile,facilityfile;
  SQM_instance *I;
  I = new SQM_instance;
  demandfile.open(Demand_nodes.c_str(),fstream::in);
  demandfile >> m;
  I->M = m;
  I->V = new point[m];
  for (int i = 0;i < m;i++) {
    demandfile >> (I->V)[i].x >> (I->V)[i].y >> (I->V)[i].demand;
  }
  demandfile.close();
  
  facilityfile.open(facility_nodes.c_str(),fstream::in);
  facilityfile >> n;
  I->N = n;
  I->W = new point[n];
  for (int j = 0;j < n;j++) {
    facilityfile >> (I->W)[j].x >> (I->W)[j].y;
  }
  facilityfile.close();
  return I;
}

SQM_instance* IC_load_instance (string Demand_nodes) {
  int m;
  double S;
  fstream demandfile;
  SQM_instance *I;
  I = new SQM_instance;
  demandfile.open(Demand_nodes.c_str(),fstream::in);
  demandfile >> m >> S;
  I->M = m;
  I->V = new point[m];
  for (int i = 0;i < m;i++) {
    demandfile >> (I->V)[i].x >> (I->V)[i].y >> (I->V)[i].demand;
  }
  demandfile.close();
  
  I->N = m;
  I->W = new point[m];
  for (int j = 0;j < m;j++) {
    (I->W)[j].x = (I->V)[j].x;
    (I->W)[j].y = (I->V)[j].y;
  }

  return I;
}

void IC_write_instance(SQM_instance* I,string demand_output,string facility_output) {
  fstream demandfile,facilityfile;
  
  demandfile.open(demand_output.c_str(),fstream::out);
  demandfile << I->M << endl;
  for (int i = 0;i < I->M;i++) {
    demandfile << (I->V)[i].x << " " 
	       << (I->V)[i].y << " "
	       << (I->V)[i].demand 
	       << endl;
  }
  demandfile.close();

  facilityfile.open(facility_output.c_str(),fstream::out);
  facilityfile << I->N << endl;
  for (int i = 0;i < I->N;i++) {
    facilityfile << (I->W)[i].x << " " 
		 << (I->W)[i].y
		 << endl;
  }
  facilityfile.close();
}

float unif(float a,float b) {
  return a + (b - a) * rand() / RAND_MAX;
}

void IC_delete_instance(SQM_instance* I) {
  delete [] I->V;
  delete [] I->W;
  delete I;
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
