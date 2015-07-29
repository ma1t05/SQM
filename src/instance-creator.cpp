
#include "instance-creator.h"

float unif(float,float);

SQM_instance* IC_create_instance(int n,int m) {
  SQM_instance *I;
  I = new SQM_instance;
  I->M = m;
  I->N = n;
  I->V = new point[m];
  I->W = new point[n];
  for (int i = 0;i < m;i++) {
    (I->V)[i].x = unif(MIN_X,MAX_X);
    (I->V)[i].y = unif(MIN_Y,MAX_Y);
  }
  for (int i = 0;i < n;i++) {
    (I->W)[i].x = unif(MIN_X,MAX_X);
    (I->W)[i].y = unif(MIN_Y,MAX_Y);
  }
  return I;
}

SQM_instance* IC_read_instance (string Demand_nodes,string facility_nodes) {
  SQM_instance *I;
  I = new SQM_instance;
  /* Pendiente */
  return I;
}

void IC_write_instance(SQM_instance* I,string output) {
  
}

float unif(float a,float b) {
  return a + (b - a) * rand() / RAND_MAX;
}
