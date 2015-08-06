
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
  SQM_instance *I;
  I = new SQM_instance;
  /* Pendiente */
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
  for (int i = 0;i < I->M;i++) {
    facilityfile << (I->W)[i].x << " " 
		 << (I->W)[i].y
		 << endl;
  }
  facilityfile.close();
}

void IC_plot_instance(string output,string demand_output,string facility_output) {
  FILE *gnuPipe = popen("gnuplot","w");
  fprintf(gnuPipe,"set term svg\n");
  fprintf(gnuPipe,"set output '%s.svg'\n",output.c_str());
  //fprintf(gnuPipe,".svg'\n");
  fprintf(gnuPipe,"unset key\n");
  fprintf(gnuPipe,"unset border\n");
  fprintf(gnuPipe,"unset yzeroaxis\n");
  fprintf(gnuPipe,"unset xtics\n");
  fprintf(gnuPipe,"unset ytics\n");
  fprintf(gnuPipe,"unset ztics\n");

  fprintf(gnuPipe,"plot ");
  fprintf(gnuPipe,"'%s' every ::1 using 1:2 with points",demand_output.c_str());
  fprintf(gnuPipe,"'%s' every ::1 using 1:2 with points",facility_output.c_str());
  fprintf(gnuPipe,"\n");
  pclose(gnuPipe);
}

float unif(float a,float b) {
  return a + (b - a) * rand() / RAND_MAX;
}
