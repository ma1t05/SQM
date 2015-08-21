
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

void IC_plot_instance(SQM_instance *I,int *Sol,string output) {
  int n = I->N,m = I->M;
  fstream demandfile,facilityfile,centersfile;
  char demand_output[32],facility_output[32],centers_output[32];
  
  sprintf(demand_output,"Tmp_demand_%d.dat",rand());
  demandfile.open(demand_output,fstream::out);
  
  for (int i = 0;i < m;i++) {
    demandfile << (I->V)[i].x << " " 
	       << (I->V)[i].y << endl;
  }
  demandfile.close();

  sprintf(facility_output,"Tmp_facility_%d.dat",rand());
  facilityfile.open(facility_output,fstream::out);
  for (int j = 0;j < n;j++) {
    facilityfile << (I->W)[j].x << " "
		 << (I->W)[j].y << endl;
  }
  facilityfile.close();

  sprintf(centers_output,"Tmp_centers_%d.dat",rand());
  centersfile.open(centers_output,fstream::out);
  for (int j = 0;j < n;j++) {
    if (Sol[j] > 0)
      centersfile << (I->W)[j].x << " "
		  << (I->W)[j].y << " "
		  << Sol[j] << endl;
  }

  FILE *gnuPipe = popen("gnuplot","w");
  fprintf(gnuPipe,"set term svg\n");
  fprintf(gnuPipe,"set output '%s.svg'\n",output.c_str());
  /*fprintf(gnuPipe,"set key outside\n");*/
  fprintf(gnuPipe,"unset key\n");
  fprintf(gnuPipe,"unset border\n");
  fprintf(gnuPipe,"unset yzeroaxis\n");
  fprintf(gnuPipe,"unset xtics\n");
  fprintf(gnuPipe,"unset ytics\n");
  fprintf(gnuPipe,"unset ztics\n");

  fprintf(gnuPipe,"plot ");
  fprintf(gnuPipe,"'%s' using 1:2 with points title 'Demand'",demand_output);
  fprintf(gnuPipe,", '%s' using 1:2 with points title 'Facility'",facility_output);
  fprintf(gnuPipe,", '%s' using 1:2:($3*10) with circles title 'Opened'",centers_output);
  fprintf(gnuPipe,"\n");
  pclose(gnuPipe);
  remove(demand_output);
  remove(facility_output);
  remove(centers_output);
}

float unif(float a,float b) {
  return a + (b - a) * rand() / RAND_MAX;
}
