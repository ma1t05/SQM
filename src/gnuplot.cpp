
#include "gnuplot.h"

void gnu_sets(FILE*);
void gnu_unsets(FILE*);

void plot_instance_solution(SQM_instance *I,int *Sol,string output) {
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
  gnu_sets(gnuPipe);
  fprintf(gnuPipe,"set output '%s.svg'\n",output.c_str());
  gnu_unsets(gnuPipe);

  fprintf(gnuPipe,"plot ");
  fprintf(gnuPipe,"'%s' using 1:2 with points title 'Demand'",demand_output);
  fprintf(gnuPipe,", '%s' using 1:2 with points title 'Facility'",facility_output);
  fprintf(gnuPipe,", '%s' using 1:2:($3*0.1) with circles title 'Opened'",centers_output);
  fprintf(gnuPipe,"\n");
  pclose(gnuPipe);
  remove(demand_output);
  remove(facility_output);
  remove(centers_output);
}

void plot_solution_allocation(SQM_instance* I,int *Sol,string output) {
  
  FILE *gnuPipe = popen("gnuplot","w");
  gnu_sets(gnuPipe);
  fprintf(gnuPipe,"set output '%s.svg'\n",output.c_str());
  gnu_unsets(gnuPipe);
  
  pclose(gnuPipe);
}

void gnu_sets(FILE *gnuPipe) {
  fprintf(gnuPipe,"set term svg\n");
  /*fprintf(gnuPipe,"set key outside\n");*/
}

void gnu_unsets(FILE *gnuPipe) {
  fprintf(gnuPipe,"unset key\n");
  fprintf(gnuPipe,"unset border\n");
  fprintf(gnuPipe,"unset yzeroaxis\n");
  fprintf(gnuPipe,"unset xtics\n");
  fprintf(gnuPipe,"unset ytics\n");
  fprintf(gnuPipe,"unset ztics\n");
}
