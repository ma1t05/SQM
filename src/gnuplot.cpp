
#include "gnuplot.h"

void gnuplot_sets(FILE*);
void gnuplot_unsets(FILE*);
void gnuplot_write_points_file(char*,point*,int);

void plot_instance_solution(SQM_instance *I,int *Sol,string output) {
  int n = I->N,m = I->M;
  char demand_output[32],facility_output[32],centers_output[32];
  fstream centersfile;
  
  sprintf(demand_output,"Tmp_demand_%d.dat",rand());
  gnuplot_write_points_file(demand_output,I->V,m);

  sprintf(facility_output,"Tmp_facility_%d.dat",rand());
  gnuplot_write_points_file(facility_output,I->W,m);

  sprintf(centers_output,"Tmp_centers_%d.dat",rand());
  centersfile.open(centers_output,fstream::out);
  for (int j = 0;j < n;j++) {
    if (Sol[j] > 0)
      centersfile << (I->W)[j].x << " "
		  << (I->W)[j].y << " "
		  << Sol[j] << endl;
  }

  FILE *gnuPipe = popen("gnuplot","w");
  gnuplot_sets(gnuPipe);
  fprintf(gnuPipe,"set output '%s.svg'\n",output.c_str());
  gnuplot_unsets(gnuPipe);

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

void plot_solution_allocation(SQM_instance* I,int p,response_unit *X,double **f,string output,string suffix) {
  int j,r;
  int edge_key;
  fstream centersfile,edges_file;
  char demand_output[32],facility_output[32],centers_output[32],edges_output[32];
  
  sprintf(demand_output,"Tmp_demand_%d.dat",rand());
  gnuplot_write_points_file(demand_output,I->V,I->M);

  sprintf(facility_output,"Tmp_facility_%d.dat",rand());
  gnuplot_write_points_file(facility_output,I->W,I->N);

  sprintf(centers_output,"Tmp_centers_%d.dat",rand());
  centersfile.open(centers_output,fstream::out);
  for (int j = 0;j < I->N;j++) {
    r = 0;
    for (int i = 0;i < p;i++) 
      if (j == X[i].location)
	r++;
    if (r > 0)
      centersfile << (I->W)[j].x << " "
		  << (I->W)[j].y << " "
		  << r << endl;
  }
  centersfile.close();

  edge_key = rand();
  for (int i = 0;i < p;i++) {
    sprintf(edges_output,"Tmp_edges_center_%d_%d.dat",i+1,edge_key);
    edges_file.open(edges_output,fstream::out);
    j = X[i].location;
    for (int k = 0;k < I->M;k++) {
      if (100 * f[i][j] > 1)
	edges_file << (I->W)[j].x << " " << (I->W)[j].y << " "
		   << (I->V)[k].x << " " << (I->V)[k].y << " "
		   << ceil(sqrt(f[i][k] * 100)) << endl;
    }
    edges_file.close();
  }

  FILE *gnuPipe = popen("gnuplot","w");
  gnuplot_sets(gnuPipe);
  gnuplot_unsets(gnuPipe);
  fprintf(gnuPipe,"set for [i=1:101] style arrow i lw i/10.0 nohead\n");

  fprintf(gnuPipe,"set output '%s_%s.svg'\n",output.c_str(),suffix.c_str());
  fprintf(gnuPipe,"plot ");
  fprintf(gnuPipe,"'%s' using 1:2 with points title 'Demand'",demand_output);
  fprintf(gnuPipe,", '%s' using 1:2 with points title 'Facility'",facility_output);
  fprintf(gnuPipe,", '%s' using 1:2:($3*0.1) with circles title 'Opened'",centers_output);
  for (int i = 0;i < p;i++) {
    sprintf(edges_output,"Tmp_edges_center_%d_%d.dat",i+1,edge_key);
    fprintf(gnuPipe,", '%s' using 1:2:($3-$1):($4-$2):5 with vectors arrowstyle variable",edges_output);
  }
  fprintf(gnuPipe,"\n");

  for (int i = 0;i < p;i++) {
    sprintf(edges_output,"Tmp_edges_center_%d_%d.dat",i+1,edge_key);
    fprintf(gnuPipe,"set output '%s_center_%02d_%s.svg'\n",output.c_str(),i+1,suffix.c_str()); 
    fprintf(gnuPipe,"plot ");
    fprintf(gnuPipe,"'%s' using 1:2 with points title 'Demand'",demand_output);
    fprintf(gnuPipe,", '%s' using 1:2:($3*0.1) with circles title 'Opened'",centers_output);
    fprintf(gnuPipe,", '%s' using 1:2:($3-$1):($4-$2):5 with vectors arrowstyle variable",edges_output);
    fprintf(gnuPipe,"\n");
  }

  pclose(gnuPipe);
  for (int i = 0;i < p;i++) {
    sprintf(edges_output,"Tmp_edges_center_%d_%d.dat",i+1,edge_key);
    remove(edges_output);
  }
  remove(demand_output);
  remove(facility_output);
  remove(centers_output);
}

void gnuplot_sets(FILE *gnuPipe) {
  fprintf(gnuPipe,"set term svg\n");
  /*fprintf(gnuPipe,"set key outside\n");*/
}

void gnuplot_unsets(FILE *gnuPipe) {
  fprintf(gnuPipe,"unset key\n");
  fprintf(gnuPipe,"unset border\n");
  fprintf(gnuPipe,"unset yzeroaxis\n");
  fprintf(gnuPipe,"unset xtics\n");
  fprintf(gnuPipe,"unset ytics\n");
  fprintf(gnuPipe,"unset ztics\n");
}

void gnuplot_write_points_file(char *output,point *Set,int k) {
  fstream outputfile;
  outputfile.open(output,fstream::out);
  for (int i = 0;i < k;i++) {
    outputfile << Set[i].x << " " 
	       << Set[i].y << endl;
  }
  outputfile.close();
}
