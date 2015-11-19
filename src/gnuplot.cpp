
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "gnuplot.h"

using namespace std;

void gnuplot_sets(FILE*);
void gnuplot_unsets(FILE*);
void gnuplot_write_points_file(char*,point*,int);

void plot_instance_solution(SQM_instance *I,int *Sol,string output) {
  int n = I->demand_points(),m = I->potential_sites();
  char demand_output[32],facility_output[32],centers_output[32];
  fstream centersfile;
  
  sprintf(demand_output,"Tmp_demand_%d.dat",rand());
  gnuplot_write_points_file(demand_output,I->demand(0),m);

  sprintf(facility_output,"Tmp_facility_%d.dat",rand());
  gnuplot_write_points_file(facility_output,I->site(0),m);

  sprintf(centers_output,"Tmp_centers_%d.dat",rand());
  centersfile.open(centers_output,fstream::out);
  for (int j = 0;j < n;j++) {
    if (Sol[j] > 0) {
      centersfile << I->site(j)->x << " "
		  << I->site(j)->y << " "
		  << Sol[j] << endl;
    }
  }

  FILE *gnuPipe = popen("gnuplot","w");
  gnuplot_sets(gnuPipe);
  fprintf(gnuPipe,"set output '%s.svg'\n",output.c_str());
  gnuplot_unsets(gnuPipe);

  fprintf(gnuPipe,"plot ");
  fprintf(gnuPipe,"'%s' using 1:2 w p pt 10 title 'Facility'",facility_output);
  fprintf(gnuPipe,", '%s' using 1:2 w p pt 7 lc rgb 'blue' title 'Demand'",demand_output);
  fprintf(gnuPipe,", '%s' using 1:2:($3*1.5) w p lt 2 pt 11 ps variable lc rgb 'dark-grey' title 'Opened'",centers_output);
  fprintf(gnuPipe,"\n");

  /*
  fprintf(gnuPipe,"'%s' using 1:2 with points title 'Demand'",demand_output);
  fprintf(gnuPipe,", '%s' using 1:2 with points title 'Facility'",facility_output);
  fprintf(gnuPipe,", '%s' using 1:2:($3*0.1) with circles title 'Opened'",centers_output);
  fprintf(gnuPipe,"\n");
  */
  pclose(gnuPipe);
  remove(demand_output);
  remove(facility_output);
  remove(centers_output);
}

void plot_solution_allocation(SQM_solution *X,double **f,string output,string suffix) {
  SQM_instance *I = X->get_instance();
  int n = I->potential_sites(),m = I->demand_points();
  int p = X->get_servers ();
  int j,r;
  int edge_key;
  fstream centersfile,edges_file;
  char demand_output[32],facility_output[32],centers_output[32],edges_output[32];
  
  sprintf(demand_output,"Tmp_demand_%d.dat",rand());
  gnuplot_write_points_file(demand_output,I->demand(0),m);

  sprintf(facility_output,"Tmp_facility_%d.dat",rand());
  gnuplot_write_points_file(facility_output,I->site(0),n);

  sprintf(centers_output,"Tmp_centers_%d.dat",rand());
  centersfile.open(centers_output,fstream::out);
  for (int j = 0;j < n;j++) {
    r = 0;
    for (int i = 0;i < p;i++) 
      if (X->get_server_location(i) == j)
	r++;
    if (r > 0)
      centersfile << I->site(j)->x << " "
		  << I->site(j)->y << " "
		  << r << endl;
  }
  centersfile.close();

  edge_key = rand();
  for (int i = 0;i < p;i++) {
    sprintf(edges_output,"Tmp_edges_center_%d_%d.dat",i+1,edge_key);
    edges_file.open(edges_output,fstream::out);
    j = X->get_server_location(i);
    for (int k = 0;k < m;k++) {
      if (100 * f[i][j] > 1)
	edges_file << I->site(j)->x << " " << I->site(j)->y << " "
		   << I->demand(k)->x << " " << I->demand(k)->y << " "
		   << ceil(sqrt(f[i][k]) * 100) << endl;
    }
    edges_file.close();
  }

  FILE *gnuPipe = popen("gnuplot","w");
  gnuplot_sets(gnuPipe);
  gnuplot_unsets(gnuPipe);
  for (int i = 1;i < 10;i++)
    fprintf(gnuPipe,"set for [i=%d:%d] style arrow i nohead lt %d lc 3 lw sqrt(i-%d)/3\n",11*(i-1)+1,11*i,10-i,11*(i-1));

  fprintf(gnuPipe,"set output '%s_%s.jpeg'\n",output.c_str(),suffix.c_str());
  fprintf(gnuPipe,"plot ");
  fprintf(gnuPipe,"'%s' using 1:2 w p lt 2 pt 10 title 'Facility'",facility_output);
  for (int i = 0;i < p;i++) {
    sprintf(edges_output,"Tmp_edges_center_%d_%d.dat",i+1,edge_key);
    fprintf(gnuPipe,", '%s' using 1:2:($3-$1):($4-$2):5 with vectors arrowstyle variable",edges_output);
  }
  fprintf(gnuPipe,", '%s' using 1:2 w p lt 2 pt 7 lc rgb 'blue' title 'Demand'",demand_output);
  fprintf(gnuPipe,", '%s' using 1:2:($3*1.5) w p lt 2 pt 11 ps variable lc rgb 'dark-grey' title 'Opened'",centers_output);
  fprintf(gnuPipe,"\n");

  /*
  for (int i = 0;i < p;i++) {
    sprintf(edges_output,"Tmp_edges_center_%d_%d.dat",i+1,edge_key);
    fprintf(gnuPipe,"set output '%s_center_%02d_%s.jpeg'\n",output.c_str(),i+1,suffix.c_str()); 
    fprintf(gnuPipe,"plot ");
    fprintf(gnuPipe,"'%s' using 1:2:($3-$1):($4-$2):5 with vectors arrowstyle variable",edges_output);
    fprintf(gnuPipe,", '%s' using 1:2:($3*1.5) w p lt 2 pt 11 ps variable lc rgb 'dark-grey' title 'Opened'",centers_output);
    fprintf(gnuPipe,", '%s' using 1:2 w p lt 2 pt 7 lc rgb 'blue' title 'Demand'",demand_output);
    fprintf(gnuPipe,"\n");
  }
  */

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
  fprintf(gnuPipe,"set term jpeg\n");
  /*fprintf(gnuPipe,"set key outside\n");*/
  fprintf(gnuPipe,"set grid ytics lt 0 lw 1 lc rgb '#bbbbbb'\n");
  fprintf(gnuPipe,"set grid xtics lt 0 lw 1 lc rgb '#bbbbbb'\n");
  fprintf(gnuPipe,"set grid\n");
}

void gnuplot_unsets(FILE *gnuPipe) {
  fprintf(gnuPipe,"unset key\n");
  fprintf(gnuPipe,"unset border\n");
  fprintf(gnuPipe,"unset yzeroaxis\n");
  /*fprintf(gnuPipe,"unset xtics\n");
  fprintf(gnuPipe,"unset ytics\n");
  fprintf(gnuPipe,"unset ztics\n");*/
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

void gnuplot_GRASP(char *filename) {
  FILE *gnuPipe = popen("gnuplot","w");
  gnuplot_sets(gnuPipe);
  fprintf(gnuPipe,"set output '%s.jpeg'\n",filename);
  gnuplot_unsets(gnuPipe); 
  fprintf(gnuPipe,"set key outside\n");
  fprintf(gnuPipe,"set xlabel 'alpha'\n");
  fprintf(gnuPipe,"set ylabel 'Objective'\n");
  
  fprintf(gnuPipe,"plot ");
  fprintf(gnuPipe,"'GRASP.dat' using 1:2 w lp title 'Best'");
  fprintf(gnuPipe,", 'GRASP.dat' using 1:3 w lp title 'Average'");
  fprintf(gnuPipe,", 'GRASP.dat' using 1:4 w lp title 'Worst'");
  pclose(gnuPipe);
}
