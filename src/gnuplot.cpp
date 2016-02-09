
#include "gnuplot.h"

using namespace std;

void gnuplot_sets(FILE*);
void gnuplot_unsets(FILE*);
void gnuplot_set_arrow_styles(FILE*);
void gnuplot_write_points_file(char*,point*,int);

void plot_instance_solution(SQM_instance &Inst,int *Sol,string output) {
  int n = Inst.demand_points(),m = Inst.potential_sites();
  char demand_output[32],facility_output[32],centers_output[32];
  fstream centersfile;
  
  sprintf(demand_output,"Tmp_demand_%d.dat",rand());
  gnuplot_write_points_file(demand_output,Inst.demand(0),m);

  sprintf(facility_output,"Tmp_facility_%d.dat",rand());
  gnuplot_write_points_file(facility_output,Inst.site(0),m);

  sprintf(centers_output,"Tmp_centers_%d.dat",rand());
  centersfile.open(centers_output,fstream::out);
  for (int j = 0;j < n;j++) {
    if (Sol[j] > 0) {
      centersfile << Inst.site(j)->x << " "
		  << Inst.site(j)->y << " "
		  << Sol[j] << endl;
    }
  }

  FILE *gnuPipe = popen("gnuplot","w");
  gnuplot_sets(gnuPipe);
  fprintf(gnuPipe,"set output '%s.jpeg'\n",output.c_str());
  gnuplot_unsets(gnuPipe);

  fprintf(gnuPipe,"plot ");
  fprintf(gnuPipe,"'%s' using 1:2 w p pt 10 title 'Facility'",facility_output);
  fprintf(gnuPipe,", '%s' using 1:2 w p pt 7 lc rgb 'blue' title 'Demand'",
	  demand_output);
  fprintf(gnuPipe,", '%s' using 1:2:($3*1.5) w p lt 2 pt 11 ",centers_output);
  fprintf(gnuPipe,"ps variable lc rgb 'dark-grey' title 'Opened'");
  fprintf(gnuPipe,"\n");

  /*
  fprintf(gnuPipe,"'%s' using 1:2 with points title 'Demand'",demand_output);
  fprintf(gnuPipe,", '%s' using 1:2 with points title 'Facility'",
	  facility_output);
  fprintf(gnuPipe,", '%s' using 1:2:($3*0.1) with circles title 'Opened'",
	  centers_output);
  fprintf(gnuPipe,"\n");
  */

  pclose(gnuPipe);
  remove(demand_output);
  remove(facility_output);
  remove(centers_output);
}

void plot_solution_allocation(SQM_solution *X,double **f,string output,
			      string suffix)
{
  SQM_instance *I = X->get_instance();
  int n = I->potential_sites(),m = I->demand_points();
  int p = X->get_servers ();
  int j,r;
  int edge_key;
  fstream centersfile,edges_file;
  char demand_output[32],facility_output[32],centers_output[32],
    edges_output[32];
  
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
	edges_file 
	  << I->site(j)->x << " " 
	  << I->site(j)->y << " "
	  << I->demand(k)->x << " " 
	  << I->demand(k)->y << " "
	  << ceil(f[i][k] * 100) << endl;
    }
    edges_file.close();
  }

  FILE *gnuPipe = popen("gnuplot","w");
  gnuplot_sets(gnuPipe);
  gnuplot_unsets(gnuPipe);
  gnuplot_set_arrow_styles(gnuPipe);
  fprintf(gnuPipe,"set output '%s_%s.jpeg'\n",output.c_str(),suffix.c_str());
  fprintf(gnuPipe,"plot ");
  fprintf(gnuPipe,"'%s' using 1:2 w p lt 2 pt 10 title 'Facility'",
	  facility_output);
  for (int i = 0;i < p;i++) {
    sprintf(edges_output,"Tmp_edges_center_%d_%d.dat",i+1,edge_key);
    fprintf(gnuPipe,", '%s' using 1:2:($3-$1):($4-$2):5 ",edges_output);
    fprintf(gnuPipe,"with vectors arrowstyle variable");
  }
  fprintf(gnuPipe,", '%s' using 1:2 w p lt 2 pt 7 lc rgb 'blue' title 'Demand'",
	  demand_output);
  fprintf(gnuPipe,", '%s' using 1:2:($3*1.5) w p lt 2 pt 11 ",centers_output);
  fprintf(gnuPipe,"ps variable lc rgb 'dark-grey' title 'Opened'");
  fprintf(gnuPipe,"\n");

  /*
  for (int i = 0;i < p;i++) {
    sprintf(edges_output,"Tmp_edges_center_%d_%d.dat",i+1,edge_key);
    fprintf(gnuPipe,"set output '%s_center_%02d_%s.jpeg'\n",output.c_str(),i+1,
	    suffix.c_str()); 
    fprintf(gnuPipe,"plot ");
    fprintf(gnuPipe,"'%s' using 1:2:($3-$1):($4-$2):5 ",edges_output);
    fprintf(gnuPipe,"with vectors arrowstyle variable");
    fprintf(gnuPipe,", '%s' using 1:2:($3*1.5) w p lt 2 pt 11 ",centers_output);
    fprintf(gnuPipe,"ps variable lc rgb 'dark-grey' title 'Opened'");
    fprintf(gnuPipe,", '%s' using 1:2 w p lt 2 pt 7 ",demand_output);
    fprintf(gnuPipe,"lc rgb 'blue' title 'Demand'");
    fprintf(gnuPipe,"\n");
  }
  */

  pclose(gnuPipe);

  /* Remove temp files */
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

void gnuplot_set_arrow_styles(FILE *gnuPipe) {
  for (int i = 1;i < 10;i++) {
    fprintf(gnuPipe,"set for [i=%d:%d] style arrow i nohead ",11*(i-1)+1,11*i);
    fprintf(gnuPipe,"lt %d lc 3 lw sqrt(i-%d)/3",10-i,11*(i-1));
    fprintf(gnuPipe,"\n");
  }
}

void gnuplot_write_points_file(char *output,point *Set,int k) {
  fstream outputfile;
  outputfile.open(output,fstream::out);
  for (int i = 0;i < k;i++) {
    outputfile << Set[i].x << " " 
	       << Set[i].y << " "
	       << i+1 << endl;
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
  fprintf(gnuPipe,", 'GRASP.dat' using 1:3 w lp title 'Average+'");
  fprintf(gnuPipe,", 'GRASP.dat' using 1:4 w lp title 'Average'");
  fprintf(gnuPipe,", 'GRASP.dat' using 1:5 w lp title 'Worst'");
  pclose(gnuPipe);
}

void gnuplot_solution(SQM_solution *X,int sub) {
  FILE *gnuPipe;
  fstream centersfile;
  char demand_output[32],facility_output[32],centers_output[32];
  int key;
  SQM_instance *I = X->get_instance();
  int m = I->demand_points();
  int n = I->potential_sites();
  int p = X->get_servers();
  
  key = rand();
  sprintf(demand_output,"Tmp_demand_%d.dat",key);
  gnuplot_write_points_file(demand_output,I->demand(0),m);

  sprintf(facility_output,"Tmp_facility_%d.dat",key);
  gnuplot_write_points_file(facility_output,I->site(0),n);

  sprintf(centers_output,"Tmp_centers_%d.dat",key);
  centersfile.open(centers_output,fstream::out);
  for (int j = 0;j < n;j++) {
    int r = 0;
    for (int i = 0;i < p;i++) 
      if (X->get_server_location(i) == j)
	r++;
    if (r > 0)
      centersfile << I->site(j)->x << " "
		  << I->site(j)->y << " "
		  << r << endl;
  }
  centersfile.close();
 
  gnuPipe = popen("gnuplot","w");
  gnuplot_sets(gnuPipe);
  fprintf(gnuPipe,"set output './plots/Local_Search_step_%2d.jpeg'\n",sub);
  gnuplot_unsets(gnuPipe);

  fprintf(gnuPipe,"plot ");
  fprintf(gnuPipe,"'%s' using 1:2 ",facility_output);
  fprintf(gnuPipe,"w p lt 2 pt 10 title 'Facility'");
  fprintf(gnuPipe,", '%s' using 1:2:3 ",facility_output);
  fprintf(gnuPipe,"w labels point offset character 0,character 1");
  fprintf(gnuPipe,", '%s' using 1:2 ",demand_output);
  fprintf(gnuPipe,"w p lt 2 pt 7 lc rgb 'blue' title 'Demand'");
  fprintf(gnuPipe,", '%s' using 1:2:($3*1.5) w p lt 2 pt 11 ",centers_output);
  fprintf(gnuPipe,"ps variable lc rgb 'dark-grey' title 'Opened'");
  fprintf(gnuPipe,"\n");

  pclose(gnuPipe);
  remove(demand_output);
  remove(facility_output);
  remove(centers_output);
}

void gnuplot_Path_Relinking
(SQM_solution &A,int *match,SQM_solution &B,string filename) {
  FILE *gnuPipe;
  fstream centersfile,edgesfile;
  char demand_output[32],facility_output[32],edges_output[32],
    centersA_output[32],centersB_output[32];
  int p = A.get_servers();
  SQM_instance *I = A.get_instance();
  int m = I->demand_points();
  int n = I->potential_sites();
  int key;
  
  key = rand();
  sprintf(demand_output,"Tmp_demand_%d.dat",key);
  gnuplot_write_points_file(demand_output,I->demand(0),m);

  sprintf(facility_output,"Tmp_facility_%d.dat",key);
  gnuplot_write_points_file(facility_output,I->site(0),n);

  sprintf(centersA_output,"Tmp_centersA_%d.dat",key);
  centersfile.open(centersA_output,fstream::out);
  for (int j = 0;j < n;j++) {
    int r = 0;
    for (int i = 0;i < p;i++) 
      if (A.get_server_location(i) == j)
	r++;
    if (r > 0)
      centersfile << I->site(j)->x << " "
		  << I->site(j)->y << " "
		  << r << endl;
  }
  centersfile.close();
  sprintf(centersB_output,"Tmp_centersB_%d.dat",key);
  centersfile.open(centersB_output,fstream::out);
  for (int j = 0;j < n;j++) {
    int r = 0;
    for (int i = 0;i < p;i++) 
      if (B.get_server_location(i) == j)
	r++;
    if (r > 0)
      centersfile
	<< I->site(j)->x << " " << I->site(j)->y << " "
	<< r << endl;
  }
  centersfile.close();
  
  sprintf(edges_output,"Tmp_edges_%d.dat",key);
  edgesfile.open(edges_output,fstream::out);
  for (int i = 0;i < p;i++) {
    int loc_a = A.get_server_location(i);
    int loc_b = B.get_server_location(match[i]);
    if (loc_a != loc_b)
      edgesfile
	<<  I->site(loc_a)->x << " " << I->site(loc_a)->y << " "
	<<  I->site(loc_b)->x << " " << I->site(loc_b)->y
	<< endl;
  }
  edgesfile.close();

  gnuPipe = popen("gnuplot","w");
  gnuplot_sets(gnuPipe);
  fprintf(gnuPipe,"set output './plots/PR_%s.jpeg'\n",filename.c_str());
  gnuplot_unsets(gnuPipe);

  fprintf(gnuPipe,"plot ");
  fprintf(gnuPipe,"'%s' using 1:2 ",facility_output);
  fprintf(gnuPipe,"w p lt 2 pt 10 title 'Facility'");
  /*
  fprintf(gnuPipe,", '%s' using 1:2:3 ",facility_output);
  fprintf(gnuPipe,"w labels point offset character 0,character 1");
  */
  fprintf(gnuPipe,", '%s' using 1:2 ",demand_output);
  fprintf(gnuPipe,"w p lt 2 pt 7 lc rgb 'blue' title 'Demand'");
  fprintf(gnuPipe,", '%s' using 1:2:($3-$1):($4-$2) ",edges_output);
  fprintf(gnuPipe,"with vectors");
  fprintf(gnuPipe,", '%s' using 1:2:($3*1.5) w p lt 2 pt 11 ",centersA_output);
  fprintf(gnuPipe,"ps variable lc rgb 'dark-grey' title 'Servers A'");
  fprintf(gnuPipe,", '%s' using 1:2:($3*1.5) w p lt 2 pt 11 ",centersB_output);
  fprintf(gnuPipe,"ps variable lc rgb 'red' title 'Servers B'");
  fprintf(gnuPipe,"\n");

  pclose(gnuPipe);
  remove(demand_output);
  remove(facility_output);
  remove(centersA_output);
  remove(centersB_output);
  remove(edges_output);
}
