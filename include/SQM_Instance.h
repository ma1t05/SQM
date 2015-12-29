
#ifndef _SQM_INSTANCE_H
#define _SQM_INSTANCE_H 1

#include "point.h"
#include <string>

using namespace std;

#define MIN_X 0.0
#define MAX_X 1024.0
#define MIN_Y 0.0
#define MAX_Y 1024.0

class SQM_instance {
private:
  int M; /* Number of demand points */
  int N; /* Number of potencial sites to locate a server */
  point *V; /* Set of demand points */
  point *W; /* Set of potencial locations sites */
  int **a; /* List of lists of preferred sites */
  double **Dist; /* Matrix of distances beetwen demand and sites */
  double **sites_dist; /* Matrix of distances beetwen sites */
  void set_distances();
  void set_sites_distances ();
  void del_sites_distances ();
  double lambda; /* total arrival rate */
  double speed; /* travel speed */
  double mu; /* Rate of service per server */
public:
  SQM_instance(int m,int n);
  SQM_instance(string);
  SQM_instance(string,string);
  ~SQM_instance();
  void write(string,string);
  int demand_points();
  int potential_sites();
  point* site(int /* site*/);
  point* demand(int /* demand */);
  /*int site_order(int i,int j);*/
  double distance(int /* site */,int /* demand */);
  double sites_distance(int /* site */,int /* site */);
  double total_demand ();
  double get_demand(int);
  double** get_distances_matrix();
};

SQM_instance* SQM_load_instance(string,int,int);
bool file_exists (const string&);

#endif

/* eof */
