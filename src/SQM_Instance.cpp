 
#include <fstream>
#include <cstdlib>
#include "SQM_Instance.h"
#include "random.h"

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
  set_sites_distances();
}

SQM_instance::SQM_instance (string nodes) {
  double S;
  fstream demandfile;

  demandfile.open(nodes.c_str(),fstream::in);
  demandfile >> M >> S;

  /* Read demand points */
  V = new point[M];
  for (int i = 0;i < M;i++)
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
  set_sites_distances();
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
  set_sites_distances();
}

SQM_instance::~SQM_instance () {
  del_sites_distances();
  for (int j = 0;j < M;j++) delete [] Dist[j];
  delete [] Dist;
  delete [] V;
  delete [] W;
}

void SQM_instance::write (string demand_output,string facility_output) {
  fstream demandfile,facilityfile;
  
  demandfile.open(demand_output.c_str(),fstream::out);
  demandfile << M << endl;
  for (int i = 0;i < M;i++) 
    demandfile << V[i].x << " " 
	       << V[i].y << " "
	       << V[i].demand 
	       << endl;
  demandfile.close();
       
  facilityfile.open(facility_output.c_str(),fstream::out);
  facilityfile << N << endl;
  for (int i = 0;i < N;i++)
    facilityfile << W[i].x << " " 
		 << W[i].y
		 << endl;
  facilityfile.close();
}

void SQM_instance::set_sites_distances () {
  sites_dist = new  double*[N];
  for (int j = 0;j < N;j++)
    sites_dist[j] = new double [N];

  for (int i = 0;i < N;i++) {
    sites_dist[i][i] = 0;
    for (int j = i+1;j < N;j++) {
      sites_dist[i][j] = dist(&(W[j]),&(W[i]));
      sites_dist[j][i] = sites_dist[i][j];
    }
  }

}

void SQM_instance::del_sites_distances () {
  for (int j = 0;j < N;j++)
    delete [] sites_dist[j];
  delete [] sites_dist;
}

void SQM_instance::set_distances () {
  Dist = new  double*[M];
  for (int j = 0;j < M;j++)
    Dist[j] = new double [N];

  for (int j = 0;j < M;j++) {
    for (int i = 0;i < N;i++) {
      Dist[j][i] = dist(&(V[j]),&(W[i]));
    }
  }
}

int SQM_instance::demand_points () {
  return M;
}

int SQM_instance::potential_sites () {
  return N;
}

point* SQM_instance::site (int i) {
  return &(W[i]);
}

point* SQM_instance::demand (int j) {
  return &(V[j]);
}

double SQM_instance::distance (int i,int j) {
  return Dist[j][i];
}

double SQM_instance::sites_distance(int i,int j) {
  return sites_dist[i][j];
}

double SQM_instance::total_demand () {
  double demand;
  demand = 0.0;
  for (int k = 0;k < M;k++) demand += V[k].demand;
  return demand;
}

double SQM_instance::get_demand (int j) {
  if (j >= 0 && j < M)
    return V[j].demand;
  return -1;
}

double** SQM_instance::get_distances_matrix() {
  return Dist;
}

SQM_instance* SQM_load_instance(string filename,int M_clients,int N_sites) {
  SQM_instance *I;
  string Path = "../git/PMCLAP/Instancias/Q_MCLP_";
  /*I = read_points(demad_file.c_str());*/
  if ((M_clients == N_sites) && (filename == "Q_MCLP")) {
    switch (M_clients) {
    case 30 : 
      filename = Path+"30.txt";
      break;
    case 324 : 
      filename = Path+"324.txt";
      break;
    case 818 : 
      filename = Path+"818.txt";
      break;
    case 3283 : 
      filename = Path+"3283.txt";
      break;
    defautl:
      return NULL;
    }
    /*cout << "Read file: " << filename << endl;*/
    I = new SQM_instance(filename);
    return I;
  }
  if (file_exists(filename+"_demand.ins") &&
      file_exists(filename+"_facility.ins")) {
    I = new SQM_instance(filename+"_demand.ins",filename+"_facility.ins");
  }
  else {
    I = new SQM_instance(M_clients,N_sites);
    I->write(filename+"_demand.ins",filename+"_facility.ins");
  }
  return I;
}

bool file_exists (const string& name) {
  if (FILE *file = fopen(name.c_str(), "r")) {
    fclose(file);
    return true;
  }
  return false;
}
