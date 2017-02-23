/*!
 * \brief Modet that determines the location of adjusters
 * 
 * Minimizing the maximum workload per agent
 */

#include "Goldberg.h"
#include "log.h"

void gnuplot_goldberg
(SQM_instance&,
 int,
 IloCplex&,
 IloBoolVarArray&,
 BoolVarArrayMatrix&
 );

void Goldberg
(SQM_instance &Inst, //!< Set of points
 int p, //!< facilities
 float mu, //!< rate parameter
 float f //!<
 )
{
  IloEnv env;
  try {
    clock_t clocks = clock();
    int i,j,k,r;
    IloInt n = Inst.potential_sites();
    IloInt m = Inst.demand_points();
    IloNum f_i;
    IloNum rho;
    IloNum M = 10000.0;
          
    // rho from Daskin
    rho = 0.0;
    for (i = 0;i < m;i++)
      rho += f * Inst.get_demand(i);
    rho /= (mu * p);
    rho /= Inst.total_demand();
    if (rho > 1) {
      logError(cout << "rho must be between (0,1)! " << endl
	       << "f = " << f << endl
	       << "mu = " << mu << endl
	       << "p = " << p << endl);
    }
    logDebug(cout << "rho = " << rho << endl);
    // rho from ReVelle & Hogan
    // pendiente

    logInfo(cout << "Comienza definicion del Modelo Goldberg" << endl);
    IloModel modelo(env);
    
    logDebug(cout << "++ Variables ++" << endl);
    IloNumVar S(env,0,IloInfinity,IloNumVar::Float,"S");
    IloBoolVarArray x(env);
    BoolVarArrayMatrix y(env,m);
    
    char VarName[16];
    logDebug(cout << "Variables x_j if an adjuster is located at potential site j" << endl);
    for(i = 0;i < n;i++){
      sprintf(VarName,"x%d",i+1);
      x.add(IloBoolVar(env,VarName));
    }

    logDebug(cout << "Variables y_ijk if the adjuster at j is the k-th in cover i" << endl);
    for (i = 0;i < m;i++) {
      BoolVarMatrix y_i(env,n);
      for (j = 0;j < n;j++) {
	IloBoolVarArray y_i_j(env);
	for (k = 0;k < p;k++) {
	  sprintf(VarName,"y_%d_%d_%d",i+1,j+1,k+1);
	  y_i_j.add(IloBoolVar(env,VarName));
	}
	y_i[j] = y_i_j;
      }
      y[i] = y_i;
    }

    logDebug(cout << "++ Constraints ++" << endl);
    logDebug(cout << "Solo una instalacion ocupa la posion k del cliente i" << endl);
    for (i = 0;i < m;i++) {
      for (k = 0;k < p;k++) {
	IloExpr Cover(env);
	for (j = 0;j < n;j++) Cover += y[i][j][k];	
	modelo.add(Cover == 1);
	Cover.end();
      }
    }

    logDebug(cout << "La instalacion j solo puede ocupar una posicion del cliente i" << endl);
    for (i = 0;i < m;i++) {
      for (j = 0;j < n;j++) {
	modelo.add(IloSum(y[i][j]) <= x[j]);
      }
    }

    logDebug(cout << "Instalaciones a abrir" << endl);
    modelo.add(IloSum(x) == p);

    logDebug(cout << "Relacionar variables de localizacion y asignacion" << endl);
    for(i = 0;i < m;i++){
      for (j = 0;j < n;j++) {
	for (k = 0;k < p;k++) 
	  modelo.add(y[i][j][k] <= x[j]);
      }
    }

    logDebug(cout << "Restricion de orden de asignaicion" << endl);
    NumMatrix Dist(env,m);
    for (i = 0;i < m;i++)
      Dist[i] = IloNumArray(env,n);
    for (i = 0;i < m;i++) {
      for (j = 0;j < n;j++) {
	Dist[i][j] = Inst.distance(i,j);
      }
    }
      
    for (i = 0;i < m;i++) {
      for (k = 1;k < p;k++) {
	for (j = 0;j < n;j++) {
	  IloExpr balance(env);
	  for (r = 0;r < n;r++) {
	    if (Dist[i][r] <= Dist[i][j])
	      balance += y[i][r][k-1];
	  }
	  modelo.add(y[i][j][k] <= balance);
	  balance.end();
	}
      }
    }

    // Cargas de trabajo
    IloNum coef; 
    IloExpr workload(env);
    for (k = 0;k < p;k++) {
      coef = (1 - rho) * pow(rho,k);
      for (i = 0;i < m;i++) {
	f_i = f * Inst.get_demand(i);
	for (j = 0;j < n;j++) {
	  workload += f_i * coef * Dist[i][j]* y[i][j][k];
	}
      }
    }
    
    logDebug(cout << "++ Funcion Objetivo ++" << endl);
    logDebug(cout << "Minimize the maximum workload" << endl);
    modelo.add(IloMinimize(env,workload));
    clocks = clock() - clocks;
    results << "," << clocks / CLOCKS_PER_SEC;

    cout << "Resuelve modelo" << endl;
    IloCplex cplex(modelo);
    stringstream ModelName;
    ModelName << "Goldberg-model_" << m << "_" << n << "_" << p << ".lp";
    cplex.exportModel(ModelName.str().c_str());
    results << "," << ModelName.str();

    LogFile << "Solve Goldberg model" << endl;
    LogFile << endl << "** Cplex Start **" << endl;
    cplex.setOut(LogFile);
    cplex.setParam(IloCplex::TiLim,TIME_MAX);

    double CplexTime = cplex.getCplexTime();
    if (cplex.solve()) {

      double gap = (cplex.getObjValue() - cplex.getBestObjValue()) / cplex.getBestObjValue();
      results << "," << cplex.getCplexTime() - CplexTime
	      << "," << cplex.getStatus()
	      << "," << cplex.getObjValue()
	      << "," << gap
	      << endl;
      
      LogFile << "Solution status: " << cplex.getStatus() << endl;
      LogFile << "Maximum profit = " << cplex.getObjValue() << endl;
      for (j = 0;j < n;j++) {
	if (cplex.getValue(x[j]) > 0.5) LogFile << j+1 << " ";
	//cout << cplex.getValue(y[j]);
      }
      cout << endl;

      /*for(i = 0;i < n;i++){
	cout << "-";
      }
      cout << endl;
      for(i = 0;i < n;i++){
	for(j = 0;j < n;j++){
	  cout << cplex.getValue(x[i][j]);
	}
	cout << endl;
      }/**/

      gnuplot_goldberg(Inst,p,cplex,x,y);
    }
    else {
      cout << "No solution found" << endl;
    }
    
  }
  catch (IloException& e) {
    cerr << "Concert exception caught: " << e << endl;
  }
  catch (...) {
    cerr << "Unknown exception caught" << endl;
  }
  env.end();
  
}

void gnuplot_goldberg
(SQM_instance &Inst,
 int p,
 IloCplex &cplex,
 IloBoolVarArray &location,
 BoolVarArrayMatrix &assignment
 )
{
  int i,j,k;
  int n = Inst.potential_sites(),m = Inst.demand_points();
  point *client = Inst.demand(0);
  point *potential_site = Inst.site(0);
  char outfilename[32],centersfilename[32],clientsfilename[32];

  sprintf(centersfilename,"Tmp_centers_%d.dat",n);
  ofstream centros(centersfilename);

  for(j = 0;j < n;j++){
    if (cplex.getValue(location[j]))
    centros << potential_site[j].x << " " 
	    << potential_site[j].y << " "
	    << cplex.getValue(location[j]) << endl;
  }
  centros.close();

  sprintf(clientsfilename,"Tmp_clients_%d.dat",n);
  ofstream clients(clientsfilename);

  for(i = 0;i < m;i++){
    clients << client[i].x << " " << client[i].y << endl;
  }
  clients.close();

  for (k = 0;k < p;k++) {
    sprintf(outfilename,"Tmp_edges_%d_%d_%d.dat",m,n,k+1);
    ofstream outfile(outfilename);

    for(i = 0;i < m;i++){
      for(j = 0;j < n;j++){
	if (cplex.getValue(assignment[i][j][k]) > 0.5)
	  outfile << client[i].x << " " 
		  << client[i].y << " " 
		  << potential_site[j].x - client[i].x << " " 
		  << potential_site[j].y - client[i].y << " "
		  << j+1 << endl;
      }
    }

    outfile.close();
  }

  FILE *gnuPipe = popen("gnuplot","w");
  fprintf(gnuPipe,"set term svg\n");
  /*fprintf(gnuPipe,"set key outside\n");*/
  fprintf(gnuPipe,"unset key\n");
  fprintf(gnuPipe,"unset border\n");
  fprintf(gnuPipe,"unset yzeroaxis\n");
  fprintf(gnuPipe,"unset xtics\n");
  fprintf(gnuPipe,"unset ytics\n");
  fprintf(gnuPipe,"unset ztics\n");

  /* plot full */
  fprintf(gnuPipe,"set output 'Goldberg_%d_%d_%d.svg'\n",m,n,p);
  fprintf(gnuPipe,"plot ");
  fprintf(gnuPipe,"'%s' using 1:2 with points lc rgb \"black\" title 'Demand'",clientsfilename);
  fprintf(gnuPipe,", '%s' using 1:2 with points lc rgb \"red\" title 'Facility'",centersfilename);
  fprintf(gnuPipe,", '%s' using 1:($3 > 0.5 ? $2 : 1/0):(10) with circles lc rgb 'blue' title 'Opened'",centersfilename);
  for (k = 0;k < p;k++) {
    sprintf(outfilename,"Tmp_edges_%d_%d_%d.dat",m,n,k+1);
    fprintf(gnuPipe,", '%s' using 1:2:3:4 with vectors nohead",outfilename);
  }
  fprintf(gnuPipe,"\n");

  for (k = 0;k < p;k++) {
    sprintf(outfilename,"Tmp_edges_%d_%d_%d.dat",m,n,k+1);
    fprintf(gnuPipe,"set output 'Goldberg_%d_%d_%d_order_%02d.svg'\n",m,n,p,k+1);
    fprintf(gnuPipe,"plot ");
    fprintf(gnuPipe,"'%s' using 1:2 with points lc rgb \"black\"",clientsfilename);
    fprintf(gnuPipe,", '%s' using 1:2 with points lc rgb \"red\"",centersfilename);
    fprintf(gnuPipe,", '%s' using 1:2:3:4 with vectors nohead",outfilename);
    fprintf(gnuPipe,"\n");
  }

  for (j = 0;j < n;j++) {
    if (cplex.getValue(location[j]) > 0.5) {
      for (k = 0;k < p;k++) {
	sprintf(outfilename,"Tmp_edges_%d_%d_%d.dat",m,n,k+1);
	fprintf(gnuPipe,"set output 'Goldberg_%d_%d_%d_center_%02d_order_%02d.svg'\n",m,n,p,j+1,k+1);
	fprintf(gnuPipe,"plot ");
	fprintf(gnuPipe,"'%s' using 1:2 with points lc rgb \"black\" title 'Demand'",clientsfilename);
	fprintf(gnuPipe,", '%s' using 1:2 with points lc rgb \"red\" title 'Facility'",centersfilename);
	fprintf(gnuPipe,", '%s' using 1:($5 == %d ? $2 : 1/0):3:4 with vectors nohead title 'Order %d'",outfilename,j+1,k+1);
	fprintf(gnuPipe,"\n");
      }
    }
  }

  pclose(gnuPipe);
  remove(clientsfilename);
  remove(centersfilename);
  for (k = 0;k < p;k++) {
    sprintf(outfilename,"Tmp_edges_%d_%d_%d.dat",m,n,k+1);
    remove(outfilename);
  }
}
