/*
 * Modelo de Aproxmiacion al problema de Berman
 * 
 */

#include "SQM_model.h"
#include "log.h"

void SQM_model
(SQM_instance &Inst, // Set of points
 int p,              // facilities
 float speed         // speed
) {
  logInfo(cout << "Start: SQM_model" << endl);
  float mu = Inst.get_service_rate (); // rate parameter
  float f = Inst.get_arrival_rate (); // portion of demand
  IloEnv env;
  try {
    int i,j,l,k,r;
    IloInt m,n = Inst.potential_sites();
    IloNum f_i;
    IloNum rho;
    IloNum M = 10000.0;
    IloNum beta = 1.5;
    k = 3;
    m = n;
    point *puntos = Inst.demand(0);
    logDebug(cout << "Comienza definicion del Modelo B" << endl);
    IloModel modelo(env);
    
    logDebug(cout << "++ Variables ++" << endl);
    IloIntVarArray x(env);
    BoolVarArrayMatrix y(env,n);
    BoolVarMatrix u(env,n);
    BoolVarMatrix v(env,n);
    
    logDebug(cout << "Nombre las variables para facil identificacion" << endl);
    char VarName[16];
    for(i = 0;i < n;i++){
      sprintf(VarName,"x%d",i+1);
      x.add(IloIntVar(env,VarName));
      logDebug(cout << VarName << "\t");
    }
    logDebug(cout << endl);

    for (i = 0;i < n;i++) {
      BoolVarMatrix y_i(env,n);
      for (j = 0;j < n;j++) {
	IloBoolVarArray y_i_j(env);
	for (l = 0;l < k;l++) {
	  sprintf(VarName,"y_%d_%d_%d",i+1,j+1,l+1);
	  y_i_j.add(IloBoolVar(env,VarName));
	  logDebug(cout << VarName << "\t");
	}
	y_i[j] = y_i_j;
	logDebug(cout << endl);
      }
      y[i] = y_i;
      logDebug(cout << endl);
    }
    logDebug(cout << endl);

    for (i = 0;i < n;i++) {
      IloBoolVarArray u_i(env);
      for (j = 0;j < n;j++) {
	sprintf(VarName,"u_%d_%d",i+1,j+1);
	u_i.add(IloBoolVar(env,VarName));
	logDebug(cout << VarName << "\t");
      }
      u[i] = u_i;
      logDebug(cout << endl);
    }
    logDebug(cout << endl);

    for (i = 0;i < n;i++) {
      IloBoolVarArray v_i(env);
      for (j = 0;j < n;j++) {
	sprintf(VarName,"v_%d_%d",i+1,j+1);
	v_i.add(IloBoolVar(env,VarName));
	logDebug(cout << VarName << "\t");
      }
      v[i] = v_i;
      logDebug(cout << endl);
    }
    logDebug(cout << endl);

    logDebug(cout << "++ Restricciones ++" << endl);
    logDebug(cout << "Instalaciones a abrir" << endl);
    modelo.add(IloSum(x) == p);

    logDebug(cout << "Relacionar variables de localizacion y asignacion" << endl);
    for(i = 0;i < n;i++){
      for (j = 0;j < n;j++) {
	for (l = 0;l < k;l++) 
	  modelo.add(y[i][j][l] <= x[j]);
      }
    }

    logDebug(cout << "Solo una instalacion ocupa la posion k del cliente i" << endl);
    for (i = 0;i < n;i++) {
      for (l = 0;l < k;l++) {
	IloExpr Cover(env);
	for (j = 0;j < n;j++) Cover += y[i][j][l];
	modelo.add(Cover == 1);
	Cover.end();
      }
    }

    logDebug(cout << "Los ajustadores en j solo puede ocupar x_j posiciones del cliente i" << endl);
    for (i = 0;i < n;i++) {
      for (j = 0;j < n;j++) {
	modelo.add(IloSum(y[i][j]) <= x[j]);
      }
    }

    logDebug(cout << "Restricion de orden de asignaicion" << endl);
    NumMatrix O(env,n);
    for (i = 0;i < n;i++)
      O[i] = IloNumArray(env,n);
    for (i = 0;i < n;i++) {
      O[i][i] = 0;
      for (j = i+1;j < n;j++) {
	O[i][j] = Inst.distance(i,j);
	O[j][i] = O[i][j];
      }
    }
      
    /*  */
    for (i = 0;i < n;i++) {
      for (l = 1;l < k;l++) {
	for (j = 0;j < m;j++) {
	  IloExpr balance(env);
	  for (r = 0;r < m;r++) {
	    if (O[i][r] <= O[i][j])
	      balance += y[i][r][l-1];
	  }
	  modelo.add(y[i][j][l] <= balance);
	  balance.end();
	}
      }
    } /* */

    logDebug(cout << "Restricciones disjuntas" << endl);
    /* Restriccion para realacionar 
       u_{ij} = 
                1: Si la suma de instalaciones mas cercanas a i hasta incluir j es menor o igual a k
		0: Si no
    */
    for (i = 0;i < n;i++) {
      for (j = 0;j < m;j++) {
	IloExpr Instalaciones(env);
	for (r = 0;r < m;r++) {
	  if (O[i][r] <= O[i][j])
	    Instalaciones += x[r];
	}
	modelo.add(Instalaciones <= p - (p - k) * u[i][j]);
	modelo.add(Instalaciones + M * u[i][j] >= k + 1);
	Instalaciones.end();
      }
    }

    /* Si u_{ij} = 1, todas las instalaciones de j deben de ser asignadas a i */
    for (i = 0;i < n;i++) {
      for (j = 0;j < m;j++) {
	modelo.add(IloSum(y[i][j]) + M * (1 - u[i][j]) >= x[j]);
      }
    }

    /* Restriccion para relacionar
       v_{ij} =
                1: Si la suma de instalaciones mas cercanas a i hasta antes de j es menor a k
		0: Si no
     */
    for (i = 0;i < n;i++) {
      for (j = 0;j < m;j++) {
	if ( i != j) {
	  IloExpr Instalaciones(env);
	  for (r = 0;r < m;r++) {
	    if (O[i][r] < O[i][j])
	      Instalaciones += x[r];
	  }
	  modelo.add(Instalaciones <= p - (p - (k - 1)) * v[i][j]);
	  modelo.add(Instalaciones + M * v[i][j] >= k);

	  /* Si v_{ij} = 1 y u_{ij} = 0 alguntas de las instalaciones de j seran asignadas a i para completar las k */
	  modelo.add(IloSum(y[i][j]) + M * (1 - v[i][j] + u[i][j]) >= k - Instalaciones);
	  modelo.add(IloSum(y[i][j]) - M * (1 - v[i][j] + u[i][j]) <= k - Instalaciones);

	  Instalaciones.end();
	}
      }
    }

    /* En caso de ser u_{ij} = 0 y v_{ij] = 0, no se asinga ninguna instalacion de j a i*/
    for (i = 0;i < n;i++) {
      for (j = 0;j < m;j++) {
	for (l = 1;l < k;l++) {
	  modelo.add(y[i][j][l] <= u[i][j] + v[i][j]);
	}
      }
    }

    /* Comienza definicion de funcion objetivo */
    // Cargas de trabajo
    IloNum coef;
    // rho from Daskin
    rho = 0.0;
    for (i = 0;i < n;i++)
      rho += f * Inst.get_demand(i);
    rho /= (mu * p);
    LogFile << "rho = " << rho << endl;
    // rho from ReVelle & Hogan
    // pendiente
    
    IloExpr workload(env);
    IloNum time_per_km = beta / speed;
    IloNum wt = 0.0;
    for (l = 0;l < k;l++) {
      coef = (1 - rho) * pow(rho,l);
      for (i = 0;i < n;i++) {
	f_i = f * Inst.get_demand(i);
	for (j = 0;j < m;j++) {
	  workload += f_i * coef * (O[i][j]*time_per_km + wt)* y[i][j][l];
	}
      }
    }
    
    logDebug(cout << "Funcion Objetivo" << endl);
    modelo.add(IloMaximize(env,workload));
    workload.end();

    logInfo(cout << "Solve model" << endl);
    IloCplex cplex(modelo);
    stringstream ModelName;
    ModelName << "SQM-model_" << n << "_" << m << "_" << p << ".lp";
    cplex.exportModel(ModelName.str().c_str());
    cplex.setParam(IloCplex::TiLim,TIME_MAX);
    if (cplex.solve()) {
      LogFile << "Solution status: " << cplex.getStatus() << endl;
      LogFile << "Maximum profit = " << cplex.getObjValue() << endl;
      for (j = 0;j < n;j++) {
	if (cplex.getValue(x[j]) > 0) {
	  LogFile << j+1 << ": "
		  << cplex.getValue(x[j]) << endl;
	} 
	//LogFile << cplex.getValue(y[j]);
      }
      LogFile << endl;
      /*for(i = 0;i < n;i++){
	LogFile << "-";
      }
      LogFile << endl;
      for(i = 0;i < n;i++){
	for(j = 0;j < n;j++){
	  LogFile << cplex.getValue(x[i][j]);
	}
	LogFile << endl;
      }/**/
    }
    else {
      LogFile << "No solution found" << endl;
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

int* SQM_model
(
 SQM_instance &Inst, // Set of points
 int p, //!< facilities
 int k, //!< Number of facilities that care
 float speed //!< speed
 )
{
  logInfo(cout << "Start: SQM_model" << endl);
  int *Sol;
  float mu = Inst.get_service_rate(); // rate parameter
  float f = Inst.get_service_rate(); // portion of demand
  IloEnv env; //!< Instanciates Ilo Enviroment
  try {
    clock_t clocks,start = clock();
    int i,j,l,r;
    IloInt m,n;
    IloNum f_i;
    IloNum rho;
    IloNum M = 10000.0;
    IloNum beta = 1.5;
    m = Inst.demand_points();
    n = Inst.potential_sites();
    logDebug(cout << "Comienza definicion del Modelo B" << endl);
    IloModel modelo(env);
    
    logDebug(cout << "++ Variables ++" << endl);
    IloNumVar S(env,0,IloInfinity,IloNumVar::Float,"S");
    IloIntVarArray x(env);
    BoolVarArrayMatrix y(env,m);
    BoolVarMatrix u(env,m);
    BoolVarMatrix v(env,m);
    
    logDebug(cout << "Nombra variables para facil identificacion" << endl);
    char VarName[16];
    for(i = 0;i < n;i++){
      sprintf(VarName,"x%d",i+1);
      x.add(IloIntVar(env,VarName));
    }

    for (i = 0;i < m;i++) {
      BoolVarMatrix y_i(env,n);
      for (j = 0;j < n;j++) {
	IloBoolVarArray y_i_j(env);
	for (l = 0;l < k;l++) {
	  sprintf(VarName,"y_%d_%d_%d",i+1,j+1,l+1);
	  y_i_j.add(IloBoolVar(env,VarName));
	}
	y_i[j] = y_i_j;
      }
      y[i] = y_i;
    }

    for (i = 0;i < m;i++) {
      IloBoolVarArray u_i(env);
      for (j = 0;j < n;j++) {
	sprintf(VarName,"u_%d_%d",i+1,j+1);
	u_i.add(IloBoolVar(env,VarName));
      }
      u[i] = u_i;
    }

    for (i = 0;i < m;i++) {
      IloBoolVarArray v_i(env);
      for (j = 0;j < n;j++) {
	sprintf(VarName,"v_%d_%d",i+1,j+1);
	v_i.add(IloBoolVar(env,VarName));
      }
      v[i] = v_i;
    }

    logDebug(cout << "++ Restricciones ++" << endl);
    logDebug(cout << "Instalaciones a abrir" << endl);
    modelo.add(IloSum(x) == p);

    logDebug(cout << "Relacionar variables de localizacion y asignacion" << endl);
    for(i = 0;i < m;i++){
      for (j = 0;j < n;j++) {
	for (l = 0;l < k;l++) 
	  modelo.add(y[i][j][l] <= x[j]);
      }
    }

    logDebug(cout << "Solo una instalacion ocupa la posion k del cliente i" << endl);
    for (i = 0;i < m;i++) {
      for (l = 0;l < k;l++) {
	IloExpr Cover(env);
	for (j = 0;j < n;j++) Cover += y[i][j][l];
	modelo.add(Cover == 1);
	Cover.end();
      }
    }

    logDebug(cout << "Los ajustadores en j solo puede ocupar x_j posiciones del cliente i" << endl);
    for (i = 0;i < m;i++) {
      for (j = 0;j < n;j++) {
	modelo.add(IloSum(y[i][j]) <= x[j]);
      }
    }

    logDebug(cout << "Restricion de orden de asignaicion" << endl);
    NumMatrix O(env,m);
    for (i = 0;i < m;i++)
      O[i] = IloNumArray(env,n);
    for (i = 0;i < m;i++)
      for (j = 0;j < n;j++)
	O[i][j] = Inst.distance(i,j);
      
    /*  */
    for (i = 0;i < m;i++) {
      for (l = 1;l < k;l++) {
	for (j = 0;j < n;j++) {
	  IloExpr balance(env);
	  for (r = 0;r < n;r++) {
	    if (O[i][r] <= O[i][j])
	      balance += y[i][r][l-1];
	  }
	  modelo.add(y[i][j][l] <= balance);
	  balance.end();
	}
      }
    } /* */

    logDebug(cout << "Restricciones disjuntas" << endl);
    /* Restriccion para realacionar 
       u_{ij} = 
                1: Si la suma de instalaciones mas cercanas a i hasta incluir j es menor o igual a k
		0: Si no
    */
    for (i = 0;i < m;i++) {
      for (j = 0;j < n;j++) {
	IloExpr Instalaciones(env);
	for (r = 0;r < n;r++) {
	  if (O[i][r] <= O[i][j])
	    Instalaciones += x[r];
	}
	modelo.add(Instalaciones <= p - (p - k) * u[i][j]);
	modelo.add(Instalaciones + M * u[i][j] >= k + 1);
	Instalaciones.end();
      }
    }

    /* Si u_{ij} = 1, todas las instalaciones de j deben de ser asignadas a i */
    for (i = 0;i < m;i++) {
      for (j = 0;j < n;j++) {
	modelo.add(IloSum(y[i][j]) + M * (1 - u[i][j]) >= x[j]);
      }
    }

    /* Restriccion para relacionar
       v_{ij} =
                1: Si la suma de instalaciones mas cercanas a i hasta antes de j es menor a k
		0: Si no
     */
    for (i = 0;i < m;i++) {
      for (j = 0;j < n;j++) {
	if ( i != j) {
	  IloExpr Instalaciones(env);
	  for (r = 0;r < n;r++) {
	    if (O[i][r] < O[i][j])
	      Instalaciones += x[r];
	  }
	  modelo.add(Instalaciones <= p - (p - (k - 1)) * v[i][j]);
	  modelo.add(Instalaciones + M * v[i][j] >= k);

	  /* Si v_{ij} = 1 y u_{ij} = 0 alguntas de las instalaciones de j seran asignadas a i para completar las k */
	  modelo.add(IloSum(y[i][j]) + M * (1 - v[i][j] + u[i][j]) >= k - Instalaciones);
	  modelo.add(IloSum(y[i][j]) - M * (1 - v[i][j] + u[i][j]) <= k - Instalaciones);

	  Instalaciones.end();
	}
      }
    }

    /* En caso de ser u_{ij} = 0 y v_{ij] = 0, no se asinga ninguna instalacion de j a i*/
    for (i = 0;i < m;i++) {
      for (j = 0;j < n;j++) {
	for (l = 1;l < k;l++) {
	  modelo.add(y[i][j][l] <= u[i][j] + v[i][j]);
	}
      }
    }

    /* Comienza definicion de funcion objetivo */
    // Cargas de trabajo
    IloNum coef;
    // rho from Daskin
    rho = 0.0;
    for (i = 0;i < m;i++)
      rho += f * Inst.get_demand(i);
    rho /= (mu * p);
    LogFile << "rho = " << rho << endl;
    // rho from ReVelle & Hogan
    // pendiente
    
    IloNum time_per_km = beta / speed;
    IloNum wt = 0.0;
    for (j = 0;j < n;j++) {
      IloExpr workload(env);
      for (l = 0;l < k;l++) {
	coef = (1 - rho) * pow(rho,l);
	for (i = 0;i < m;i++) {
	  f_i = f * Inst.get_demand(i);
	  workload += f_i * coef * (O[i][j]*time_per_km + wt)* y[i][j][l];
	}
      }
      modelo.add(S >= workload);
    }
    
    logDebug(cout << "++ Funcion Objetivo ++" << endl);
    modelo.add(IloMinimize(env,S));
    clocks = clock() - start;
    results << "," << clocks / CLOCKS_PER_SEC;

    IloCplex cplex(modelo);
    stringstream ModelName;
    ModelName << "SQM-model_" << m << "_" << n << "_" << p << "_" << k << ".lp";
    cplex.exportModel(ModelName.str().c_str());
    results << "," << ModelName.str();

    LogFile << "Solve SQM model" << endl;
    LogFile << endl << "** Cplex Start **" << endl;
    cplex.setOut(LogFile);
    cplex.setParam(IloCplex::TiLim,TIME_MAX);

    double CplexTime = cplex.getCplexTime();
    if(cplex.solve()){

      double gap = (cplex.getObjValue() - cplex.getBestObjValue()) / cplex.getBestObjValue();
      results << "," << cplex.getCplexTime() - CplexTime
	      << "," << cplex.getStatus()
	      << "," << cplex.getObjValue()
	      << "," << gap
	      << endl;

      Sol = new int[n];
      LogFile << "** Cplex Ends **" << endl << endl;
      LogFile << "Solution status: " << cplex.getStatus() << endl;
      LogFile << "Maximum profit = " << cplex.getObjValue() << endl;
      for(j = 0;j < n;j++){
	Sol[j] = 0;
	if(cplex.getValue(x[j]) > 0.5) {
	  Sol[j] = cplex.getValue(x[j]);
	  LogFile << j+1 << ": "
	       << Sol[j] << endl;
	    } 
	//LogFile << cplex.getValue(y[j]);
      }
      LogFile << endl;
      /*for(i = 0;i < n;i++){
	LogFile << "-";
      }
      LogFile << endl;
      for(i = 0;i < n;i++){
	for(j = 0;j < n;j++){
	  LogFile << cplex.getValue(x[i][j]);
	}
	LogFile << endl;
      }/**/
    }
    else {
      LogFile << "No solution found" << endl;
      LogFile << "** Cplex Ends **" << endl << endl;
      Sol = NULL;
    }
    
  }
  catch (IloException& e) {
    cerr << "Concert exception caught: " << e << endl;
    Sol = NULL;
  }
  catch (...) {
    cerr << "Unknown exception caught" << endl;
    Sol = NULL;
  }
  env.end();
  return Sol;
}

void gnuplot_goldberg_
(SQM_instance &Inst,
 int p,
 IloCplex *cplex,
 IloBoolVarArray *x,
 BoolVarArrayMatrix *y
 )
{
  int i,j,k;
  int n = Inst.potential_sites();
  point *puntos = Inst.demand(0);

  char outfilename[32],centersfilename[32];

  sprintf(centersfilename,"../Instancias/Centros_%d.dat",n);
  ofstream centros(centersfilename);

  for(i = 0;i < n;i++){
    if (cplex->getValue((*x)[i]))
      centros << puntos[i].x << " " << puntos[i].y << endl;
  }
  centros.close();

  for (k = 0;k < p;k++) {
    sprintf(outfilename,"../Instancias/Ejes_%d_%d.dat",n,k+1);
    ofstream outfile(outfilename);

    for(i = 0;i < n;i++){
      for(j = 0;j < n;j++){
	if (i != j && cplex->getValue((*y)[i][j][k]))
	  outfile << puntos[i].x << " " 
		  << puntos[i].y << " " 
		  << puntos[j].x - puntos[i].x << " " 
		  << puntos[j].y - puntos[i].y << endl;
      }
    }

    outfile.close();
  }

  FILE *gnuPipe = popen("gnuplot","w");
  fprintf(gnuPipe,"set term svg\n");
  fprintf(gnuPipe,"set output '../gnuplot/Goldberg_%d_%d.svg'\n",n,p);
  fprintf(gnuPipe,"unset key\n");
  fprintf(gnuPipe,"unset border\n");
  fprintf(gnuPipe,"unset yzeroaxis\n");
  fprintf(gnuPipe,"unset xtics\n");
  fprintf(gnuPipe,"unset ytics\n");
  fprintf(gnuPipe,"unset ztics\n");
  //fprintf(gnuPipe,"set title \"Servicio de %.0f\n",cplex->getObjValue());
  //fprintf(gnuPipe,"set style arrow 1 nohead lw 2\n");
  //fprintf(gnuPipe,"set arrow arrowstyle 1\n");
  fprintf(gnuPipe,"plot ");
  fprintf(gnuPipe,"'../Instancias/Q_MCLP_%d.txt' every ::1 using 1:2 with points lc rgb \"black\"",n);
  fprintf(gnuPipe,", '../Instancias/Centros_%d.dat' using 1:2 with points lc rgb \"red\"",n);
  for (k = 0;k < p;k++) {
    fprintf(gnuPipe,", '../Instancias/Ejes_%d_%d.dat' using 1:2:3:4 with vectors nohead" /* linecolor rgb \"dark-blue\"" */ ,n,k+1);
  }

  fprintf(gnuPipe,"\n");
  pclose(gnuPipe);

}
