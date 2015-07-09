/*
 * Modelo de Aproxmiacion al problema de Berman
 * 
 */


void SQM_model
(instance* I, // Set of points
 int p, // facilities
 float mu, // rate parameter
 float f) // 
{
  IloEnv env;
  try {
    IloInt i,j,k,r;
    IloInt n = I->n;
    IloNum f_i;
    IloNum rho;
    IloNum M = 10000.0;
    point *puntos = I->points;
    cout << "Comienza definicion del Modelo" << endl;
    IloModel modelo(env);
    
    cout << "++ Variables ++" << endl;
    IloIntVarArray x(env);
    BoolVarArrayMatrix y(env,n);
    BoolVarMatrix u(env,n);
    BoolVarMatrix v(env,n);
    
    cout << "Nombre las variables para facil identificacion" << endl;
    char VarName[16];
    for(i = 0;i < n;i++){
      sprintf(VarName,"x%d",i+1);
      x.add(IloBoolVar(env,VarName));
    }

    for (i = 0;i < n;i++) {
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

    for (i = 0;i < n;i++) {
      IloBoolVarArray u_i(env);
      for (j = 0;j < n;j++) {
	sprintf(VarName,"u_%d_%d",i+1,j+1);
	u_i.add(IloBoolVar(env,VarName));
      }
      u[i] = u_i;
    }

    for (i = 0;i < n;i++) {
      IloBoolVarArray v_i(env);
      for (j = 0;j < n;j++) {
	sprintf(VarName,"v_%d_%d",i+1,j+1);
	v_i.add(IloBoolVar(env,VarName));
      }
      v[i] = v_i;
    }

    cout << "++ Restricciones ++" << endl;
    cout << "Instalaciones a abrir" << endl;
    modelo.add(IloSum(x) == p);

    cout << "Relacionar variables de localizacion y asignacion" << endl;
    for(i = 0;i < n;i++){
      for (j = 0;j < n;j++) {
	for (k = 0;k < p;k++) 
	  modelo.add(y[i][j][k] <= x[j]);
      }
    }

    cout << "Solo una instalacion ocupa la posion k del cliente i" << endl;
    for (i = 0;i < n;i++) {
      for (l = 0;l < k;l++) {
	IloExpr Cover(env);
	for (j = 0;j < n;j++) Cover += y[i][j][l];
	modelo.add(Cover == 1);
	Cover.end();
      }
    }

    cout << "Los ajustadores en j solo puede ocupar x_j posiciones del cliente i" << endl;
    for (i = 0;i < n;i++) {
      for (j = 0;j < n;j++) {
	modelo.add(IloSum(y[i][j]) <= x[j]);
      }
    }

    cout << "Restricion de orden de asignaicion" << endl;
    NumMatrix O(env,n);
    for (i = 0;i < n;i++)
      O[i] = IloNumArray(env,n);
    for (i = 0;i < n;i++) {
      O[i][i] = 0;
      for (j = i+1;j < n;j++) {
	O[i][j] = dist(&(puntos[i]),&(puntos[j]));
	O[j][i] = dist(&(puntos[i]),&(puntos[j]));
      }
    }
      
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
    }

    cout << "Restricciones disjuntas" << endl;
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
	IloExpr Instalaciones(env);
	for (r = 0;r < m;r++) {
	  if (O[i][r] < O[i][j])
	    Instalaciones += x[r];
	}
	modelo.add(Instalaciones <= p - (p - (k - 1)) * v[i][j]);
	modelo.add(Instalaciones + M * v[i][j] >= k + 1);

	/* Si v_{ij} = 1 y u_{ij} = 0 alguntas de las instalaciones de j seran asignadas a i para completar las k */
	modelo.add(IloSum(y[i][j]) + M * (1 - v[i][j] + u[i][j]) >= k - Instalaciones);
	modelo.add(IloSum(y[i][j]) + M * (1 - v[i][j] + u[i][j]) <= k - Instalaciones);

	Instalaciones.end();
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

    // Cargas de trabajo
    IloNum coef;
    // rho from Daskin
    rho = 0.0;
    for (i = 0;i < n;i++)
      rho += f * puntos[i].demand;
    rho /= (mu * p);
    cout << "rho = " << rho << endl;
    // rho from ReVelle & Hogan
    // pendiente
    
    for (j = 0;j < n;j++) {
      IloExpr workload(env);
      for (k = 0;k < p;k++) {
	coef = (1 - rho) * pow(rho,k);
	for (i = 0;i < n;i++) {
	  f_i = f * puntos[i].demand;
	  workload += f_i * coef * O[i][j]* y[i][j][k];
	}
      }
      modelo.add(S + M*x[j] <= workload + M);
    }
    
    cout << "Funcion Objetivo" << endl;
    modelo.add(IloMaximize(env,S));

    cout << "Resuelver modelo" << endl;
    IloCplex cplex(modelo);
    cplex.exportModel("Problema-P.lp");
    if(cplex.solve()){
      cout << "Solution status: " << cplex.getStatus() << endl;
      cout << "Maximum profit = " << cplex.getObjValue() << endl;
      for(j = 0;j < n;j++){
	if(cplex.getValue(x[j]) > 0) cout << j+1 << " ";
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
      gnuplot_goldberg(I,p,&cplex,&x,&y);
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

void gnuplot_goldberg(instance *I,int p,IloCplex *cplex, IloBoolVarArray *x, BoolVarArrayMatrix *y) {
  IloInt i,j,k;
  IloInt n = I->n;
  point *puntos = I->points;

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