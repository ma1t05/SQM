/*
 * Modelo que determina ubicaciones de ajustadores
 * 
 */

#include "Goldberg.h"

void Goldberg
(SQM_instance* I, // Set of points
 int p, // facilities
 float mu, // rate parameter
 float f) //
{
  IloEnv env;
  try {
    int i,j,k,r;
    IloInt n = I->N,m= I->M;
    IloNum f_i;
    IloNum rho;
    IloNum M = 10000.0;
    point *client,*potencial_site;
          
    cout << "Comienza definicion del Modelo" << endl;
    IloModel modelo(env);
    
    cout << "++ Variables ++" << endl;
    IloNumVar S(env,0,IloInfinity,IloNumVar::Float,"S");
    IloBoolVarArray x(env);
    BoolVarArrayMatrix y(env,m);
    
    char VarName[16];
    for(i = 0;i < n;i++){
      sprintf(VarName,"x%d",i+1);
      x.add(IloBoolVar(env,VarName));
    }

    cout << "Nombre las variables para facil identificacion" << endl;
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

    cout << "++ Restricciones ++" << endl;
    cout << "Solo una instalacion ocupa la posion k del cliente i" << endl;
    for (i = 0;i < m;i++) {
      for (k = 0;k < p;k++) {
	IloExpr Cover(env);
	for (j = 0;j < n;j++) Cover += y[i][j][k];	
	modelo.add(Cover == 1);
	Cover.end();
      }
    }

    cout << "La instalacion j solo puede ocupar una posicion del cliente i" << endl;
    for (i = 0;i < m;i++) {
      for (j = 0;j < n;j++) {
	modelo.add(IloSum(y[i][j]) <= 1);
      }
    }

    cout << "Instalaciones a abrir" << endl;
    modelo.add(IloSum(x) == p);

    cout << "Relacionar variables de localizacion y asignacion" << endl;
    for(i = 0;i < m;i++){
      for (j = 0;j < n;j++) {
	for (k = 0;k < p;k++) 
	  modelo.add(y[i][j][k] <= x[j]);
      }
    }

    cout << "Restricion de orden de asignaicion" << endl;
    client = I->V;
    potencial_site = I->W;
    NumMatrix O(env,m);
    for (i = 0;i < m;i++)
      O[i] = IloNumArray(env,n);
    for (i = 0;i < m;i++) {
      for (j = 0;j < n;j++) {
	O[i][j] = dist(&(client[i]),&(potencial_site[j]));
      }
    }
      
    for (i = 0;i < m;i++) {
      for (k = 1;k < p;k++) {
	for (j = 0;j < n;j++) {
	  IloExpr balance(env);
	  for (r = 0;r < n;r++) {
	    if (O[i][r] <= O[i][j])
	      balance += y[i][r][k-1];
	  }
	  modelo.add(y[i][j][k] <= balance);
	  balance.end();
	}
      }
    }

    // Cargas de trabajo
    IloNum coef;
    // rho from Daskin
    rho = 0.0;
    for (i = 0;i < m;i++)
      rho += f * client[i].demand;
    rho /= (mu * p);
    cout << "rho = " << rho << endl;
    // rho from ReVelle & Hogan
    // pendiente
    
    for (j = 0;j < n;j++) {
      IloExpr workload(env);
      for (k = 0;k < p;k++) {
	coef = (1 - rho) * pow(rho,k);
	for (i = 0;i < m;i++) {
	  f_i = f * client[i].demand;
	  workload += f_i * coef * O[i][j]* y[i][j][k];
	}
      }
      modelo.add(S + M*x[j] <= workload + M);
    }
    
    cout << "++ Funcion Objetivo ++" << endl;
    modelo.add(IloMaximize(env,S));

    cout << "Resuelve modelo" << endl;
    IloCplex cplex(modelo);
    cplex.setOut(LogFile);
    stringstream ModelName;
    ModelName << "Goldberg-model_" << m << "_" << n << "_" << p << ".lp";
    cplex.exportModel(ModelName.str().c_str());
    if (cplex.solve()) {
      cout << "Solution status: " << cplex.getStatus() << endl;
      cout << "Maximum profit = " << cplex.getObjValue() << endl;
      for (j = 0;j < n;j++) {
	if (cplex.getValue(x[j]) > 0.5) cout << j+1 << " ";
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

void gnuplot_goldberg(SQM_instance *I,int p,IloCplex *cplex, IloBoolVarArray *x, BoolVarArrayMatrix *y) {
  int i,j,k;
  int n = I->N,m = I->M;
  point *client = I->V;
  point *potencial_site = I->W;
  char outfilename[32],centersfilename[32],clientsfilename[32];

  sprintf(centersfilename,"Tmp_centers_%d.dat",n);
  ofstream centros(centersfilename);

  for(j = 0;j < n;j++){
    if (cplex->getValue((*x)[j]))
      centros << potencial_site[j].x << " " << potencial_site[j].y << endl;
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
	cout << "i:" << i << "\tj:" << j << "\tk:" << k << endl;
	if (cplex->getValue((*y)[i][j][k]) > 0.5)
	  outfile << client[i].x << " " 
		  << client[i].y << " " 
		  << potencial_site[j].x - client[i].x << " " 
		  << potencial_site[j].y - client[i].y << endl;
      }
    }

    outfile.close();
  }

  FILE *gnuPipe = popen("gnuplot","w");
  fprintf(gnuPipe,"set term svg\n");
  fprintf(gnuPipe,"unset key\n");
  fprintf(gnuPipe,"unset border\n");
  fprintf(gnuPipe,"unset yzeroaxis\n");
  fprintf(gnuPipe,"unset xtics\n");
  fprintf(gnuPipe,"unset ytics\n");
  fprintf(gnuPipe,"unset ztics\n");
  //fprintf(gnuPipe,"set title \"Servicio de %.0f\n",cplex->getObjValue());
  //fprintf(gnuPipe,"set style arrow 1 nohead lw 2\n");
  //fprintf(gnuPipe,"set arrow arrowstyle 1\n");
  fprintf(gnuPipe,"set output 'Goldberg_%d_%d_%d.svg'\n",m,n,p);
  fprintf(gnuPipe,"plot ");
  fprintf(gnuPipe,"'%s' using 1:2 with points lc rgb \"black\"",clientsfilename);
  fprintf(gnuPipe,", '%s' using 1:2 with points lc rgb \"red\"",centersfilename);
  for (k = 0;k < p;k++) {
    sprintf(outfilename,"Tmp_edges_%d_%d_%d.dat",m,n,k+1);
    fprintf(gnuPipe,", '%s' using 1:2:3:4 with vectors nohead" /* linecolor rgb \"dark-blue\"" */ ,outfilename);
  }
  fprintf(gnuPipe,"\n");

  for (k = 0;k < p;k++) {
    sprintf(outfilename,"Tmp_edges_%d_%d_%d.dat",m,n,k+1);
    fprintf(gnuPipe,"set output 'Goldberg_%d_%d_%d_%d.svg'\n",m,n,p,k+1);
    fprintf(gnuPipe,"plot ");
    fprintf(gnuPipe,"'%s' using 1:2 with points lc rgb \"black\"",clientsfilename);
    fprintf(gnuPipe,", '%s' using 1:2 with points lc rgb \"red\"",centersfilename);
    fprintf(gnuPipe,", '%s' using 1:2:3:4 with vectors nohead" /* linecolor rgb \"dark-blue\"" */ ,outfilename);
    fprintf(gnuPipe,"\n");
  }

  pclose(gnuPipe);

}
