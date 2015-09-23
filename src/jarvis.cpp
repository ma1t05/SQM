/***************************************************************************
 * Sourece code developed by Luis Maltos
 * mailto: maltosla@gmail.com
 * Implementation of the:
 * Approximating the equilibrium behaivor of multi-server loss systems
 * J. P. Jarvis, Managment Science, February 1985
 **************************************************************************/

/***************************************************************************
 * Given:
 *  \lambda_m, arrival rate of costumer of type m
 *  \Tao_{im}, expected service time for server i and customer of type m
 *  a_{mk} the kth preferred server for customers of type m
 *     from k = 1,...,C; i = 1,...,N; k = 1,...,N
 * Returns:
 *  f_{ij}, matrix of probabilities of i attends j
 * Initialize:
 *  \rho_i = \sum_{m:a_{m1=i}{\lambda_m * \Tao_{im}}
 *  \Tao = \sum_{m=1}^{C}{(lambda_m/lambda)*\Tao_{a_{m1},m}
 * Iteration:
 * (1) Compute Q(N,\rho,k) for k = 0,1,...,N-1 
 *   where \rho = \lambda * \tao / N
 * (2) For i = 1,...,N, the new aproximaion of
 *   \rho_i = V_i / (1 + V_i)
 * (3) Stop if max change in \rho_i is less than convergence criterion.
 * (4) Else compute P_N, \Tao, and f_{im}
 * (5) Return to step 1
 **************************************************************************/

#include "jarvis.h"

num correction_factor_Q(int N,num rho,int k);
int factorial(int k);
num P_0,P_N; /* Pendiente: encontrar valores de P_0 y P_N */

num** jarvis_hypercube_approximation
(int C, /* Number of types of customers*/
 int N, /* Number of servers */
 num *lambda, /* arrive rate according to a Poisson process per type */
 num **Tao, /* expected service time for service i and customer of node m */
 int **a, /* for customers of type m, the list of preferred servers */
 num **f) {

  num Lambda;
  num Rho;
  num tao;
  num *rho;
  num *Q_N_rho;
  num *new_rho;

  /* INITIALIZE: */
  Lambda = 0.0;
  for (int m = 0;m < C;m++) Lambda += lambda[m];

  rho = new num[N];
  for (int i = 0; i < N;i++) {
    rho[i] = 0.0;
    for (int m = 0; m < C;m++) {
      if (a[m][0] == i) rho[i] += lambda[m] * Tao[i][m];
    }
  }

  /* mean service time 'tao' */
  tao = 0.0;
  for (int m = 0;m < C;m++) tao += lambda[m] * Tao[a[m][0]][m];
  tao /= Lambda;

  /* Define intial values for P_0 and P_N */
  P_0 = 1.0;
  for (int i = 0;i < N;i++)
    P_0 *= (1 - rho[i]);

  P_N = 1.0;
  for (int i = 0;i < N;i++)
    P_N *= rho[i];

  /* ITERATION: */
  Q_N_rho = new num[N];
  new_rho = new num[N];
  do {

    /*
    cout << "'rho_i':";
    for (int i = 0;i < N;i++) cout << " " << rho[i];
    cout << "\r";
    */

    /* traffic intensity */
    Rho = Lambda * tao / N;

    /* Compute Q(N,'rho',k) */
    for (int k = 0;k < N;k++)
      Q_N_rho[k] = correction_factor_Q(N,Rho,k);
  
    /* Aproximation of 'rho'_i */
    for (int i = 0;i < N;i++) {
      num Vi = 0.0;
      for (int m = 0;m < C;m++) {
	num rho_a_ml = 1.0;
	int k;
	for (k = 0;k < N;k++) {
	  if (a[m][k] == i)
	    break;
	  rho_a_ml *= rho[a[m][k]];
	}
	Vi += lambda[m] * Tao[i][m] * Q_N_rho[k] * rho_a_ml;
      }
      new_rho[i] = Vi / (1 + Vi);
    }

    /* Convergence criterion */
    num max_change = 0.0;
    for (int i = 0;i < N;i++)
      if (abs(rho[i] - new_rho[i]) > max_change)
	max_change = abs(rho[i] - new_rho[i]);
    /*
    cout << "max change = " << max_change;
    if (max_change < epsilon)
      cout << " **STOP**" << endl;
    else 
      cout << "\r";
    */
    if (max_change < epsilon) break; /* STOP */

    for (int i = 0;i < N;i++)
      rho[i] = new_rho[i];

    /* Compute P_0 */
    P_0 = 1.0;
    for (int i = 0;i < N;i++) P_0 *= (1 - rho[i]);

    /* Compute P_N */
    num s_rho = 0.0;
    for (int i = 0;i < N;i++) s_rho += rho[i];
    P_N = 1.0 - s_rho / (N * Rho);

    /* Compute mean service time 'tao' */
    tao = 0.0;
    for (int m = 0;m < C;m++) {
      num tmp = 0.0;
      for (int i = 0;i < N;i++)
	tmp += Tao[i][m] * f[i][m];
      tao += lambda[m] * tmp / (Lambda * (1 - P_N));
    }
 
    /* Compute f_{im} */
    for (int i = 0;i < N;i++) {
      for (int m = 0;m < C;m++) {
	num rho_a_ml = 1.0;
	int k;
	for (k = 0;k < N;k++) {
	  if (a[m][k] == i)
	    break;
	  rho_a_ml *= rho[a[m][k]];
	}
	f[i][m] = Q_N_rho[k] * (1 - rho[i]) * rho_a_ml;
      }
    }

  } while (1);
  delete [] new_rho;
  delete [] Q_N_rho;

  /*
  cout << "finsh jarvis" << endl;
  for (int i = 0;i < N;i++) {
    for (int m = 0;m < C;m++) {
      cout << f[i][m] << " ";
    }
    cout << endl;
  }
  */
  return f;
}

/* 
   Using equation (3)
   Q(N,rho,k) = \sum_{j = k}^{N-1}{\frac{(N-j)(N^j)(\rho^{j-k})P_0(N-k-1)!}{(j-k)!(1-P_N)^kN!(1-\rho(1-P_N))}}
   P_0 probability all servers are idle
   P_N probability all servers are busy
   \rho(1-P_N) average server workload or utilization (Sevast'yonov 1957)
*/
num correction_factor_Q(int N,num rho,int k){
  num Q;
  num N_k_1 = factorial(N - k - 1);
  num PNK_N_1_rho_1_PN = pow(1 - P_N,k) * factorial(N) * (1 - rho * (1 - P_N));
  if (PNK_N_1_rho_1_PN == 0) 
    cout << "Algun factor es 0 N = " << N 
	 << " rho = " << rho 
	 << " k = " << k 
	 << " P_0 = " << P_0
	 << " P_N = " << P_N << endl;
  Q = 0.0;
  for (int j = k;j < N;j++) {
    Q += (N - j) * pow(N,j) * pow(rho,j - k) * P_0 * N_k_1 / (factorial(j - k) * PNK_N_1_rho_1_PN);
  }
  return Q;
}
					  
int factorial(int n) {
  if (n < 2)
    return 1;
  return n * factorial(n - 1);
}

/* eof */
