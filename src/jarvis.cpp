/***************************************************************************
 * Given:
 *  \lambda_m, arrival rate of costumer of type m
 *  \Tao_{im}, expected service time for server i and customer of type m
 *  a_{mk} the kth preferred server for customers of type m
 *     from k = 1,...,C; i = 1,...,N; k = 1,...,N
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

double correction_factor_Q(int N,double rho,int k);
int factorial(int k);

double jarvis_hypercube_approximation
(int C, /* Number of types of customers*/
 int N, /* Number of servers */
 double *lambda, /* */
 double **Tao, /* */
 int **a /* */) {

  double rho;
  double *Q_N_rho;
  Q_N_rho = new double[N];
  for (int k = 0;k < N;k++)
    Q_N_rho[k] = correction_factor_Q(N,rho,k);
  
}

double correction_factor_Q(int N,double rho,int k){
  double Q;
  double P_0,P_N;
  for (int j = k;j < N;j++) {
    Q += (N - j) * pow(N,j) * pow(rho,j - k) / factorial(j - k);
  }
  return Q * P_0 * factorial(N - k - 1) / (pow(1 - P_N,k) * factorial(N) * (1 - rho * (1 - P_N)));
}
					  
int factorial(int n) {
  if (n < 2)
    return 1;
  return n * factorial(n - 1);
}

/* eof */
