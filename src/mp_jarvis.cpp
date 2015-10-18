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
 * (4) Else compute P__N, \Tao, and f_{im}
 * (5) Return to step 1
 **************************************************************************/

#include "mp_jarvis.h"

void correction_factor_Q(mpf_t Q,int N,mpf_t rho,int k);
unsigned long int factorial(unsigned long int n);
mpf_t P__0,P__N;

void jarvis_hypercube_approximation
(int C, /* Number of types of customers*/
 int N, /* Number of servers */
 mpf_t *lambda, /* arrive rate according to a Poisson process per type */
 mpf_t **Tao, /* expected service time for service i and customer of node m */
 int **a, /* for customers of type m, the list of preferred servers */
 mpf_t **f) {

  mpf_t Lambda;
  mpf_t Rho;
  mpf_t tao;
  mpf_t *rho;
  mpf_t *new_rho;
  mpf_t *Q_N_rho;

  mpf_t max_change;
  mpf_t tmp;

  mpf_init(Lambda);
  rho = new mpf_t [N];
  for (int i = 0; i < N;i++)
    mpf_init(rho[i]);
  mpf_init(tao);
  mpf_init_set_ui(P__0,1);
  mpf_init_set_ui(P__N,1);
  Q_N_rho = new mpf_t [N];
  for (int i = 0;i < N;i++)
    mpf_init(Q_N_rho[i]);
  new_rho = new mpf_t [N];
  for (int i = 0;i < N;i++)
    mpf_init(new_rho[i]);
  mpf_init(Rho);
  mpf_init(max_change);
  mpf_init(tmp);

  /* INITIALIZE: */
  for (int m = 0;m < C;m++) 
    mpf_add(Lambda,Lambda,lambda[m]);

  /* Busy probability rho_i */
  for (int i = 0; i < N;i++) {
    for (int m = 0; m < C;m++) {
      if (a[m][0] == i) {
	mpf_mul(tmp,lambda[m],Tao[i][m]);
	mpf_add(rho[i],rho[i],tmp);
      }
    }
  }

  /* mean service time 'tao' */
  for (int m = 0;m < C;m++) {
    mpf_mul(tmp,lambda[m],Tao[a[m][0]][m]);
    mpf_add(tao,tao,tmp);
  }
  mpf_div(tao,tao,Lambda);

  /* Define intial values for P__0 and P__N */
  for (int i = 0;i < N;i++) {
    mpf_ui_sub(tmp,1,rho[i]);
    mpf_mul(P__0,P__0,tmp);
  }

  for (int i = 0;i < N;i++)
    mpf_mul(P__N,P__N,rho[i]);

  /* ITERATION: */
  do {

    cout << "'rho_i':";
    for (int i = 0;i < N;i++) cout << " " << mpf_get_d(rho[i]);
    cout << endl;

    /* traffic intensity */
    mpf_mul(Rho,Lambda,tao);
    mpf_div_ui(Rho,Rho,N);

    /* Compute Q(N,'rho',k) */
    for (int k = 0;k < N;k++)
      correction_factor_Q(Q_N_rho[k],N,Rho,k);
  
    /* Aproximation of 'rho'_i */
    mpf_t Vi,rho_a_ml;
    mpf_init(Vi);
    mpf_init(rho_a_ml);
    for (int i = 0;i < N;i++) {
      mpf_set_ui(Vi,0);
      for (int m = 0;m < C;m++) {
	mpf_set_ui(rho_a_ml,1);
	for (int k = 0;k < N;k++) {
	  if (a[m][k] == i) {
	    mpf_mul(tmp,lambda[m],Tao[i][m]);
	    mpf_mul(tmp,tmp,Q_N_rho[k]);
	    mpf_mul(tmp,tmp,rho_a_ml);
	    mpf_add(Vi,Vi,tmp);
	    /* Vi += lambda[m] * Tao[i][m] * Q_N_rho[k] * rho_a_ml */
	    break;
	  }
	  mpf_mul(rho_a_ml,rho_a_ml,rho[a[m][k]]);
	}
      }
      mpf_add_ui(tmp,Vi,1);
      mpf_div(new_rho[i],Vi,tmp);
    }
    mpf_clear(rho_a_ml);
    mpf_clear(Vi);

    /* Convergence criterion */
    mpf_set_ui(max_change,0);
    for (int i = 0;i < N;i++) {
      mpf_sub(tmp,rho[i],new_rho[i]);
      mpf_abs(tmp,tmp);
      if (mpf_cmp(tmp,max_change) > 0)
	mpf_set(max_change,tmp);
    }
    
    cout << "max change = " << mpf_get_d(max_change);
    if (mpf_cmp_d(max_change,epsilon) < 0)
      cout << " **STOP**" << endl;
    else 
      cout << "\r";

    if (mpf_cmp_d(max_change,epsilon) < 0) break; /* STOP */

    for (int i = 0;i < N;i++)
      mpf_set(rho[i],new_rho[i]);

    /* Compute P__0 */
    mpf_set_ui(P__0,1);
    for (int i = 0;i < N;i++) {
      mpf_ui_sub(tmp,1,rho[i]);
      mpf_mul(P__0,P__0,tmp);
    }

    /* Compute P__N */
    mpf_t s_rho;
    mpf_init(s_rho);
    for (int i = 0;i < N;i++) 
      mpf_add(s_rho,s_rho,rho[i]);
    mpf_div_ui(tmp,s_rho,N);
    mpf_div(tmp,tmp,Rho);
    mpf_ui_sub(P__N,1,tmp);

    /* Compute f_{im} */
    mpf_init(rho_a_ml);
    for (int i = 0;i < N;i++) {
      for (int m = 0;m < C;m++) {
	mpf_set_ui(rho_a_ml,1);
	for (int k = 0;k < N;k++) {
	  if (a[m][k] == i) {
	    mpf_ui_sub(tmp,1,rho[i]);
	    mpf_mul(tmp,tmp,rho_a_ml);
	    mpf_mul(f[i][m],tmp,Q_N_rho[k]);
	    break;
	  }
	  mpf_mul(rho_a_ml,rho_a_ml,rho[a[m][k]]);
	}
      }
    }
    mpf_clear(rho_a_ml);

    /* Compute mean service time 'tao' */
    mpf_t aux;
    mpf_set_ui(tao,0);
    mpf_init(aux);
    for (int m = 0;m < C;m++) {
      mpf_set_ui(tmp,0);
      for (int i = 0;i < N;i++) {
	mpf_mul(aux,Tao[i][m],f[i][m]);
	mpf_add(tmp,tmp,aux);
      }
      mpf_mul(tmp,lambda[m],tmp);
      mpf_add(tao,tao,tmp);
    }
    mpf_div(tao,tao,Lambda); // / Lambda
    mpf_ui_sub(aux,1,P__N); // 1 - P__N
    mpf_div(tao,tao,aux); // / (1 - P__N)
    mpf_clear(aux);
 
  } while (1);

  /*
    cout << "finsh jarvis" << endl;
    for (int i = 0;i < N;i++) {
    for (int m = 0;m < C;m++) {
    cout << f[i][m] << " ";
    }
    cout << endl;
    }
  */

  mpf_clear(tmp);
  mpf_clear(max_change);
  mpf_clear(Rho);
  for (int i = 0;i < N;i++)
    mpf_clear(new_rho[i]);
  delete [] new_rho;
  for (int i = 0;i < N;i++)
    mpf_clear(Q_N_rho[i]);
  delete [] Q_N_rho;
  mpf_clear(P__N);
  mpf_clear(P__0);
  mpf_clear(tao);
  for (int i = 0;i < N;i++)
    mpf_clear(rho[i]);
  delete [] rho;
  mpf_clear(Lambda);
}

/* 
   Using equation (3)
   Q(N,rho,k) = \sum_{j = k}^{N-1}{\frac{(N-j)(N^j)(\rho^{j-k})P__0(N-k-1)!}{(j-k)!(1-P__N)^kN!(1-\rho(1-P__N))}}
   P__0 probability all servers are idle
   P__N probability all servers are busy
   \rho(1-P__N) average server workload or utilization (Sevast'yonov 1957)
*/
			    
void correction_factor_Q(mpf_t Q, int N,mpf_t rho,int k){
  unsigned long int f;
  mpf_t tmp,tmp1,tmp2;
  /*
  if (PNK_N_1_rho_1_PN == 0) 
    cout << "Algun factor es 0 N = " << N 
	 << " rho = " << rho 
	 << " k = " << k 
	 << " P__0 = " << P__0
	 << " P__N = " << P__N << endl;
  */
  mpf_init(tmp);
  mpf_init_set_d(tmp1,pow(N,k-1)); // N^(k-1)
  f = 1;
  mpf_set_ui(Q,0);
  for (unsigned long int j = k;j < N;j++) {
    /*Q += (N - j) * pow(N,j) * pow(rho,j - k) * P__0 * N_k_1 / (factorial(j - k) * PNK_N_1_rho_1_PN); */
    /* Q += [(N - j) * N^j * rho^(j - k) / (j - k)!] * P__0 * N_k_1 / PNK_N_1_rho_1_PN; */
    mpf_mul_ui(tmp1,tmp1,N); // * N
    if (j - k > 1) mpf_div_ui(tmp1,tmp1,j - k); // / (j - k)
    mpf_mul_ui(tmp,tmp1,N - j); // * (N - j)
    mpf_add(Q,Q,tmp); // Q += (N - j) * N^j * rho^(j - k) / (j - k)!
    mpf_mul(tmp1,tmp1,rho); // * rho
  }
  mpf_init(tmp2);
  mpf_ui_sub(tmp2,1,P__N); // 1 - P__N
  mpf_pow_ui(tmp1,tmp2,k); // (1 - P__N)^k
  mpf_mul(tmp2,tmp2,rho); // rho * (1 - P__N)
  mpf_ui_sub(tmp2,1,tmp2); // 1 - rho * (1 - P__N)
  
  mpf_mul(Q,Q,P__0); // * P__0
  mpf_mul_ui(Q,Q,factorial(N-k-1)); // * (N - k - 1)!

  /* pow(1 - P__N,k) * factorial(N) * (1 - rho * (1 - P__N)) */
  mpf_mul(tmp2,tmp1,tmp2); // (1 - P__N)^k * (1 - rho * (1 - P__N))
  mpf_mul_ui(tmp2,tmp2,factorial(N)); // (1 - P__N)^k * (1 - rho * (1 - P__N)) * N!
  mpf_div(Q,Q,tmp2); // / [(1 - P__N)^k * (1 - rho * (1 - P__N)) * N!]

  mpf_clear(tmp);
  mpf_clear(tmp1);
  mpf_clear(tmp2);
}
			    
unsigned long int factorial(unsigned long int n) {
  unsigned long int N = 1;
  for (unsigned long int i = 2;i <= n;i++)
    N *= i;
  return N;
}

/* eof */
