
#ifndef _SQM_H
#define _SQM_H

#include <ctime>

using namespace std;
#include "config.h"
#include "Goldberg.h"
#include "SQM_model.h"
#include "SQM_GRASP.h"
#include "PathRelinking.h"
#include "Local_Search.h"

#define BETA 1.5

void Log_Start_SQMH(int M_clients,int N_sites,int p,double mu,double f);
void Call_SQM_model(SQM_instance*,int,int,double,double,double,string);
void Call_SQM_heuristic(SQM_instance* I,int p,double f,double mu,double v);
void Call_SQM_GRASP(SQM_instance *I,int p,double lambda,double Mu_NT,double v);
void Call_SQM_random(SQM_instance *I,int p,double lambda,double Mu_NT,double v);
void Call_SQM_Path_Relinking(SQM_instance *I,int p,double lambda,double Mu_NT,double v);
void Call_SQM_Local_Search(SQM_instance *I,int p,double lambda,double Mu_NT,double v);

#endif

/* eof */

