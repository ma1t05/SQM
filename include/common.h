
#ifndef _COMMON_H
#define _COMMON_H

#include <ctime>
#include <iostream>
#include <unistd.h> /* getopt */
#include <getopt.h> /* getopt_long */

using namespace std;
#include "config.h"
#include "Goldberg.h"
#include "SQM_model.h"
#include "SQM_GRASP.h"
#include "PathRelinking.h"
#include "Reference_Set.h"
#include "Local_Search.h"

extern string Instance_Name;
extern int M_clients,N_sites;
extern int p,l;
extern double v;

void print_usage ();
void process_command_line(int,char**);
void read_config_file(string configFile);
void Log_Start_SQMH(int M_clients,int N_sites,int p,double mu,double f);
typedef list<SQM_solution*> Solutions;

/* Commands */
typedef void (command)(SQM_instance&, /* The instance*/
		       int,           /* The number of servers */
		       double);       /* The speed */
command Test_SQM_model;
command Test_SQM_multistart;
command Test_SQM_heuristic;
command Test_SQM_GRASP;
command Test_SQM_random;
command Test_SQM_Path_Relinking;
command Test_SQM_Local_Search;
command Tune_Path_Relinking;
/*  */
extern command *Test_Function;

/* Global variables read from config */
extern double BETA;

#endif

/* eof */
