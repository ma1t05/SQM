
#ifndef _LOCAL_SEARCH_H
#define _LOCAL_SEARCH_H 1

#include "SQM_Solution.h"

typedef list<int> Sites;
void Local_Search(SQM_solution&);
int LS_get_server_with_less_workload(SQM_solution&);
int LS_get_server_with_more_workload(SQM_solution&);
Sites* LS_get_adjacent_sites(SQM_solution&,int);

#endif

/* eof */
