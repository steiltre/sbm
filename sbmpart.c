
#include <stdio.h>
#include <stdlib.h>
#include "omp.h"
#include "sbm_file.h"
#include "sbm_dyn_csr.h"
#include "sbm_partition.h"
#include "sbm_util.h"

double total_time = 0;
/*
double flushing_time = 0;
double block_merge_time = 0;
double nodal_update_time = 0;
*/

int main(
    int argc,
    char ** argv)
{
  char * infile = argv[1];
  int nfiles = atoi(argv[2]);

  sbm_stream_partition(infile, nfiles, 3.0);
}
