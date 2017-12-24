
#include <stdio.h>
#include <stdlib.h>
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
  char * truthfile = argv[2];
  char * outfile = argv[3];

  //sbm_dyn_csr * mat = sbm_read_tsv( "test.tsv" );
  //sbm_dyn_csr * mat = sbm_read_tsv("/home/HGST4TB/graph_challenge_data/synthetic/graph500-scale18-ef16/graph500-scale18-ef16_adj.tsv");

  sbm_dyn_csr * mat = sbm_read_tsv(infile, 1);
  int32_t * blk_truth = sbm_read_truth_tsv(truthfile, mat->nrows);

  double start = monotonic_seconds();
  int32_t * blk_assn = sbm_partition_nblks(mat, 3);
  total_time += monotonic_seconds() - start;

  sbm_partition_metrics * prec_recall = sbm_pairwise_prec_recall(blk_assn, blk_truth, mat->nrows);

  printf("\nNumber of vertices: %d\nNumber of edges: %d\n", mat->nrows, mat->nnz);
  printf("\nPrecision: %0.04f\n", prec_recall->prec);
  printf("Recall: %0.04f\n", prec_recall->recall);

  printf("\nTotal time: %0.04f\n", total_time);
  sbm_write_partitions(outfile, blk_assn, mat->nrows);

  free(blk_truth);
  free(blk_assn);
  free(prec_recall);

  sbm_dyn_csr_free(mat);
}
