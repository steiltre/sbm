
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sbm_file.h"
#include "graph.h"

void sbm_read_tsv_append(
    char const * const fname,
		sbm_dyn_csr * adj_mat,
    const int32_t indexing)
{
  /* Open file */
  FILE * fin;
  if ((fin = fopen(fname, "r")) == NULL) {
    fprintf(stderr, "unable to open '%s' for reading.\n", fname);
    exit(EXIT_FAILURE);
  }

  char * line = malloc(1024 * 1024);
  size_t len = 0;
  ssize_t read = getline(&line, &len, fin);

  int nedges = 0;

  while (read >= 0) {
    nedges++;

    read = getline(&line, &len, fin);
  }

  edge * edges = malloc(nedges * sizeof(*edges));

  fin = fopen(fname, "r");
  read = getline(&line, &len, fin);

  int32_t i = 0;
  while (read >= 0) {
    char * ptr = strtok(line, "\t");
    char * end = NULL;

    int32_t row = strtol(ptr, &end, 10) - indexing;


    ptr = strtok(NULL, "\t");
    end = NULL;

    int32_t col = strtol(ptr, &end, 10) - indexing;

    ptr = strtok(NULL, "\t");

    int32_t val = strtol(ptr, &end, 10);

    (edges+i)->from = row;
    (edges+i)->to = col;
    (edges+i)->wgt = val;

    read = getline(&line, &len, fin);
    i++;
  }

  qsort( edges, nedges, sizeof(*edges), sbm_lt_from );

  for (int e=0; e<nedges; e++) {
    sbm_dyn_csr_add(adj_mat, (edges+e)->from, (edges+e)->to, (edges+e)->wgt);
  }

  sbm_dyn_csr_flush_adj_lists(adj_mat);

  fclose(fin);

  free(line);
  free(edges);
}

sbm_dyn_csr * sbm_read_tsv(
		char const * const fname,
		const int32_t indexing)
{
  sbm_dyn_csr * mat = sbm_dyn_csr_alloc();

	sbm_read_tsv_append(fname, mat, indexing);

	return mat;
}

int32_t sbm_count_nodes_tsv(
		char const * const fname,
		const int32_t indexing)
{
  /* Open file */
  FILE * fin;
  if ((fin = fopen(fname, "r")) == NULL) {
    fprintf(stderr, "unable to open '%s' for reading.\n", fname);
    exit(EXIT_FAILURE);
  }

  char * line = malloc(1024 * 1024);
  size_t len = 0;
  ssize_t read = getline(&line, &len, fin);

  int32_t max_lbl = 0;

  while (read >= 0) {
    char * ptr = strtok(line, "\t");
    char * end = NULL;

    int32_t row = strtol(ptr, &end, 10);


    ptr = strtok(NULL, "\t");
    end = NULL;

    int32_t col = strtol(ptr, &end, 10);

		if (row > max_lbl) {
			max_lbl = row;
		}
		if (col > max_lbl) {
			max_lbl = col;
		}
  	read = getline(&line, &len, fin);
  }

	fclose(fin);

	free(line);

	return max_lbl - indexing + 1;    /* Think of # of rows with 0-indexing to see necessity of +1 */
}

int32_t * sbm_read_truth_tsv(
    char const * const fname,
    const int32_t nnodes,
		const int32_t indexing)
{
  /* Open file */
  FILE * fin;
  if ((fin = fopen(fname, "r")) == NULL) {
    fprintf(stderr, "unable to open '%s' for reading.\n", fname);
    exit(EXIT_FAILURE);
  }

  int32_t * truth = malloc(nnodes * sizeof(*truth));

  char * line = malloc(1024 * 1024);
  size_t len = 0;
  ssize_t read = getline(&line, &len, fin);

  while (read >= 0) {
    char * ptr = strtok(line, "\t");
    char * end = NULL;

		int32_t row = strtol(ptr, &end, 10) - indexing;

    ptr = strtok(NULL, "\t");
    truth[row] = strtol(ptr, &end, 10);

    read = getline(&line, &len, fin);
  }

  fclose(fin);

  free(line);

  return truth;
}

void sbm_write_partitions(
    char const * const fname,
    const int32_t * const blk_assn,
    const int32_t nnodes)
{
  FILE * fout = fopen(fname, "w");

  for (int i=0; i<nnodes; i++) {
    fprintf(fout, "%d\n", blk_assn[i]);
  }
  fclose(fout);
}
