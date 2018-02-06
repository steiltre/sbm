
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sbm_dyn_csr.h"
#include "sbm_util.h"

sbm_adj_list_node * sbm_adj_list_node_init()
{
  sbm_adj_list_node * node = malloc(sizeof(*node));

  node->ind = -1;
  node->deg = -1;
  node->next = NULL;

  return node;
}

void sbm_adj_list_add_node(
    int32_t ind,
    int32_t deg,
    sbm_adj_list_node * prev)
{
  sbm_adj_list_node * node = sbm_adj_list_node_init();

  node->ind = ind;
  node->deg = deg;
  node->next = prev->next;
  prev->next = node;
}

sbm_dyn_csr * sbm_dyn_csr_alloc()
{
  sbm_dyn_csr * mat = malloc( sizeof(*mat) );

  mat->row_ptr = malloc( sizeof(*mat->row_ptr) );
  mat->row_ind = NULL;
  mat->row_val = NULL;

  mat->col_ptr = malloc( sizeof(*mat->col_ptr) );
  mat->col_ind = NULL;
  mat->col_val = NULL;

  mat->row_adj_lists = malloc(sizeof(*mat->row_adj_lists));
  mat->col_adj_lists = malloc(sizeof(*mat->col_adj_lists));

  mat->nrows = 0;
  mat->ncols = 0;
  mat->nnz = 0;
  mat->adj_list_size = 0;
  mat->row_adj_list_new_nodes = 0;
  mat->col_adj_list_new_nodes = 0;

  mat->row_ptr[0] = 0;
  mat->col_ptr[0] = 0;

  return mat;
}

void sbm_dyn_csr_free(
    sbm_dyn_csr * mat)
{
  if (mat->row_ptr != NULL) {
    free(mat->row_ptr);
  }
  if (mat->row_ind != NULL) {
    free(mat->row_ind);
  }
  if (mat->row_val != NULL) {
    free(mat->row_val);
  }
  if (mat->col_ptr != NULL) {
    free(mat->col_ptr);
  }
  if (mat->col_ind != NULL) {
    free(mat->col_ind);
  }
  if (mat->col_val != NULL) {
    free(mat->col_val);
  }

  for (int32_t i=0; i<mat->nrows; i++) {
    sbm_adj_list_node * prev = mat->row_adj_lists[i];
    sbm_adj_list_node * curr = mat->row_adj_lists[i]->next;
    while (curr != NULL) {
      free(prev);
      prev = curr;
    }
    free(prev);
  }
  free(mat->row_adj_lists);

  for (int32_t i=0; i<mat->ncols; i++) {
    sbm_adj_list_node * prev = mat->col_adj_lists[i];
    sbm_adj_list_node * curr = mat->col_adj_lists[i]->next;
    while (curr != NULL) {
      free(prev);
      prev = curr;
    }
    free(prev);
  }
  free(mat->col_adj_lists);

  free(mat);
}

/* Assumes adjacency lists have been flushed */
sbm_dyn_csr * sbm_dyn_csr_deep_cpy(
    sbm_dyn_csr * mat)
{
  sbm_dyn_csr * ret = malloc(sizeof(*ret));;

  ret->row_ptr = malloc((mat->nrows+1)*sizeof(*ret->row_ptr));
  memcpy(ret->row_ptr, mat->row_ptr, (mat->nrows+1)*sizeof(*ret->row_ptr));
  ret->col_ptr = malloc((mat->ncols+1)*sizeof(*ret->col_ptr));
  memcpy(ret->col_ptr, mat->col_ptr, (mat->ncols+1)*sizeof(*ret->col_ptr));

  ret->row_ind = malloc(mat->nnz*sizeof(*ret->row_ind));
  memcpy(ret->row_ind, mat->row_ind, mat->nnz*sizeof(*ret->row_ind));
  ret->col_ind = malloc(mat->nnz*sizeof(*ret->col_ind));
  memcpy(ret->col_ind, mat->col_ind, mat->nnz*sizeof(*ret->col_ind));

  ret->row_val = malloc(mat->nnz*sizeof(*ret->row_val));
  memcpy(ret->row_val, mat->row_val, mat->nnz*sizeof(*ret->row_val));
  ret->col_val = malloc(mat->nnz*sizeof(*ret->col_val));
  memcpy(ret->col_val, mat->col_val, mat->nnz*sizeof(*ret->col_val));

  ret->row_adj_lists = malloc(mat->nrows*sizeof(*ret->row_adj_lists));
  for (int32_t i=0; i<mat->nrows; i++) {
    ret->row_adj_lists[i] = sbm_adj_list_node_init();
  }
  ret->col_adj_lists = malloc(mat->nrows*sizeof(*ret->col_adj_lists));
  for (int32_t i=0; i<mat->ncols; i++) {
    ret->col_adj_lists[i] = sbm_adj_list_node_init();
  }

  ret->nrows = mat->nrows;
  ret->ncols = mat->ncols;
  ret->nnz = mat->nnz;
  ret->adj_list_size = mat->adj_list_size;
  ret->row_adj_list_new_nodes = mat->row_adj_list_new_nodes;
  ret->col_adj_list_new_nodes = mat->col_adj_list_new_nodes;

  return ret;
}

void sbm_dyn_csr_add(
    sbm_dyn_csr * mat,
    int32_t row,
    int32_t col,
    int32_t val)
{
	/* Determine largest index (row or col) involved to make sure matrix remains square */
	int32_t max_row_col;
	if (row > col) {
		max_row_col = row;
	}
	else {
		max_row_col = col;
	}

  /* First add entry to row representation */
  if (max_row_col >= mat->nrows + mat->row_adj_list_new_nodes) {
    mat->row_adj_lists = (sbm_adj_list_node **) realloc(mat->row_adj_lists, (max_row_col+1) * sizeof(*mat->row_adj_lists));
    for (int32_t i = mat->nrows + mat->row_adj_list_new_nodes; i<max_row_col+1; i++) {
      mat->row_adj_lists[i] = sbm_adj_list_node_init();
    }
    mat->row_adj_list_new_nodes = max_row_col - mat->nrows + 1;    /* +1 accounts for 0-indexing */
  }
  sbm_adj_list_node * curr = mat->row_adj_lists[row];

  if (curr->next != NULL) {
    while (curr->next->ind <= col) {    /* Sort adjacency list in ascending order for ease of flushing */
      curr = curr->next;

      if (curr->next == NULL) {
        break;
      }
    }
  }

	/* Updates values by adding to current value. Replacing may be desired. */
  if (curr->ind == col) {     /* Same entry in matrix */
    curr->deg += val;
  }
  else {
    sbm_adj_list_add_node(col, val, curr);
  }

  /* Then add entry to column representation */
  if (max_row_col >= mat->ncols + mat->col_adj_list_new_nodes) {
    mat->col_adj_lists = (sbm_adj_list_node **) realloc(mat->col_adj_lists, (max_row_col+1) * sizeof(*mat->col_adj_lists));
    for (int32_t i = mat->ncols + mat->col_adj_list_new_nodes; i<max_row_col+1; i++) {
      mat->col_adj_lists[i] = sbm_adj_list_node_init();
    }
    mat->col_adj_list_new_nodes = max_row_col - mat->ncols + 1;    /* +1 accounts for 0-indexing */
  }
  curr = mat->col_adj_lists[col];

  if (curr->next != NULL) {
    while (curr->next->ind <= row) {    /* Sort adjacency list in ascending order for ease of flushing */
      curr = curr->next;

      if (curr->next == NULL) {
        break;
      }
    }
  }

	/* Updates values by adding to current value. Replacing may be desired. */
  if (curr->ind == row) {     /* Same entry in matrix */
    curr->deg += val;
  }
  else {
    sbm_adj_list_add_node(row, val, curr);
    mat->adj_list_size++;
  }
}

void sbm_dyn_csr_flush_adj_lists(
    sbm_dyn_csr * mat)
{
  /* Allocate new csr arrays to fill */
  int32_t * new_row_ind = malloc( (mat->nnz + mat->adj_list_size) * sizeof(*new_row_ind) );
  int32_t * new_row_val = malloc( (mat->nnz + mat->adj_list_size) * sizeof(*new_row_val) );
  int32_t * new_col_ind = malloc( (mat->nnz + mat->adj_list_size) * sizeof(*new_col_ind) );
  int32_t * new_col_val = malloc( (mat->nnz + mat->adj_list_size) * sizeof(*new_col_val) );

  mat->row_ptr = (int32_t *) realloc(mat->row_ptr, (mat->nrows + mat->row_adj_list_new_nodes + 1) * sizeof(*mat->row_ptr));
  for (int32_t i=0; i<mat->row_adj_list_new_nodes; i++) {
    mat->row_ptr[ mat->nrows + i + 1] = mat->row_ptr[mat->nrows];
  }
  mat->col_ptr = (int32_t *) realloc(mat->col_ptr, (mat->ncols + mat->col_adj_list_new_nodes + 1) * sizeof(*mat->col_ptr));
  for (int32_t i=0; i<mat->col_adj_list_new_nodes; i++) {
    mat->col_ptr[ mat->ncols + i + 1] = mat->col_ptr[mat->ncols];
  }
  mat->nrows += mat->row_adj_list_new_nodes;
  mat->ncols += mat->col_adj_list_new_nodes;
  mat->row_adj_list_new_nodes = 0;
  mat->col_adj_list_new_nodes = 0;
  mat->adj_list_size = 0;

  /* Flush adjacency list for row representation */
  int32_t edges_processed = 0;
  for (int32_t i=0; i<mat->nrows; i++) {
    sbm_adj_list_node * prev;
    sbm_adj_list_node * curr = mat->row_adj_lists[i]->next;
    mat->row_adj_lists[i]->next = NULL;

    int32_t prev_row_ptr = mat->row_ptr[i];
    mat->row_ptr[i] = edges_processed;

    for (int32_t j=prev_row_ptr; j<mat->row_ptr[i+1]; j++) {     /* Walk along row of old csr */
      if (curr == NULL) {    /* No more entries in adjacency list. Copy rest of original row */
        memmove(new_row_val + edges_processed, mat->row_val + j, (mat->row_ptr[i+1] - j)*sizeof(*mat->row_val));
        memmove(new_row_ind + edges_processed, mat->row_ind + j, (mat->row_ptr[i+1] - j)*sizeof(*mat->row_ind));
        edges_processed += mat->row_ptr[i+1] - j;
        break;
      }

      while (curr->ind < mat->row_ind[j]) {
        new_row_val[edges_processed] = curr->deg;
        new_row_ind[edges_processed++] = curr->ind;
        prev = curr;
        curr = curr->next;
        free(prev);
        if (curr == NULL) {
          break;
        }
      }

      if (curr != NULL) {
        while (curr->ind == mat->row_ind[j]) {   /* Add degrees to value in original array when row_ind already exists */
          mat->row_val[j] += curr->deg;
          prev = curr;
          curr = curr->next;
          free(prev);
          if (curr == NULL) {
            break;
          }
        }
      }

      if (mat->row_val[j] != 0) {    /* Don't want to add values that have become zero after adding adjacency list */
        new_row_val[edges_processed] = mat->row_val[j];
        new_row_ind[edges_processed++] = mat->row_ind[j];
      }
    }
    /* All values from original row have been added. Walk to end of adjacency list */
    int32_t prev_row_ind = -1;
    while (curr != NULL) {
      if (curr->ind == prev_row_ind) {
        new_row_val[edges_processed-1] += curr->deg;
      }
      else {
        new_row_val[edges_processed] = curr->deg;
        new_row_ind[edges_processed++] = curr->ind;
      }
      prev_row_ind = curr->ind;
      sbm_adj_list_node * prev = curr;
      curr = curr->next;
      free(prev);
    }
  }
  mat->row_ptr[mat->nrows] = edges_processed;

  /* Flush adjacency list for column representation */
  edges_processed = 0;
  for (int32_t i=0; i<mat->ncols; i++) {
    sbm_adj_list_node * prev;
    sbm_adj_list_node * curr = mat->col_adj_lists[i]->next;
    mat->col_adj_lists[i]->next = NULL;

    int32_t prev_col_ptr = mat->col_ptr[i];
    mat->col_ptr[i] = edges_processed;

    for (int32_t j=prev_col_ptr; j<mat->col_ptr[i+1]; j++) {     /* Walk along col of old csr */
      if (curr == NULL) {    /* No more entries in adjacency list. Copy rest of original col */
        memmove(new_col_val + edges_processed, mat->col_val + j, (mat->col_ptr[i+1] - j)*sizeof(*mat->col_val));
        memmove(new_col_ind + edges_processed, mat->col_ind + j, (mat->col_ptr[i+1] - j)*sizeof(*mat->col_ind));
        edges_processed += mat->col_ptr[i+1] - j;
        break;
      }

      while (curr->ind < mat->col_ind[j]) {
        new_col_val[edges_processed] = curr->deg;
        new_col_ind[edges_processed++] = curr->ind;
        prev = curr;
        curr = curr->next;
        free(prev);
        if (curr == NULL) {
          break;
        }
      }

      if (curr != NULL) {
        while (curr->ind == mat->col_ind[j]) {   /* Add degrees to value in original array when col_ind already exists */
          mat->col_val[j] += curr->deg;
          prev = curr;
          curr = curr->next;
          free(prev);
          if (curr == NULL) {
            break;
          }
        }
      }

      if (mat->col_val[j] != 0) {    /* Don't want to add values that have become zero after adding adjacency list */
        new_col_val[edges_processed] = mat->col_val[j];
        new_col_ind[edges_processed++] = mat->col_ind[j];
      }
    }
    /* All values from original col have been added. Walk to end of adjacency list */
    int32_t prev_col_ind = -1;
    while (curr != NULL) {
      if (curr->ind == prev_col_ind) {
        new_col_val[edges_processed-1] += curr->deg;
      }
      else {
        new_col_val[edges_processed] = curr->deg;
        new_col_ind[edges_processed++] = curr->ind;
      }
      prev_col_ind = curr->ind;
      sbm_adj_list_node * prev = curr;
      curr = curr->next;
      free(prev);
    }
  }
  mat->col_ptr[mat->ncols] = edges_processed;
  mat->nnz = edges_processed;

  free(mat->row_ind);
  free(mat->col_ind);
  free(mat->row_val);
  free(mat->col_val);
  mat->row_ind = new_row_ind;
  mat->col_ind = new_col_ind;
  mat->row_val = new_row_val;
  mat->col_val = new_col_val;
}
