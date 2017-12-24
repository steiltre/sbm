
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sbm_dyn_csr.h"
#include "sbm_util.h"

sbm_adj_list_node * sbm_adj_list_node_init()
{
  sbm_adj_list_node * node = malloc(sizeof(*node));

  node->row_ind = -1;
  node->deg = -1;
  node->next = NULL;

  return node;
}

void sbm_adj_list_add_node(
    int32_t row_ind,
    int32_t deg,
    sbm_adj_list_node * prev)
{
  sbm_adj_list_node * node = sbm_adj_list_node_init();

  node->row_ind = row_ind;
  node->deg = deg;
  node->next = prev->next;
  prev->next = node;
}

sbm_dyn_csr * sbm_dyn_csr_alloc()
{
  sbm_dyn_csr * mat = malloc( sizeof(*mat) );

  mat->row_ptr = malloc( sizeof(*mat->row_ptr) );
  mat->row_ind = NULL;
  mat->val = NULL;

  mat->adj_lists = malloc(sizeof(*mat->adj_lists));

  mat->nrows = 0;
  mat->nnz = 0;
  mat->adj_list_size = 0;
  mat->adj_list_new_nodes = 0;

  mat->row_ptr[0] = 0;

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
  if (mat->val != NULL) {
    free(mat->val);
  }

  for (int32_t i=0; i<mat->nrows; i++) {
    sbm_adj_list_node * prev = mat->adj_lists[i];
    sbm_adj_list_node * curr = mat->adj_lists[i]->next;
    while (curr != NULL) {
      free(prev);
      prev = curr;
    }
    free(prev);
  }
  free(mat->adj_lists);

  free(mat);
}

/* Assumes adjacency lists have been flushed */
sbm_dyn_csr * sbm_dyn_csr_deep_cpy(
    sbm_dyn_csr * mat)
{
  sbm_dyn_csr * ret = malloc(sizeof(*ret));;

  ret->row_ptr = malloc((mat->nrows+1)*sizeof(*ret->row_ptr));
  memcpy(ret->row_ptr, mat->row_ptr, (mat->nrows+1)*sizeof(*ret->row_ptr));

  ret->row_ind = malloc(mat->nnz*sizeof(*ret->row_ind));
  memcpy(ret->row_ind, mat->row_ind, mat->nnz*sizeof(*ret->row_ind));

  ret->val = malloc(mat->nnz*sizeof(*ret->val));
  memcpy(ret->val, mat->val, mat->nnz*sizeof(*ret->val));

  ret->adj_lists = malloc(mat->nrows*sizeof(*ret->adj_lists));
  for (int32_t i=0; i<mat->nrows; i++) {
    ret->adj_lists[i] = sbm_adj_list_node_init();
  }

  ret->nrows = mat->nrows;
  ret->nnz = mat->nnz;
  ret->adj_list_size = mat->adj_list_size;
  ret->adj_list_new_nodes = mat->adj_list_new_nodes;

  return ret;
}

void sbm_dyn_csr_add(
    sbm_dyn_csr * mat,
    int32_t row,
    int32_t col,
    int32_t val)
{
  if (row >= mat->nrows + mat->adj_list_new_nodes) {
    mat->adj_lists = (sbm_adj_list_node **) realloc(mat->adj_lists, (row+1) * sizeof(*mat->adj_lists));
    for (int32_t i = mat->nrows + mat->adj_list_new_nodes; i<row+1; i++) {
      mat->adj_lists[i] = sbm_adj_list_node_init();
    }
    mat->adj_list_new_nodes = row - mat->nrows + 1;
  }
  sbm_adj_list_node * curr = mat->adj_lists[row];

  if (curr->next != NULL) {
    while (curr->next->row_ind <= col) {    /* Sort adjacency list in ascending order for ease of flushing */
      curr = curr->next;

      if (curr->next == NULL) {
        break;
      }
    }
  }

  /*
  if (curr == NULL) {
    sbm_ll_add_node(col, val, curr);
    mat->adj_list_size++;
  }
  else { */
    if (curr->row_ind == col) {     /* Same entry in matrix */
      curr->deg += val;
    }
    else {
      sbm_adj_list_add_node(col, val, curr);
      mat->adj_list_size++;
    }
  //}
}

void sbm_dyn_csr_flush_adj_lists(
    sbm_dyn_csr * mat)
{
  /* Allocate new csr arrays to fill */
  int32_t * new_row_ind = malloc( (mat->nnz + mat->adj_list_size) * sizeof(*new_row_ind) );
  int32_t * new_val = malloc( (mat->nnz + mat->adj_list_size) * sizeof(*new_val) );

  /*
  mat->nnz += mat->adj_list_size;
  mat->adj_list_size = 0;
  mat->val = (int32_t *) realloc(mat->val, mat->nnz * sizeof(*mat->val));
  mat->row_ind = (int32_t *) realloc(mat->row_ind, mat->nnz * sizeof(*mat->row_ind));
  */
  mat->row_ptr = (int32_t *) realloc(mat->row_ptr, (mat->nrows + mat->adj_list_new_nodes + 1) * sizeof(*mat->row_ptr));
  for (int32_t i=0; i<mat->adj_list_new_nodes; i++) {
    mat->row_ptr[ mat->nrows + i + 1] = mat->row_ptr[mat->nrows];
  }
  mat->nrows += mat->adj_list_new_nodes;
  mat->adj_list_new_nodes = 0;
  mat->adj_list_size = 0;

  int32_t edges_processed = 0;
  for (int32_t i=0; i<mat->nrows; i++) {
    sbm_adj_list_node * prev;
    sbm_adj_list_node * curr = mat->adj_lists[i]->next;
    mat->adj_lists[i]->next = NULL;

    int32_t prev_row_ptr = mat->row_ptr[i];
    mat->row_ptr[i] = edges_processed;

    for (int32_t j=prev_row_ptr; j<mat->row_ptr[i+1]; j++) {     /* Walk along row of old csr */
      if (curr == NULL) {    /* No more entries in adjacency list. Copy rest of original row */
        memmove(new_val + edges_processed, mat->val + j, (mat->row_ptr[i+1] - j)*sizeof(*mat->val));
        memmove(new_row_ind + edges_processed, mat->row_ind + j, (mat->row_ptr[i+1] - j)*sizeof(*mat->row_ind));
        edges_processed += mat->row_ptr[i+1] - j;
        break;
      }

      while (curr->row_ind < mat->row_ind[j]) {
        new_val[edges_processed] = curr->deg;
        new_row_ind[edges_processed++] = curr->row_ind;
        prev = curr;
        curr = curr->next;
        free(prev);
        if (curr == NULL) {
          break;
        }
      }

      if (curr != NULL) {
        while (curr->row_ind == mat->row_ind[j]) {   /* Add degrees to value in original array when row_ind already exists */
          mat->val[j] += curr->deg;
          prev = curr;
          curr = curr->next;
          free(prev);
          if (curr == NULL) {
            break;
          }
        }
      }

      if (mat->val[j] != 0) {    /* Don't want to add values that have become zero after adding adjacency list */
        new_val[edges_processed] = mat->val[j];
        new_row_ind[edges_processed++] = mat->row_ind[j];
      }
    }
    /* All values from original row have been added. Walk to end of adjacency list */
    int32_t prev_row_ind = -1;
    while (curr != NULL) {
      if (curr->row_ind == prev_row_ind) {
        new_val[edges_processed-1] += curr->deg;
      }
      else {
        new_val[edges_processed] = curr->deg;
        new_row_ind[edges_processed++] = curr->row_ind;
      }
      prev_row_ind = curr->row_ind;
      sbm_adj_list_node * prev = curr;
      curr = curr->next;
      free(prev);
    }
  }
  mat->row_ptr[mat->nrows] = edges_processed;
  mat->nnz = edges_processed;

  free(mat->row_ind);
  free(mat->val);
  mat->row_ind = new_row_ind;
  mat->val = new_val;
}
