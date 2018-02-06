
#ifndef SBM_DYN_CSR_H
#define SBM_DYN_CSR_H

#include <stdint.h>
#include "sbm_dyn_arr.h"

/*
 * Adjacency list node for vertex
 */
typedef struct sbm_adj_list_node {
  /* Index of column or row in array (column or row depending on whether adding to row or column representation of data) */
  int32_t ind;

  /* Degree of edge */
  int32_t deg;

  /* Neighbor in list */
  struct sbm_adj_list_node * next;
} sbm_adj_list_node;

/*
 * Sparse CSR dynamic array
 */
typedef struct {
  /* Number of rows in array */
  int32_t nrows;

  /* Number of columns in array */
  int32_t ncols;

  /* Number of nonzeroes in array */
  int32_t nnz;

  /* Array of pointers to beginnings of rows and columns */
  int32_t * row_ptr;
  int32_t * col_ptr;

  /* Array of indices */
  int32_t * row_ind;
  int32_t * col_ind;

  /* Array of matrix values */
  int32_t * row_val;
  int32_t * col_val;

  /* Adjacency list buffers for each vertex */
  sbm_adj_list_node ** row_adj_lists;
  sbm_adj_list_node ** col_adj_lists;

  /* Number of edges in adjacency lists */
  int32_t adj_list_size;

  /* Number of added nodes in adjacency lists */
  int32_t row_adj_list_new_nodes;
  int32_t col_adj_list_new_nodes;

} sbm_dyn_csr;

/*
 * @brief Allocate space for adjacency list node and initialize values
 *
 * @return Initialized node
 */
sbm_adj_list_node * sbm_adj_list_node_init();

/*
 * @brief Add node to adjacency list
 *
 * @param row_ind Row index of new node
 * @param deg Degree of edge
 * @param prev Previous neighbor
 */
void sbm_adj_list_add_node(
    int32_t ind,
    int32_t deg,
    sbm_adj_list_node * prev);

/*
 * @brief Allocate space for empty CSR
 *
 * @return Dynamic CSR with empty arrays
 */
sbm_dyn_csr * sbm_dyn_csr_alloc();

/*
 * @brief Free dynamic CSR
 *
 * @param mat Structure to free
 */
void sbm_dyn_csr_free();

/*
 * @brief Create a deep copy of a csr
 *
 * @param mat Matrix to copy
 *
 * @return CSR with same contents as original
 */
sbm_dyn_csr * sbm_dyn_csr_deep_cpy(
    sbm_dyn_csr * mat);

/*
 * @brief Add a new value to the linked list buffers for adding to csr later
 *
 * @param mat Matrix to add value to
 * @param row Row to add value to
 * @param col Column to add value to
 * @param val Value to add
 */
void sbm_dyn_csr_add(
    sbm_dyn_csr * mat,
    int32_t row,
    int32_t col,
    int32_t val);

/*
 * @brief Flush adjacency lists to static csr structure
 *
 * @param mat CSR matrix
 */
void sbm_dyn_csr_flush_adj_lists(
    sbm_dyn_csr * mat);

#endif
