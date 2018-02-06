
#ifndef SBM_FILE_H
#define SBM_FILE_H

#include "sbm_dyn_csr.h"

/*
 * @brief Read .tsv file and append edges to existing graph
 *
 * @param fname File to read
 * @param adj_mat Adjacency matrix to append edges to
 * @param indexing Lowest index of nodes in file
 */
void sbm_read_tsv_append(
		char const * const fname,
		sbm_dyn_csr * adj_mat,
		const int32_t indexing);

/*
 * @brief Read .tsv file of graph
 *
 * @param fname File to read
 * @param indexing Lowest index of nodes in file
 *
 * @return CSR array of graph
 */
sbm_dyn_csr * sbm_read_tsv(
    char const * const fname,
    const int32_t indexing);

/*
 * @brief Read .tsv file to find the number of nodes
 *
 * @param fname Name of file
 * @param indexing Lowest index of nodes in file
 *
 * @return Number of nodes in graph
 */
int32_t sbm_count_nodes_tsv(
		char const * const fname,
		const int32_t nnodes);

/*
 * @brief Read true partition labels from .tsv file
 *
 * @param fname File to read
 * @param nnodes Number of nodes
 * @param indexing Lowest index of nodes in file
 *
 * @return Array of true labels
 */
int32_t * sbm_read_truth_tsv(
    char const * const fname,
    const int32_t nnodes,
		const int32_t indexing);

/*
 * @brief Output results of partitioning
 *
 * @param fname File to write to
 * @param blk_assn Vector of block assignments
 * @param nnodes Number of nodes
 */
void sbm_write_partitions(
    char const * const fname,
    const int32_t * const blk_assn,
    const int32_t nnodes);

#endif
