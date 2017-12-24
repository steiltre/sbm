
#ifndef SBM_FILE_H
#define SBM_FILE_H

#include "sbm_dyn_csr.h"

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
    int32_t indexing);

/*
 * @brief Read true partition labels from .tsv file
 *
 * @param fname File to read
 * @param nnodes Number of nodes
 *
 * @return Array of true labels
 */
int32_t * sbm_read_truth_tsv(
    char const * const fname,
    int32_t nnodes);

/*
 * @brief Output results of partitioning
 *
 * @param fname File to write to
 * @param blk_assn Vector of block assignments
 * @param nnodes Number of nodes
 */
void sbm_write_partitions(
    char const * const fname,
    int32_t * blk_assn,
    int32_t nnodes);

#endif
