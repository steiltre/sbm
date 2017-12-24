
#ifndef SBM_PARTITION_H
#define SBM_PARTITION_H
#include "sbm_dyn_csr.h"

typedef struct {
  double prec;
  double recall;
} sbm_partition_metrics;

/*
 * @brief Compute the pair-wise precision and recall of partitioning
 *
 * @param blk_assn Block assignments predicted by algorithm
 * @param blk_truth True block assignments
 * @param nnodes Number of nodes
 *
 * @return Pairwise precision and recall
 */
sbm_partition_metrics * sbm_pairwise_prec_recall(
    const int32_t * const blk_assn,
    const int32_t * const blk_truth,
    const int32_t nnodes);

/*
 * @brief Initialize block assignments. Only used for initial testing.
 *
 * @param blk_assn Block assignment array
 * @param nnodes Number of nodes
 * @param nblks Number of blocks
 */
void sbm_init_assn_rand(
    int32_t * blk_assn,
    const int32_t nnodes,
    const int32_t nblks);

/*
 * @brief Compute interblock edge counts
 *
 * @param adj_mat Adjacency matrix of graph
 * @param blk_assn Block assignments for nodes
 * @param nblks Number of blocks
 */
sbm_dyn_csr * sbm_compute_iblk_ec(
    sbm_dyn_csr * adj_mat,
    int32_t * blk_assn,
    int32_t nblks);

/*
 * @brief Compute the next number of blocks to use for golden ratio search
 *
 * @param desc_len Description length of current number of blocks
 * @param nblks Current number of blocks
 * @param blk_assn Current block assignments
 * @param iblk_ec Current inter-block edge counts
 * @param best_desc_lens Array for description lengths of search
 * @param best_nblks Number of blocks associated to search configurations
 * @param best_blk_assn Array of block assignments corresponding to each description length
 * @param best_iblk_ec Array of inter-block edge count matrices associated to each description length
 * @param contraction_factor Multiplicative factor for shrinking number of blocks before golden ratio bracket is established
 *
 * @return Next number of blocks to try
 */
int32_t sbm_new_nblks(
    const double desc_len,
    const int32_t nblks,
    int32_t * const blk_assn,
    sbm_dyn_csr * const iblk_ec,
    double * best_desc_lens,
    int32_t * best_nblks,
    int32_t ** best_blk_assn,
    sbm_dyn_csr ** best_iblk_ec,
    const int32_t nnodes,
    const double contraction_factor);

/*
 * @brief Compute the total number of edges incident to a given block
 *
 * @param blk_incident Matrix for holding block incidence counts
 * @param iblk_ec Matrix of interblock edge counts
 * @param nblks Number of blocks
 */
void sbm_compute_blk_incident_edges(
    int32_t * blk_incident,
    sbm_dyn_csr * iblk_ec,
    int32_t nblks);

/*
 * @brief Count number of neighbors from each block for given node
 *
 * @param adj_mat Adjacency matrix of graph
 * @param node ID of node to count for
 * @param blk_assn Block assignments for nodes
 * @param nblks Number of blocks
 * @param node_ngbrs Array to put counts in
 */
void sbm_node_ngbr_blk_count(
    const sbm_dyn_csr * const adj_mat,
    const int32_t node,
    const int32_t * const blk_assn,
    const int32_t nblks,
    int32_t * node_ngbrs);

/*
 * @brief Propose a block to move a node to
 *
 * @param node ID of node to move
 * @param adj_mat Adjacency matrix of graph
 * @param blk_assn Matrix of block assignments
 * @param iblk_ec Matrix of interblock edge counts
 * @param blk_incident Matrix of block incidence counts
 * @param nblks Number of blocks
 * @param force_unique Flag to force new block to be different than original
 *
 * @return ID of block to move node to
 */
int32_t sbm_propose_new_blk(
    const int32_t node,
    const sbm_dyn_csr * const adj_mat,
    const int32_t * const blk_assn,
    const sbm_dyn_csr * const iblk_ec,
    const int32_t * const blk_incident,
    const int32_t nblks,
    const int force_unique);

/*
 * @brief Update edge count matrix from proposed update
 *
 * @param iblk_ec Interbock edge count matrix
 * @param blk_assn Block assignments of nodes
 * @param nblks Number of blocks
 * @param edge_int Array of vertex IDs for terminal end of edges on node moving blocks
 * @param edge_val Array of edge weights for edges attached to moving node
 * @param nedges Number of edges on moving node
 * @param from_blk ID of block node is moving from
 * @param to_blk ID of block node is moving to
 */
void sbm_update_iblk_ec(
    sbm_dyn_csr * iblk_ec,
    const int32_t * const blk_assn,
    const int32_t nblks,
    const int32_t * const edge_ind,
    const int32_t * const edge_val,
    const int32_t nedges,
    const int32_t from_blk,
    const int32_t to_blk);

/*
 * @brief Update count of edges incident to blocks after proposed move
 *
 * @param blk_incident Array of edge counts incident to each block
 * @param blk_assn Block assignments of nodes
 * @param nblks Number of blocks
 * @param edge_val Array of edge weights for edges attached to moving node
 * @param nedges Number of edges on moving node
 * @param from_blk ID of block node is moving from
 * @param to_blk ID of block node is moving to
 */
void sbm_update_blk_incident(
    int32_t * blk_incident,
    const int32_t nblks,
    const int32_t * const edge_val,
    const int32_t nedges,
    const int32_t from_blk,
    const int32_t to_blk);

/*
 * @brief Perform block merges
 *
 * @param blk_moves Best merge for each block
 * @param blk_del_ent Change in entropy associated to best move for each block
 * @param blks_to_merge Number of blocks to merge
 * @param nblks Total number of blocks
 * @param blk_map Array to hold new label for each block
 *
 * @return New number of blocks
 */
int32_t sbm_merge_blks(
    const int32_t * const blk_moves,
    const double * const blk_del_ent,
    const int32_t blks_to_merge,
    const int32_t nblks,
    int32_t * blk_map);

/*
 * @brief Relabel block assignments
 *
 * @param blk_assn Block assignments for nodes
 * @param blk_map Map giving new labels
 * @param nnodes Number of nodes
 */
void sbm_relabel_blks(
    int32_t * blk_assn,
    const int32_t * const blk_map,
    const int32_t nnodes);

/*
 * @brief Compute the entropy of the current partitioning
 *
 * @param iblk_ec Matrix of interblock edge counts
 * @param blk_incident Count of edges incident to blocks
 * @param nblks Number of blocks
 * @param nnodes Number of nodes
 *
 * @return Entropy of partitioning
 */
double sbm_desc_len(
    const sbm_dyn_csr * const iblk_ec,
    const int32_t * const blk_incident,
    const int32_t nblks,
    const int32_t nnodes);

/*
 * @brief Compute the change in entropy from proposed move
 *
 * @param iblk_ec Interblock edge count matrix
 * @param blk_incident Count of edges incident to blocks
 * @param node_ngbrs Array with number of neighbors of current node in each block
 * @param self_edges Number of edges node has with self
 * @param nblks Number of blocks
 * @param from_blk Block node is moving out of
 * @param to blk Block node is moving into
 *
 * @return Change in entropy
 */
double sbm_delta_entropy(
    const sbm_dyn_csr * const iblk_ec_new,
    const int32_t * const blk_incident,
    const int32_t * const node_ngbrs,
    const int32_t self_edges,
    const int32_t nblks,
    const int32_t from_blk,
    const int32_t to_blk);

/*
 * @brief Compute the Metropolis-Hastings correction
 *
 * @param iblk_ec Interblock edge-count matrix
 * @param blk_incident Array of edge counts incident to each block
 * @param blk_assn Block assignments of nodes
 * @param node_ngbrs Array with number of neighbors of current node in each block
 * @param nblks Number of blocks
 * @param from_blk ID of block node is moving from
 * @param to_blk ID of block node is moving to
 *
 * @return Metropolis-Hastings correction
 */
double sbm_met_has(
    const sbm_dyn_csr * const iblk_ec,
    const int32_t * const blk_incident,
    const int32_t * const blk_assn,
    const int32_t * const node_ngbrs,
    const int32_t nblks,
    const int32_t from_blk,
    const int32_t to_blk);

/*
 * @brief Compute partitioning
 *
 * @param adj_mat Adjacency matrix of graph
 * @param beta Parameter for exponential weight of entropy change in Metropolis-Hastings correction
 *
 * @return Array of block labels
 */
int32_t * sbm_partition_nblks(
    sbm_dyn_csr * adj_mat,
    const double beta);

#endif
