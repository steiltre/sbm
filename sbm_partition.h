
#ifndef SBM_PARTITION_H
#define SBM_PARTITION_H
#include "sbm_dyn_csr.h"

typedef struct {
  double prec;
  double recall;
} sbm_partition_metrics;

/*
 * @brief Structure for holding preferred block merge for sorting
 */
typedef struct {
  double del_ent;
  int32_t blk_lbl;
  int32_t blk_move;
} sbm_blk_merge_info;

/*
 * @brief Structure for holding block labels and sizes
 */
typedef struct {
	int32_t blk_lbl;
	int32_t blk_size;
} sbm_blk_size;

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
 * @brief Comparison operator used for sorting sbm_blk_merge_info objects
 *        based on entropy change
 *
 * @param a Pointer to first object
 * @param b Pointer to second object
 *
 * @return Integer specifying which element is larger
 */
int sbm_gt_sbm_blk_merge_info(
    const void * a_ptr,
    const void * b_ptr);

/*
 * @brief Comparison operator used for sorting sbm_blk_size objects
 *
 * @param a Pointer to first object
 * @param b Pointer to second object
 *
 * @return Integer specifying which element is larger
 */
int sbm_lt_sbm_blk_size(
		const void * a_ptr,
		const void * b_ptr);

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
 * @brief Randomly bisect a block
 *
 * @param blk_lbl Label of block to bisect
 * @param blk_assn Current block assignments of nodes
 * @param blk_sizes Sizes of each block
 * @param nblks Current number of blocks
 */
void sbm_random_bisect(
		const int32_t blk_lbl,
		int32_t * blk_assn,
    int32_t * blk_sizes,
		const int32_t nblks,
		const int32_t nnodes);

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
 * @brief Compute the degree of each block
 *
 * @param blk_deg_in Matrix for holding block in-degrees
 * @param blk_deg_out Matrix for holding block out-degrees
 * @param iblk_ec Matrix of interblock edge counts
 * @param nblks Number of blocks
 */
void sbm_compute_blk_deg(
    int32_t * blk_in_deg,
    int32_t * blk_out_deg,
    sbm_dyn_csr * iblk_ec,
    int32_t nblks);

/*
 * @brief Count number of in- and out-neighbors from each block for given node
 *
 * @param adj_mat Adjacency matrix of graph
 * @param node ID of node to count for
 * @param blk_assn Block assignments for nodes
 * @param nblks Number of blocks
 * @param node_ngbrs_in Array to hold in-neighbors
 * @param node_ngbrs_out Array to hold out-neighbors
 */
void sbm_node_ngbr_blk_count(
    const sbm_dyn_csr * const adj_mat,
    const int32_t node,
    const int32_t * const blk_assn,
    const int32_t nblks,
    int32_t * node_ngbrs_in,
    int32_t * node_ngbrs_out);

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
    const int32_t * const blk_in_deg,
    const int32_t * const blk_out_deg,
    const int32_t nblks,
    const int force_unique);

/*
 * @brief Update edge count matrix from proposed update
 *
 * @param adj_mat Adjacency matrix
 * @param iblk_ec Interbock edge count matrix
 * @param blk_assn Block assignments of nodes
 * @param nblks Number of blocks
 * @param node_id Index of node being moved
 * @param from_blk ID of block node is moving from
 * @param to_blk ID of block node is moving to
 */
void sbm_update_iblk_ec(
    const sbm_dyn_csr * const adj_mat,
    sbm_dyn_csr * iblk_ec,
    const int32_t * const blk_assn,
    const int32_t nblks,
    const int32_t node_id,
    const int32_t from_blk,
    const int32_t to_blk);

/*
 * @brief Update count of edges incident to blocks after proposed move
 *
 * @param blk_in_deg Array of incoming degree for each block
 * @param blk_out_deg Array of outgoing degree for each block
 * @param adj_mat Adjacency matrix
 * @param blk_assn Block assignments of nodes
 * @param nblks Number of blocks
 * @param node_id Index of node that is moving
 * @param from_blk ID of block node is moving from
 * @param to_blk ID of block node is moving to
 */
void sbm_update_blk_deg(
    int32_t * blk_in_deg,
    int32_t * blk_out_deg,
    const sbm_dyn_csr * const adj_mat,
    const int32_t nblks,
    const int32_t node_id,
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
 * @brief Recalculate block sizes after block merges
 *
 * @param new_blk_sizes Array for holding new block sizes
 * @param old_blk_sizes Array of block sizes before merging
 * @param blk_map Map giving new labels
 * @param old_nblks Number of blocks before merging
 */
void sbm_merged_blk_sizes(
    int32_t * new_blk_sizes,
    const int32_t * const old_blk_sizes,
    const int32_t * const blk_map,
    const int32_t old_nblks);

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
    const int32_t * const blk_in_deg,
    const int32_t * const blk_out_deg,
    const int32_t nblks,
    const int32_t nnodes);

/*
 * @brief Compute the change in entropy from proposed move
 *
 * @param adj_mat Adjacency matrix
 * @param iblk_ec Interblock edge count matrix
 * @param blk_in_deg Total degree of edges going to each block
 * @param blk_out_deg Total degree of edges going from each block
 * @param blk_assn Block assignments for nodes
 * @param node_id Index of node being moved
 * @param nblks Number of blocks
 * @param from_blk Block node is moving out of
 * @param to blk Block node is moving into
 *
 * @return Change in entropy
 */
double sbm_delta_entropy(
    const sbm_dyn_csr * const adj_mat,
    const sbm_dyn_csr * const iblk_ec,
    const int32_t * const blk_in_deg,
    const int32_t * const blk_out_deg,
    const int32_t * const blk_assn,
    const int32_t node_id,
    const int32_t nblks,
    const int32_t from_blk,
    const int32_t to_blk);

/*
 * @brief Compute the Metropolis-Hastings correction
 *
 * @param iblk_ec Interblock edge-count matrix
 * @param blk_in_deg Array of total in-degree for each block
 * @param blk_out_deg Array of total out-degree for each block
 * @param blk_assn Block assignments of nodes
 * @param in_ngbrs Array with number of incoming edges from each block to given node
 * @param out_ngbrs Array with number of outgoing edges to each block from given node
 * @param nblks Number of blocks
 * @param from_blk ID of block node is moving from
 * @param to_blk ID of block node is moving to
 *
 * @return Metropolis-Hastings correction
 */
double sbm_met_has(
    const sbm_dyn_csr * const iblk_ec,
    const int32_t * const blk_in_deg,
    const int32_t * const blk_out_deg,
    const int32_t * const blk_assn,
    const int32_t * const in_ngbrs,
    const int32_t * const out_ngbrs,
    const int32_t nblks,
    const int32_t from_blk,
    const int32_t to_blk);

/*
 * @brief Compute partitions for each time-step of a dynamic graph
 *
 * @param ifile_prefix Beginning of filenames for graph data
 * @param nfiles Number of data files
 * @param beta Parameter for exponential weight of entropy change in Metropolis-Hastings correction
 */
void sbm_stream_partition(
		char * ifile_prefix,
		int nfiles,
		double beta);

/*
 * @brief Compute partitioning
 *
 * @param adj_mat Adjacency matrix of graph
 * @param beta Parameter for exponential weight of entropy change in Metropolis-Hastings correction
 * @param blk_assn Initial guess at partitioning (will be updated)
 * @param nblks Number of blocks in initial guess
 *
 * @return Number of blocks in optimal partitioning
 */
int32_t sbm_partition(
    sbm_dyn_csr * adj_mat,
    const double beta,
		int32_t * blk_assn,
		int32_t nblks);

#endif
