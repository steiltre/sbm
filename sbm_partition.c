
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "sbm_partition.h"
#include "sbm_stack.h"
#include "sbm_util.h"

#define MAX_TRIALS 10
#define CONTRACTION_FACTOR 0.5

double block_merge_time = 0;
double nodal_update_time = 0;
double flushing_time = 0;

sbm_partition_metrics * sbm_pairwise_prec_recall(
    const int32_t * const blk_assn,
    const int32_t * const blk_truth,
    const int32_t nnodes)
{
  int tmpm = 0;   /* Truth Match Prediction Match */
  int tmpd = 0;   /* Truth Match Prediction Differ */
  int tdpm = 0;   /* Truth Differ Prediction Match */
  sbm_partition_metrics * res = malloc(sizeof(*res));

  for (int32_t i=0; i<nnodes-1; i++) {
    for (int32_t j=i+1; j<nnodes; j++) {
      if (blk_truth[i] == blk_truth[j] && blk_assn[i] == blk_assn[j]) {
        tmpm++;
      }
      if (blk_truth[i] == blk_truth[j] && blk_assn[i] != blk_assn[j]) {
        tmpd++;
      }
      if (blk_truth[i] != blk_truth[j] && blk_assn[i] == blk_assn[j]) {
        tdpm++;
      }
    }
  }
  res->prec = ( (double) tmpm ) / (tmpm + tdpm);
  res->recall = ( (double) tmpm ) / (tmpm + tmpd);

  return res;
}

void sbm_init_assn_rand(
    int32_t * blk_assn,
    const int32_t nnodes,
    const int32_t nblks)
{
  for (int32_t i=0; i<nnodes; i++) {
    blk_assn[i] = rand() % nblks;
  }
}

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
    const double contraction_factor)
{
  //if (desc_len < best_desc_lens[1] || best_desc_lens[1] == -1) {    /* New desc len is better than previous best */
  if ( (desc_len - best_desc_lens[1])/best_desc_lens[1] < .01 || best_desc_lens[1] == -1) {
    int slide_pos;
    if (nblks > best_nblks[1]) {    /* New desc len corresponds to a larger number of blocks */
      slide_pos = 0;
    }
    else {
      slide_pos = 2;
    }
    best_nblks[slide_pos] = best_nblks[1];
    best_nblks[1] = nblks;

    best_desc_lens[slide_pos] = best_desc_lens[1];
    best_desc_lens[1] = desc_len;

    int32_t * temp_int = best_blk_assn[slide_pos];
    best_blk_assn[slide_pos] = best_blk_assn[1];
    best_blk_assn[1] = temp_int;
    if (best_blk_assn[1] == NULL) {
      best_blk_assn[1] = malloc(nnodes*sizeof(*best_blk_assn[1]));
    }
    memcpy(best_blk_assn[1], blk_assn, nnodes * sizeof(*best_blk_assn[1]));

    if (best_iblk_ec[slide_pos] != NULL) {
      sbm_dyn_csr_free(best_iblk_ec[slide_pos]);
    }
    best_iblk_ec[slide_pos] = best_iblk_ec[1];
    best_iblk_ec[1] = sbm_dyn_csr_deep_cpy(iblk_ec);
  }
  else {
    int new_pos;
    if (nblks > best_nblks[1]) {    /* More blocks than optimal */
      new_pos = 2;
    }
    else {
      new_pos = 0;
    }
    best_nblks[new_pos] = nblks;
    best_desc_lens[new_pos] = desc_len;

    if (best_blk_assn[new_pos] == NULL) {
      best_blk_assn[new_pos] = malloc(nnodes*sizeof(*best_blk_assn[new_pos]));
    }
    memcpy(best_blk_assn[new_pos], blk_assn, nnodes*sizeof(*blk_assn));

    if (best_iblk_ec[new_pos] != NULL) {
      sbm_dyn_csr_free(best_iblk_ec[new_pos]);
    }
    best_iblk_ec[new_pos] = sbm_dyn_csr_deep_cpy(iblk_ec);
  }

  /* Determine next number of blocks to use */
  if (best_desc_lens[0] == -1) {    /* Still don't have bracket established */
    return ( (int) nblks * contraction_factor);
  }
  else if (best_desc_lens[2] == -1) {    /* Still don't have bracket established (in other direction) */
    return ( (int) nblks / contraction_factor);
  }
  else {
    if (best_nblks[2] - best_nblks[0] == 2) {    /* Search has concluded */
      return 0;
    }
    else if (best_nblks[2] - best_nblks[1] < best_nblks[1] - best_nblks[0]) {
      return (int) (best_nblks[1] - ( 0.618 * (best_nblks[1] - best_nblks[0]) ) );
    }
    else {
      return (int) (best_nblks[1] + ( 0.618 * (best_nblks[2] - best_nblks[1]) ) );
    }
  }
}

/* Storing edge counts as dense matrix. Will need to switch to sparse for when each node is its own block (O(n^2) memory). */
sbm_dyn_csr *  sbm_compute_iblk_ec(
    sbm_dyn_csr * adj_mat,
    int32_t * blk_assn,
    int32_t nblks)
{
  /*
  for (int32_t i=0; i<nblks*nblks; i++) {
    iblk_ec[i] = 0;
  }
  */
  sbm_dyn_csr * iblk_ec = sbm_dyn_csr_alloc();

  for (int32_t i=0; i<adj_mat->nrows; i++) {
    for (int32_t j=adj_mat->row_ptr[i]; j<adj_mat->row_ptr[i+1]; j++) {
      //iblk_ec[ blk_assn[i]*nblks + blk_assn[adj_mat->row_ind[j]-1] ]++;          /* Need to subtract 1 from adj_mat->row_ind[j] when nodes are 1-indexed */
      sbm_dyn_csr_add( iblk_ec, blk_assn[i], blk_assn[adj_mat->row_ind[j]], adj_mat->val[j] );
    }
  }
  double start = monotonic_seconds();
  sbm_dyn_csr_flush_adj_lists(iblk_ec);
  flushing_time += monotonic_seconds() - start;
  return iblk_ec;
}

void sbm_compute_blk_incident_edges(
    int32_t * blk_incident,
    sbm_dyn_csr * iblk_ec,
    int32_t nblks)
{
  for (int32_t i=0; i<nblks; i++) {
    blk_incident[i] = 0;
  }

  /* Sum up rows of iblk_ic to give blk_incident */
  for (int32_t i=0; i<nblks; i++) {
    for (int32_t j=iblk_ec->row_ptr[i]; j<iblk_ec->row_ptr[i+1]; j++) {
      blk_incident[i] += iblk_ec->val[j];
    }
  }
}

void sbm_node_ngbr_blk_count(
    const sbm_dyn_csr * const adj_mat,
    const int32_t node,
    const int32_t * const blk_assn,
    const int32_t nblks,
    int32_t * node_ngbrs)
{
  for (int32_t i=adj_mat->row_ptr[node]; i<adj_mat->row_ptr[node+1]; i++) {
    node_ngbrs[ blk_assn[adj_mat->row_ind[i]] ] += adj_mat->val[i];
    node_ngbrs[nblks] += adj_mat->val[i];
  }
}

/* Method of proposing block can be found in "Efficient Monte Carlo and greedy heuristic for the inference of stochastic block models" by Peixoto.
 * We have fixed \epsilon = 1 in this implementation */
int32_t sbm_propose_new_blk(
    const int32_t node,
    const sbm_dyn_csr * const adj_mat,
    const int32_t * const blk_assn,
    const sbm_dyn_csr * const iblk_ec,
    const int32_t * const blk_incident,
    const int32_t nblks,
    const int force_unique)
{
  /* Select a random neighbor of node and get its block */
  int32_t nngbrs = adj_mat->row_ptr[node+1] - adj_mat->row_ptr[node];          /* Subtract 1 from node when nodes are 1-indexed */

  if (nngbrs == 0) {    /* Return random proposal if node has no neighbors */
    int32_t prop;
    if (force_unique) {
      prop = rand() % (nblks-1);
      if (prop >= blk_assn[node]) {  /* Don't allow proposal to be same as original block */
        prop++;
      }
    }
    else {
      prop = rand() % nblks;
    }
    return prop;
  }

  int32_t ngbr = adj_mat->row_ind[ adj_mat->row_ptr[node] + rand()%nngbrs ];       /* Subtract 1 when nodes are 1-indexed */
  int32_t blk_ngbr = blk_assn[ngbr];       /* Subtract 1 when nodes are 1-indexed */
  int32_t ngbr_blk_nngbrs = blk_incident[blk_ngbr];
  if (force_unique) {     /* Make sure old block is not considered */
    for (int32_t idx = iblk_ec->row_ptr[blk_ngbr]; idx < iblk_ec->row_ptr[blk_ngbr+1]; idx++) {
      if (iblk_ec->row_ind[idx] == blk_assn[node]) {
        ngbr_blk_nngbrs -= iblk_ec->val[idx];
      }
    }
  }


  /* Accept random block? */
  double accept_rand_blk_prob = ( (double) nblks) / (blk_incident[blk_ngbr] + nblks);
  double r = sbm_rand_dbl(0, 1);

  if (r < accept_rand_blk_prob || ngbr_blk_nngbrs == 0) {
    int32_t prop = -1;
    if (force_unique) {
      prop = rand() % (nblks-1);
      if (prop >= blk_assn[node]) {  /* Don't allow proposal to be same as original block */
        prop++;
      }
    }
    else {
      prop = rand() % nblks;
    }
    return prop;
  }
  else {
    int32_t counter;
    //for (int j=0; j<20; j++) {
    counter = rand() % ngbr_blk_nngbrs;
    int32_t prop = -1;
    //printf("%d, %d\n", counter, blk_ngbrs);
    for (int32_t i=iblk_ec->row_ptr[blk_ngbr]; i<iblk_ec->row_ptr[blk_ngbr+1]; i++) {
      if (iblk_ec->row_ind[i] != blk_assn[node] || !force_unique) {
        counter -= iblk_ec->val[i];
      }
      if (counter < 0) {
        prop = iblk_ec->row_ind[i];
        break;
        //printf("%d\n", iblk_ec->row_ind[i]);
        //return iblk_ec->row_ind[i];
      }
    //}
    }
    //exit(EXIT_FAILURE);
    //printf("Shouldn't get here\nCounter: %d\n", counter);
    return prop;
  }
}

void sbm_update_iblk_ec(
    sbm_dyn_csr * iblk_ec,
    const int32_t * const blk_assn,
    const int32_t nblks,
    const int32_t * const edge_ind,
    const int32_t * const edge_val,
    const int32_t nedges,
    const int32_t from_blk,
    const int32_t to_blk)
{
  for (int32_t i=0; i<nedges; i++) {
    int32_t ngbr_blk = blk_assn[edge_ind[i]];
    sbm_dyn_csr_add( iblk_ec, ngbr_blk, from_blk, -1 * edge_val[i] );
    sbm_dyn_csr_add( iblk_ec, ngbr_blk, to_blk, edge_val[i] );

    sbm_dyn_csr_add( iblk_ec, from_blk, ngbr_blk, -1 * edge_val[i] );
    sbm_dyn_csr_add( iblk_ec, to_blk, ngbr_blk, edge_val[i] );
  }
  double start = monotonic_seconds();
  sbm_dyn_csr_flush_adj_lists(iblk_ec);
  flushing_time += monotonic_seconds() - start;
}

void sbm_update_blk_incident(
    int32_t * blk_incident,
    const int32_t nblks,
    const int32_t * const edge_val,
    const int32_t nedges,
    const int32_t from_blk,
    const int32_t to_blk)
{
  for (int i=0; i<nedges; i++) {
    blk_incident[from_blk] -= edge_val[i];
    blk_incident[to_blk] += edge_val[i];
  }
}

int32_t sbm_merge_blks(
    const int32_t * const blk_moves,
    const double * const blk_del_ent,
    const int32_t blks_to_merge,
    const int32_t nblks,
    int32_t * blk_map)
{
  double * blk_del_ent_sort = malloc(nblks * sizeof(*blk_del_ent_sort));
  memcpy(blk_del_ent_sort, blk_del_ent, nblks * sizeof(*blk_del_ent));
  qsort(blk_del_ent_sort, nblks, sizeof(*blk_del_ent), sbm_gt_dbl);

  double ent_threshold = blk_del_ent_sort[blks_to_merge];

  free(blk_del_ent_sort);

  /* Create array for map between old and new labels */
  for (int32_t i=0; i<nblks; i++) {
    blk_map[i] = i;
  }

  /* Create stacks for moving already merged blocks */
  sbm_ll_node ** blk_stacks = malloc(nblks*sizeof(*blk_stacks));
  for (int32_t i=0; i<nblks; i++) {
    *(blk_stacks + i) = sbm_ll_node_alloc();
  }

  /* Create array for tracking which block labels are no longer used for relabeling at end */
  int * merged_lbls = calloc(nblks, sizeof(*merged_lbls));

  int32_t removed_blks = 0;  /* Counts number of actual merges */
  for (int32_t i=0; i<nblks; i++) {
    if (blk_del_ent[i] <= ent_threshold) {  /* Merge blocks */
      blk_map[i] = blk_map[blk_moves[i]];
      while (blk_stacks[i]->next != NULL) {    /* Drag along blocks that have already merged */
        int32_t already_moved_blk = sbm_ll_node_pop(blk_stacks[i]);
        blk_map[already_moved_blk] = blk_map[i];
      }
      sbm_ll_node_push(i, *(blk_stacks + blk_moves[i]));

      /* Prevent pairs of blocks that both want to merge with each other from being double counted */
      if (i>blk_moves[i] || i != blk_moves[blk_moves[i]]) {
        merged_lbls[i] = 1;    /* Mark label as being unused */
        removed_blks++;
      }
    }
  }

  /* Perform a scan on merged labels */
  for (int32_t i=1; i<nblks; i++) {
    merged_lbls[i] += merged_lbls[i-1];
  }

  /* Compress block labels toward zero */
  for (int32_t i=0; i<nblks; i++) {
    blk_map[i] -= merged_lbls[blk_map[i]];
  }

  for (int32_t i=0; i<nblks; i++) {
    sbm_ll_node * curr = blk_stacks[i]->next;
    sbm_ll_node * prev = blk_stacks[i];
    while (curr != NULL) {
      free(prev);
      prev = curr;
      curr = curr->next;
    }
    free(prev);
  }
  free(blk_stacks);

  free(merged_lbls);

  return nblks - removed_blks;
}

void sbm_relabel_blks(
    int32_t * blk_assn,
    const int32_t * const blk_map,
    const int32_t nnodes)
{
  for (int32_t i=0; i<nnodes; i++) {
    blk_assn[i] = blk_map[blk_assn[i]];
  }
}

double sbm_desc_len(
    const sbm_dyn_csr * const iblk_ec,
    const int32_t * const blk_incident,
    const int32_t nblks,
    const int32_t nnodes)
{
  int32_t nedges = 0;
  for (int i=0; i<nblks; i++) {
    nedges += blk_incident[i];
  }
  nedges = nedges / 2;

  double sum = 0;
  double ratio = ( (double) nblks * nblks ) / nedges;

  sum = nedges * ( (1 + ratio)*log(1+ratio) - ratio*log(ratio) ) + nnodes * log(nblks);

  for (int i=0; i<nblks; i++) {
    for (int j=iblk_ec->row_ptr[i]; j<iblk_ec->row_ptr[i+1]; j++) {
      double log1 = log( (double) iblk_ec->val[j] );
      double log2 = log( (double) blk_incident[i] * blk_incident[ iblk_ec->row_ind[j] ] );
      //double factor = ( (double) iblk_ec[i*nblks + j]) / (blk_incident[i]*blk_incident[j]);
      sum -= iblk_ec->val[j] * ( log1 - log2 );
    }
  }
  return sum;
}

/* Assumes graph is undirected. Need to store iblk_ec by row and column to efficiently use directed edges. */
double sbm_delta_entropy(
    const sbm_dyn_csr * const iblk_ec,
    const int32_t * const blk_incident,
    const int32_t * const node_ngbrs,
    const int32_t self_edges,
    const int32_t nblks,
    const int32_t from_blk,
    const int32_t to_blk)
{
  double init_ent = 0;   /* Contribution to initial entropy from affected rows and columns */

  for (int32_t i=iblk_ec->row_ptr[from_blk]; i<iblk_ec->row_ptr[from_blk+1]; i++) {
    double ent = iblk_ec->val[i] * log(( (double) iblk_ec->val[i] ) / (blk_incident[from_blk] * blk_incident[ iblk_ec->row_ind[i] ]));
    if (iblk_ec->row_ind[i] == from_blk || iblk_ec->row_ind[i] == to_blk) {
      ent = ent/2;
    }
    init_ent -= ent;
  }
  for (int32_t i=iblk_ec->row_ptr[to_blk]; i<iblk_ec->row_ptr[to_blk+1]; i++) {
    double ent = iblk_ec->val[i] * log(( (double) iblk_ec->val[i] ) / (blk_incident[to_blk] * blk_incident[ iblk_ec->row_ind[i] ]));
    if (iblk_ec->row_ind[i] == from_blk || iblk_ec->row_ind[i] == to_blk) {
      ent = ent/2;
    }
    init_ent -= ent;
  }

  double final_ent = 0;   /* Contribution to final entropy from affected rows and columns */

  /* Can walk along row for from_row as before because block must have nonzero interblock edge-count if the moving node
   * has any edges to a block */
  for (int32_t i=iblk_ec->row_ptr[from_blk]; i<iblk_ec->row_ptr[from_blk+1]; i++) {
    double ent = 0;
    int32_t iblk_edges;
    if (iblk_ec->row_ind[i] == to_blk) {
      iblk_edges = iblk_ec->val[i] - node_ngbrs[to_blk] + node_ngbrs[from_blk];
    }
    else if (iblk_ec->row_ind[i] == from_blk) {
      iblk_edges = iblk_ec->val[i] - 2 * node_ngbrs[from_blk] - self_edges;
    }
    else {
      iblk_edges = iblk_ec->val[i] - node_ngbrs[ iblk_ec->row_ind[i] ];
    }

    if (iblk_edges > 0) {
      if (iblk_ec->row_ind[i] == from_blk) {
        if (blk_incident[from_blk] != node_ngbrs[nblks]) {
          ent = 0.5 * iblk_edges * log(( (double) iblk_edges ) / ( (blk_incident[from_blk] - node_ngbrs[nblks]) * ( blk_incident[from_blk] - node_ngbrs[nblks]) ));
        }
      }
      else if (iblk_ec->row_ind[i] == to_blk) {
        if (blk_incident[from_blk] != node_ngbrs[nblks]) {
          ent = 0.5 * iblk_edges * log(( (double) iblk_edges ) / ( (blk_incident[from_blk] - node_ngbrs[nblks]) * (blk_incident[to_blk] + node_ngbrs[nblks]) ));
        }
      }
      else if (blk_incident[from_blk] != node_ngbrs[nblks]) {
        ent = iblk_edges * log(( (double) iblk_edges ) / ( (blk_incident[from_blk] - node_ngbrs[nblks]) * blk_incident[ iblk_ec->row_ind[i] ] ));
      }
    }

    final_ent -= ent;
  }

  /* Can't walk along row of iblk_ec as before because zeroes may become nonzero after moving node */
  int32_t ptr = iblk_ec->row_ptr[to_blk];
  for (int32_t i=0; i<nblks; i++) {
    int32_t iblk_edges;
    if (i == to_blk) {
      iblk_edges = 2 * node_ngbrs[i] + self_edges;
    }
    else if (i == from_blk) {
      iblk_edges = node_ngbrs[from_blk] - node_ngbrs[to_blk];
    }
    else {
      iblk_edges = node_ngbrs[i];
    }
    if (ptr < iblk_ec->row_ptr[to_blk+1]) {
      if (iblk_ec->row_ind[ptr] == i) {
        iblk_edges += iblk_ec->val[ptr++];
      }
    }

    double ent = 0;
    if (iblk_edges > 0) {
      if (i == from_blk) {
        if (blk_incident[i] != node_ngbrs[nblks]) {
          ent = 0.5 * iblk_edges * log(( (double) iblk_edges ) / ( (blk_incident[to_blk] + node_ngbrs[nblks]) * (blk_incident[i] - node_ngbrs[nblks]) ));
        }
      }
      else if (i == to_blk) {
        ent = 0.5 * iblk_edges * log(( (double) iblk_edges ) / ( (blk_incident[to_blk] + node_ngbrs[nblks]) * (blk_incident[to_blk] + node_ngbrs[nblks]) ));
      }
      else {
        ent = iblk_edges * log(( (double) iblk_edges ) / ( (blk_incident[to_blk] + node_ngbrs[nblks]) * blk_incident[i] ));
      }
    }
    final_ent -= ent;
  }

  return final_ent - init_ent;
}

/* Method for computing Metropolis-Hastings given in "Efficient Monte Carlo and greedy heuristic for the inference of stochastic block models" by Peixoto. */
double sbm_met_has(
    const sbm_dyn_csr * const iblk_ec,
    const int32_t * const blk_incident,
    const int32_t * const blk_assn,
    const int32_t * const node_ngbrs,
    const int32_t nblks,
    const int32_t from_blk,
    const int32_t to_blk)
{
  double num = 0;
  double den = 0;

  ///* Store transition probabilities to avoid recomputing */
  //double * trans_prob_old = calloc( nblks, sizeof(*trans_prob_old) );
  //double * trans_prob_new = calloc( nblks, sizeof(*trans_prob_old) );

  ///* Add terms to numerator and denominator by edge, rather than by block */
  //for (int i=0; i<nedges; i++) {
  //  int32_t ngbr_blk = blk_assn[ edge_ind[i]-1 ];    /* Need to subtract 1 when nodes are 1-indexed */

  //  if (trans_prob_old[ngbr_blk] == 0) {   /* Compute value if not already computed */
  //    trans_prob_old[ngbr_blk] = ( (double) iblk_ec->val[ ngbr_blk * nblks + to_blk ] + 1 ) / ( blk_incident_old[ ngbr_blk ] + nblks );
  //  }

  //  if (trans_prob_new[ngbr_blk] == 0) {   /* Compute value if not already computed */
  //    trans_prob_new[ngbr_blk] = ( (double) iblk_ec_new[ ngbr_blk * nblks + from_blk ] + 1) / ( blk_incident_new[ ngbr_blk ] + nblks );
  //  }

  //  /* Add values to sums */
  //  num += edge_val[i] * trans_prob_new[ngbr_blk];
  //  den += edge_val[i] * trans_prob_old[ngbr_blk];
  //}

  /*
  for (int32_t i=iblk_ec->row_ptr[to_blk]; i<iblk_ec->row_ptr[to_blk+1]; i++) {
    den += node_ngbrs[iblk_ec->row_ind[i]] * (iblk_ec->val[i] + 1) / ( (double) blk_incident[ iblk_ec->row_ind[i] ] + nblks );
  }
  */

  int32_t ptr = iblk_ec->row_ptr[to_blk];
  for (int32_t i=0; i<nblks; i++) {
    int32_t iblk_edges = 0;
    if (iblk_ec->row_ind[ptr] == i) {
      iblk_edges += iblk_ec->val[ptr];
      if (ptr < iblk_ec->row_ptr[to_blk+1] - 1) {
        ptr++;
      }
    }
    den += node_ngbrs[i] * (iblk_edges + 1) / ( (double) blk_incident[i] + nblks );
  }

  ptr = iblk_ec->row_ptr[from_blk];
  for (int32_t i=0; i<nblks; i++) {
    int32_t iblk_edges;
    if (i == from_blk) {
      iblk_edges = -2 * node_ngbrs[from_blk];
    }
    else if (i == to_blk) {
      iblk_edges = node_ngbrs[from_blk] - node_ngbrs[to_blk];
    }
    else {
      iblk_edges = node_ngbrs[to_blk];
    }
    if (ptr < iblk_ec->row_ptr[from_blk+1]) {
      if (iblk_ec->row_ind[ptr] == i) {
        iblk_edges += iblk_ec->val[ptr++];
      }
    }

    if (i == to_blk) {
      num += node_ngbrs[i] * (iblk_edges + 1) / ( (double) blk_incident[i] + node_ngbrs[nblks] + nblks );
    }
    else if (i == from_blk) {
      num += node_ngbrs[i] * (iblk_edges + 1) / ( (double) blk_incident[i] - node_ngbrs[nblks] + nblks );
    }
    else {
      num += node_ngbrs[i] * (node_ngbrs[i] + 1) / ( (double) blk_incident[i] + nblks );
    }
  }

  return num/den;
}

int32_t * sbm_partition_nblks(
    sbm_dyn_csr * adj_mat,
    const double beta)
{
  srand(1);

  int32_t nblks = adj_mat->nrows;
  int32_t * blk_assn = malloc(adj_mat->nrows * sizeof(*blk_assn));
  int32_t * opt_blk_assn = malloc(adj_mat->nrows * sizeof(*blk_assn));

  int32_t * blk_incident = malloc(nblks * sizeof(*blk_assn));

  for (int32_t i=0; i<nblks; i++) {
    blk_assn[i] = i;
  }
  int32_t * blk_moves = malloc(nblks * sizeof(*blk_moves));
  double * blk_ent_change = malloc(nblks * sizeof(*blk_ent_change));
  int32_t * blk_map = malloc(nblks * sizeof(*blk_map));
  int32_t new_nblks;

  int32_t * blk_ngbrs = malloc((nblks+1) * sizeof(*blk_ngbrs));

  double moving_delta_ent[3];

  double delta_ent = 0;
  double avg_delta_ent;

  int32_t * blk_lbl = malloc(nblks * sizeof(*blk_lbl));
  for (int32_t i=0; i<nblks; i++) {
    blk_lbl[i] = i;
  }

  double best_desc_lens[3] = {-1, -1, -1};
  int32_t best_nblks[3] = {-1, -1, -1};
  int32_t * best_blk_assn[3] = {NULL, NULL, NULL};
  sbm_dyn_csr * best_iblk_ec[3] = {NULL, NULL, NULL};

  sbm_dyn_csr * iblk_ec = sbm_compute_iblk_ec(adj_mat, blk_assn, nblks);
  sbm_compute_blk_incident_edges(blk_incident, iblk_ec, nblks);

  double desc_len = sbm_desc_len(iblk_ec, blk_incident, nblks, adj_mat->nrows);
  printf("Description Length: %0.04f\n", desc_len);
  printf("\tNumber of blocks: %d\n", nblks);

  new_nblks = nblks * CONTRACTION_FACTOR;

  while (new_nblks > 0) {

    /* Block Merges */
    double start = monotonic_seconds();
    int32_t * iblk_diag = malloc(nblks * sizeof(*iblk_diag));
    for (int32_t i=0; i<nblks; i++) {
      iblk_diag[i] = 0;
      for (int32_t j=iblk_ec->row_ptr[i]; j<iblk_ec->row_ptr[i+1]; j++) {
        if (iblk_ec->row_ind[j] == i) {
          iblk_diag[i] = iblk_ec->val[j];
        }
      }
    }

    for (int32_t i=0; i<nblks; i++) {
      blk_moves[i] = -1;
      blk_ent_change[i] = 0;
      memset(blk_ngbrs, 0, (nblks+1) * sizeof(*blk_ngbrs));
      sbm_node_ngbr_blk_count(iblk_ec, i, blk_lbl, nblks, blk_ngbrs);
      for (int32_t trial=0; trial<MAX_TRIALS; trial++) {

        int32_t merge_blk = sbm_propose_new_blk(i, iblk_ec, blk_lbl, iblk_ec, blk_incident, nblks, 1);

        delta_ent = sbm_delta_entropy(iblk_ec, blk_incident, blk_ngbrs, iblk_diag[i], nblks, blk_lbl[i], merge_blk);

        if (delta_ent < blk_ent_change[i] || trial == 0) {
          blk_ent_change[i] = delta_ent;
          blk_moves[i] = merge_blk;
        }
      }
    }

    printf("\nAttempting to reduce from %d blocks to %d blocks...\n", nblks, new_nblks);

    nblks = sbm_merge_blks(blk_moves, blk_ent_change, nblks - new_nblks, nblks, blk_map);    /* Returns actual number of blocks left */

    printf("Actual number of blocks: %d\n", nblks);

    sbm_relabel_blks(blk_assn, blk_map, adj_mat->nrows);

    sbm_dyn_csr_free(iblk_ec);
    iblk_ec = sbm_compute_iblk_ec(adj_mat, blk_assn, nblks);
    sbm_compute_blk_incident_edges(blk_incident, iblk_ec, nblks);
    block_merge_time += monotonic_seconds() - start;

    free(iblk_diag);

    memcpy(opt_blk_assn, blk_assn, adj_mat->nrows * sizeof(*blk_assn));    /* Copy initial block assignments in case no improvements are made with nodal updates */

    /* Nodal Updates */
    start = monotonic_seconds();
    int iter = 0;
    int moves = 0;
    do {
      delta_ent = 0;
      double opt_del_ent = 0;

      for (int32_t i=0; i<adj_mat->nrows; i++) {
        int32_t new_blk;

        /* Assumes node has no self-edges */
        new_blk = sbm_propose_new_blk(i, adj_mat, blk_assn, iblk_ec, blk_incident, nblks, 0);  /* Need to add 1 to first argument when nodes are 1-indexed */
        if (new_blk == blk_assn[i]) {    /* Trying to move node to same block */
          continue;
        }

        int32_t * node_ngbrs = calloc( (nblks+1), sizeof(*node_ngbrs) );
        sbm_node_ngbr_blk_count(adj_mat, i, blk_assn, nblks, node_ngbrs);

        double met_has = sbm_met_has(iblk_ec, blk_incident, blk_assn, node_ngbrs, nblks, blk_assn[i], new_blk);

        double d_ent = sbm_delta_entropy(iblk_ec, blk_incident, node_ngbrs, 0, nblks, blk_assn[i], new_blk);    /* Assumes no self-edges (0 in argument 4) */

        double rand_dbl = sbm_rand_dbl(0, 1);

        if (rand_dbl <= met_has * exp(-1*beta*d_ent)) {     /* Move is accepted */
          int32_t nedges = adj_mat->row_ptr[i+1] - adj_mat->row_ptr[i];
          sbm_update_iblk_ec(iblk_ec, blk_assn, nblks, adj_mat->row_ind + adj_mat->row_ptr[i], adj_mat->val + adj_mat->row_ptr[i], nedges, blk_assn[i], new_blk);
          sbm_update_blk_incident(blk_incident, nblks, adj_mat->val + adj_mat->row_ptr[i], nedges, blk_assn[i], new_blk);

          blk_assn[i] = new_blk;
          delta_ent += d_ent;
          moves++;
        }

        free(node_ngbrs);

      }

      if (delta_ent < opt_del_ent) {    /* Current assignment is better than previous best */
        opt_del_ent = delta_ent;
        memcpy(opt_blk_assn, blk_assn, adj_mat->nrows * sizeof(*blk_assn));
      }

      //printf("Delta entropy: %0.04f\n", delta_ent);
      moving_delta_ent[iter%3] = delta_ent;
      desc_len += delta_ent;
      iter++;

      avg_delta_ent = 0;
      if (iter >=3) {
        for (int k=0; k<3; k++) {
          avg_delta_ent += moving_delta_ent[k];
        }
        avg_delta_ent = avg_delta_ent / 3;
        //printf("Average delta entropy: %0.04f\n", avg_delta_ent);
      }
    } while ( ( (avg_delta_ent < -0.00001 * desc_len) || (iter < 3) ) && iter<20 );
    nodal_update_time += monotonic_seconds() - start;

    double desc_len = sbm_desc_len(iblk_ec, blk_incident, nblks, adj_mat->nrows);
    printf("Total nodal moves: %d\nDescription Length: %0.04f\n", moves, desc_len);

    new_nblks = sbm_new_nblks(desc_len, nblks, opt_blk_assn, iblk_ec, best_desc_lens, best_nblks, best_blk_assn, best_iblk_ec, adj_mat->nrows, CONTRACTION_FACTOR);
    if (new_nblks >= adj_mat->nrows) {
      new_nblks = adj_mat->nrows / 10;
    }
    /*
    if (new_nblks == best_nblks[0] || new_nblks == best_nblks[1]) {
      new_nblks = 0;
    }
    */
    printf("[ %d, %d, %d ]\n[ %0.02f, %0.02f, %0.02f ]\n", best_nblks[0], best_nblks[1], best_nblks[2], best_desc_lens[0], best_desc_lens[1], best_desc_lens[2]);

    /* Start with best assignments with larger number of blocks */
    if (new_nblks <= best_nblks[1] || best_nblks[0]==-1) {
      nblks = best_nblks[1];
      memcpy(blk_assn, best_blk_assn[1], adj_mat->nrows*sizeof(*blk_assn));
      if (iblk_ec != NULL) {
        sbm_dyn_csr_free(iblk_ec);
      }
      iblk_ec = sbm_dyn_csr_deep_cpy(best_iblk_ec[1]);
      sbm_compute_blk_incident_edges(blk_incident, iblk_ec, nblks);    /* Could store instead of recomputing */
    }
    else {
      nblks = best_nblks[2];
      memcpy(blk_assn, best_blk_assn[2], adj_mat->nrows*sizeof(*blk_assn));
      if (iblk_ec != NULL) {
        sbm_dyn_csr_free(iblk_ec);
      }
      iblk_ec = sbm_dyn_csr_deep_cpy(best_iblk_ec[2]);
      sbm_compute_blk_incident_edges(blk_incident, iblk_ec, nblks);
    }
  }

  memcpy(opt_blk_assn, best_blk_assn[1], adj_mat->nrows*sizeof(*opt_blk_assn));

  free(blk_assn);
  free(blk_incident);
  free(blk_moves);
  free(blk_ent_change);
  free(blk_map);
  free(blk_lbl);
  free(blk_ngbrs);

  for (int i=0; i<3; i++) {
    if (best_blk_assn[i] != NULL) {
      free(best_blk_assn[i]);
    }
    if (best_iblk_ec[i] != NULL ) {
      sbm_dyn_csr_free(best_iblk_ec[i]);
    }
  }

  sbm_dyn_csr_free(iblk_ec);

  printf("\nFlushing_time: %0.04f\n", flushing_time);
  printf("\nBlock Merge time: %0.04f\n", block_merge_time);
  printf("\nNodal Update time: %0.04f\n", nodal_update_time);

  return opt_blk_assn;
}
