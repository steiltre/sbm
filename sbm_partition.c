
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "sbm_partition.h"
#include "sbm_stack.h"
#include "sbm_util.h"
#include "sbm_file.h"
#include "pcg_variants.h"

#define MAX_TRIALS 10
#define CONTRACTION_FACTOR 0.5

double block_merge_time = 0;
double nodal_update_time = 0;
double flushing_time = 0;
double del_ent_time = 0;

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

int sbm_gt_sbm_blk_merge_info(
    const void * a_ptr,
    const void * b_ptr)
{
  sbm_blk_merge_info a = *(sbm_blk_merge_info *) a_ptr;
  sbm_blk_merge_info b = *(sbm_blk_merge_info *) b_ptr;

  if (a.del_ent < b.del_ent) {
    return -1;
  }
  else if (a.del_ent > b.del_ent) {
    return 1;
  }
  else {
    return 0;
  }
}

int sbm_lt_sbm_blk_size(
		const void * a_ptr,
		const void * b_ptr)
{
	sbm_blk_size a = *(sbm_blk_size *) a_ptr;
	sbm_blk_size b = *(sbm_blk_size *) b_ptr;

	if (a.blk_size > b.blk_size) {
		return -1;
	}
	else if (a.blk_size < b.blk_size) {
		return 1;
	}
	else {
		return 0;
	}
}

void sbm_init_assn_rand(
    int32_t * blk_assn,
    const int32_t nnodes,
    const int32_t nblks)
{
  for (int32_t i=0; i<nnodes; i++) {
    //blk_assn[i] = randint(nblks);
    blk_assn[i] = pcg32_boundedrand(nblks);
  }
}

void sbm_random_bisect(
		const int32_t blk_lbl,
		int32_t * blk_assn,
    int32_t * blk_sizes,
		const int32_t nblks,
		const int32_t nnodes)
{
	/* Scan all nodes for those in block being split. Randomly place nodes in
	 * current block or move to new block */
	for (int32_t i=0; i<nnodes; i++) {
		if (blk_assn[i] == blk_lbl) {
			double move_rnd = sbm_rand_dbl(0,1);
			if (move_rnd > 0.5) {
				blk_assn[i] = nblks;
        blk_sizes[nblks]++;
        blk_sizes[ blk_assn[i] ]--;
			}
		}
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
  if (desc_len < best_desc_lens[1] || best_desc_lens[1] == -1) {    /* New desc len is better than previous best */
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
    return ( (int) nblks * contraction_factor + 0.5);
  }
  else if (best_desc_lens[2] == -1) {    /* Still don't have bracket established (in other direction) */
    return ( (int) best_nblks[1] / contraction_factor + 0.5);
  }
  else {
    if (best_nblks[2] - best_nblks[0] == 2) {    /* Search has concluded */
      return 0;
    }
    else if (best_nblks[2] - best_nblks[1] < best_nblks[1] - best_nblks[0]) {
      return (int) (best_nblks[1] - ( 0.618 * (best_nblks[1] - best_nblks[0]) ) + 0.5 );
    }
    else {
      return (int) (best_nblks[1] + ( 0.618 * (best_nblks[2] - best_nblks[1]) ) + 0.5 );
    }
  }
}

sbm_dyn_csr *  sbm_compute_iblk_ec(
    sbm_dyn_csr * adj_mat,
    int32_t * blk_assn,
    int32_t nblks)
{
  sbm_dyn_csr * iblk_ec = sbm_dyn_csr_alloc();

  for (int32_t i=0; i<adj_mat->nrows; i++) {
    for (int32_t j=adj_mat->row_ptr[i]; j<adj_mat->row_ptr[i+1]; j++) {
      sbm_dyn_csr_add( iblk_ec, blk_assn[i], blk_assn[adj_mat->row_ind[j]], adj_mat->row_val[j] );
    }
  }
  double start = monotonic_seconds();
  sbm_dyn_csr_flush_adj_lists(iblk_ec);
  flushing_time += monotonic_seconds() - start;
  return iblk_ec;
}

void sbm_compute_blk_deg(
    int32_t * blk_in_deg,
    int32_t * blk_out_deg,
    sbm_dyn_csr * iblk_ec,
    int32_t nblks)
{
  for (int32_t i=0; i<nblks; i++) {
    blk_in_deg[i] = 0;
    blk_out_deg[i] = 0;
  }

  /* Sum up cols of iblk_ic to give blk_in_deg */
  for (int32_t i=0; i<nblks; i++) {
    for (int32_t j=iblk_ec->col_ptr[i]; j<iblk_ec->col_ptr[i+1]; j++) {
      blk_in_deg[i] += iblk_ec->col_val[j];
    }
  }

  /* Sum up rows of iblk_ic to give blk_out_deg */
  for (int32_t i=0; i<nblks; i++) {
    for (int32_t j=iblk_ec->row_ptr[i]; j<iblk_ec->row_ptr[i+1]; j++) {
      blk_out_deg[i] += iblk_ec->row_val[j];
    }
  }
}

void sbm_node_ngbr_blk_count(
    const sbm_dyn_csr * const adj_mat,
    const int32_t node,
    const int32_t * const blk_assn,
    const int32_t nblks,
    int32_t * node_ngbrs_in,
    int32_t * node_ngbrs_out)
{
  for (int32_t i=adj_mat->row_ptr[node]; i<adj_mat->row_ptr[node+1]; i++) {
    node_ngbrs_out[ blk_assn[adj_mat->row_ind[i]] ] += adj_mat->row_val[i];
    node_ngbrs_out[nblks] += adj_mat->row_val[i];
  }

  for (int32_t i=adj_mat->col_ptr[node]; i<adj_mat->col_ptr[node+1]; i++) {
    node_ngbrs_in[ blk_assn[adj_mat->col_ind[i]] ] += adj_mat->col_val[i];
    node_ngbrs_in[nblks] += adj_mat->col_val[i];
  }
}

/* Method of proposing block can be found in "Efficient Monte Carlo and greedy heuristic for the inference of stochastic block models" by Peixoto.
 * We have fixed \epsilon = 1 in this implementation */
int32_t sbm_propose_new_blk(
    const int32_t node,
    const sbm_dyn_csr * const adj_mat,
    const int32_t * const blk_assn,
    const sbm_dyn_csr * const iblk_ec,
    const int32_t * const blk_in_deg,
    const int32_t * const blk_out_deg,
    const int32_t nblks,
    const int force_unique)
{
  if (nblks == 1) {    /* Avoid dividing by zero later */
    return 0;
  }

  /* Select a random neighbor of node and get its block */
  int32_t in_nngbrs = adj_mat->row_ptr[node+1] - adj_mat->row_ptr[node];          /* Subtract 1 from node when nodes are 1-indexed */
  int32_t out_nngbrs = adj_mat->col_ptr[node+1] - adj_mat->col_ptr[node];
  int32_t nngbrs = in_nngbrs + out_nngbrs;

  if (nngbrs == 0) {    /* Return random proposal if node has no neighbors */
    int32_t prop;
    if (force_unique) {
      prop = pcg32_boundedrand(nblks-1);
      if (prop >= blk_assn[node]) {  /* Don't allow proposal to be same as original block */
        prop++;
      }
    }
    else {
      prop = pcg32_boundedrand(nblks);
    }
    return prop;
  }

  int32_t ngbr;
  int32_t ngbr_ind = pcg32_boundedrand(nngbrs);
  if (ngbr_ind < in_nngbrs) {
    ngbr = adj_mat->row_ind[ adj_mat->row_ptr[node] + ngbr_ind ];
  }
  else {
    ngbr = adj_mat->col_ind[ adj_mat->col_ptr[node] + ngbr_ind - in_nngbrs ];
  }

  int32_t blk_ngbr = blk_assn[ngbr];       /* Subtract 1 when nodes are 1-indexed */
  int32_t ngbr_blk_nngbrs = blk_in_deg[blk_ngbr] + blk_out_deg[blk_ngbr];    /* Assuming only outgoing edges matter */
  if (force_unique) {     /* Make sure old block is not considered */
    for (int32_t idx = iblk_ec->row_ptr[blk_ngbr]; idx < iblk_ec->row_ptr[blk_ngbr+1]; idx++) {
      if (iblk_ec->row_ind[idx] == blk_assn[node]) {
        ngbr_blk_nngbrs -= iblk_ec->row_val[idx];
      }
    }
    for (int32_t idx = iblk_ec->col_ptr[blk_ngbr]; idx < iblk_ec->col_ptr[blk_ngbr+1]; idx++) {
      if (iblk_ec->col_ind[idx] == blk_assn[node]) {
        ngbr_blk_nngbrs -= iblk_ec->col_val[idx];
      }
    }
  }


  /* Accept random block? */
  double accept_rand_blk_prob = ( (double) nblks) / (blk_out_deg[blk_ngbr] + blk_in_deg[blk_ngbr] + nblks);    /* Assuming only outgoing edges matter */
  double r = sbm_rand_dbl(0, 1);

  if (r < accept_rand_blk_prob || ngbr_blk_nngbrs == 0) {
    int32_t prop = -1;
    if (force_unique) {
      prop = pcg32_boundedrand(nblks-1);
      if (prop >= blk_assn[node]) {  /* Don't allow proposal to be same as original block */
        prop++;
      }
    }
    else {
      prop = pcg32_boundedrand(nblks);
    }
    return prop;
  }
  else {
    int32_t counter;
    counter = pcg32_boundedrand(ngbr_blk_nngbrs);
    int32_t prop = -1;
    for (int32_t i=iblk_ec->row_ptr[blk_ngbr]; i<iblk_ec->row_ptr[blk_ngbr+1]; i++) {
      if (iblk_ec->row_ind[i] != blk_assn[node] || !force_unique) {
        counter -= iblk_ec->row_val[i];
      }
      if (counter < 0) {
        prop = iblk_ec->row_ind[i];
        break;
      }
    }
    if (counter >= 0) {
      for (int32_t i=iblk_ec->col_ptr[blk_ngbr]; i<iblk_ec->col_ptr[blk_ngbr+1]; i++) {
        if (iblk_ec->col_ind[i] != blk_assn[node] || !force_unique) {
          counter -= iblk_ec->col_val[i];
        }
        if (counter < 0) {
          prop = iblk_ec->col_ind[i];
          break;
        }
      }
    }
    assert(prop != -1);
    return prop;
  }
}

void sbm_update_iblk_ec(
    const sbm_dyn_csr * const adj_mat,
    sbm_dyn_csr * iblk_ec,
    const int32_t * const blk_assn,
    const int32_t nblks,
    const int32_t node_id,
    const int32_t from_blk,
    const int32_t to_blk)
{
  /* Update rows of iblk_ec */
  for (int32_t i=adj_mat->row_ptr[node_id]; i<adj_mat->row_ptr[node_id+1]; i++) {
    int32_t ngbr_blk = blk_assn[ adj_mat->row_ind[i] ];
    sbm_dyn_csr_add( iblk_ec, from_blk, ngbr_blk, -1 * adj_mat->row_val[i] );
    sbm_dyn_csr_add( iblk_ec, to_blk, ngbr_blk, adj_mat->row_val[i] );
  }

  /* Update columns of iblk_ec */
  for (int32_t i=adj_mat->col_ptr[node_id]; i<adj_mat->col_ptr[node_id+1]; i++) {
    int32_t ngbr_blk = blk_assn[ adj_mat->col_ind[i] ];
    sbm_dyn_csr_add( iblk_ec, ngbr_blk, from_blk, -1 * adj_mat->col_val[i] );
    sbm_dyn_csr_add( iblk_ec, ngbr_blk, to_blk, adj_mat->col_val[i] );
  }

  double start = monotonic_seconds();
  sbm_dyn_csr_flush_adj_lists(iblk_ec);
  flushing_time += monotonic_seconds() - start;
}

void sbm_update_blk_deg(
    int32_t * blk_in_deg,
    int32_t * blk_out_deg,
    const sbm_dyn_csr * const adj_mat,
    const int32_t nblks,
    const int32_t node_id,
    const int32_t from_blk,
    const int32_t to_blk)
{
  /* Update outgoing degrees */
  for (int32_t i=adj_mat->row_ptr[node_id]; i<adj_mat->row_ptr[node_id+1]; i++) {
    blk_out_deg[from_blk] -= adj_mat->row_val[i];
    blk_out_deg[to_blk] += adj_mat->row_val[i];
  }

  /* Update incoming degrees */
  for (int32_t i=adj_mat->col_ptr[node_id]; i<adj_mat->col_ptr[node_id+1]; i++) {
    blk_in_deg[from_blk] -= adj_mat->col_val[i];
    blk_in_deg[to_blk] += adj_mat->col_val[i];
  }
}

int32_t sbm_merge_blks(
    const int32_t * const blk_moves,
    const double * const blk_del_ent,
    const int32_t blks_to_merge,
    const int32_t nblks,
    int32_t * blk_map)
{
  sbm_blk_merge_info * blk_merges = malloc(nblks * sizeof(*blk_merges));

  for (int32_t i=0; i<nblks; i++) {
    blk_merges[i].del_ent = blk_del_ent[i];
    blk_merges[i].blk_move = blk_moves[i];
    blk_merges[i].blk_lbl = i;
    blk_map[i] = i;
  }

  qsort(blk_merges, nblks, sizeof(*blk_merges), sbm_gt_sbm_blk_merge_info);

  /* Create stacks for moving already merged blocks */
  sbm_ll_node ** blk_stacks = malloc(nblks*sizeof(*blk_stacks));
  for (int32_t i=0; i<nblks; i++) {
    *(blk_stacks + i) = sbm_ll_node_alloc();
  }

  /* Create array for tracking which block labels are no longer used for relabeling at end */
  int * merged_lbls = calloc(nblks, sizeof(*merged_lbls));

  int32_t removed_blks = 0;  /* Counts number of actual merges */
  for (int32_t i=0; i<nblks; i++) {
    blk_map[ blk_merges[i].blk_lbl ] = blk_map[blk_merges[i].blk_move];
    while (blk_stacks[ blk_merges[i].blk_lbl ]->next != NULL) {    /* Drag along blocks that have already merged */
      int32_t already_moved_blk = sbm_ll_node_pop(blk_stacks[ blk_merges[i].blk_lbl ]);
      blk_map[already_moved_blk] = blk_map[ blk_merges[i].blk_lbl ];
      if ( blk_map[ blk_merges[i].blk_move ] != blk_merges[i].blk_lbl ) {    /* Don't add vertex to same stack it came from */
        sbm_ll_node_push(already_moved_blk, *(blk_stacks + blk_map[ blk_merges[i].blk_move ]));
      }
    }
    sbm_ll_node_push(blk_merges[i].blk_lbl, *(blk_stacks + blk_map[ blk_merges[i].blk_move ]));

    /* Prevent pairs of blocks that both want to merge with each other from being double counted */
    if ( merged_lbls[ blk_merges[i].blk_move ] != 1 || blk_map[ blk_merges[i].blk_move ] != blk_merges[i].blk_lbl ) {
      merged_lbls[ blk_merges[i].blk_lbl ] = 1;    /* Mark label as being unused */
      removed_blks++;
    }
    //else {    /* Fix second merge when pair are mutually trying to merge */
    //}

    if (removed_blks >= blks_to_merge) {
      break;
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

  free(blk_merges);

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

void sbm_merged_blk_sizes(
    int32_t * new_blk_sizes,
    const int32_t * const old_blk_sizes,
    const int32_t * const blk_map,
    const int32_t old_nblks)
{
  for (int32_t i=0; i<old_nblks; i++) {
    new_blk_sizes[ blk_map[i] ] += old_blk_sizes[i];
  }
}

double sbm_desc_len(
    const sbm_dyn_csr * const iblk_ec,
    const int32_t * const blk_in_deg,
    const int32_t * const blk_out_deg,
    const int32_t nblks,
    const int32_t nnodes)
{
  int32_t nedges = 0;
  for (int i=0; i<nblks; i++) {
    nedges += blk_in_deg[i];
  }

  double sum = 0;
  double ratio = ( (double) nblks * nblks ) / nedges;

  /* This is as given in Streaming Graph Challenge, but the results aren't very good. It doesn't appear
   * to penalize large numbers of blocks enough */
  sum = nedges * ( (1 + ratio)*log(1+ratio) - ratio*log(ratio) ) + nnodes * log(nblks);
  /* This is as given in the sample Python code */
  //sum = nedges * (1+ratio)*log(1+ratio) - ratio*log(ratio) + nnodes * log(nblks);

  for (int i=0; i<nblks; i++) {
    for (int j=iblk_ec->row_ptr[i]; j<iblk_ec->row_ptr[i+1]; j++) {
      double log1 = log( (double) iblk_ec->row_val[j] );
      double log2 = log( (double) blk_out_deg[i] );
      double log3 = log( (double) blk_in_deg[ iblk_ec->row_ind[j] ] );
      sum -= iblk_ec->row_val[j] * ( log1 - log2 - log3 );
    }
  }
  return sum;
}

double sbm_delta_entropy(
    const sbm_dyn_csr * const adj_mat,
    const sbm_dyn_csr * const iblk_ec,
    const int32_t * const blk_in_deg,
    const int32_t * const blk_out_deg,
    const int32_t * blk_assn,
    const int32_t node_id,
    const int32_t nblks,
    const int32_t from_blk,
    const int32_t to_blk)
{
  double start = monotonic_seconds();
  int32_t * node_blk_in = calloc(nblks, sizeof(*node_blk_in));
  int32_t * node_blk_out = calloc(nblks, sizeof(*node_blk_out));
  int32_t node_total_in = 0;
  int32_t node_total_out = 0;
  int32_t self_edges_out = 0;    /* Self-edges must be handled differently when moving edges between blocks */
  int32_t self_edges_in = 0;

  for (int32_t i=adj_mat->row_ptr[node_id]; i<adj_mat->row_ptr[node_id+1]; i++) {
    if (adj_mat->row_ind[i] == node_id) {
      self_edges_out += adj_mat->row_val[i];
    }
    else {
      node_blk_out[ blk_assn[adj_mat->row_ind[i]] ] += adj_mat->row_val[i];
      node_total_out += adj_mat->row_val[i];
    }
  }
  for (int32_t i=adj_mat->col_ptr[node_id]; i<adj_mat->col_ptr[node_id+1]; i++) {
    if (adj_mat->col_ind[i] == node_id) {
      self_edges_in += adj_mat->col_val[i];
    }
    else {
      node_blk_in[ blk_assn[adj_mat->col_ind[i]] ] += adj_mat->col_val[i];
      node_total_in += adj_mat->col_val[i];
    }
  }
  assert( self_edges_out == self_edges_in );

  int32_t to_in_deg_new = blk_in_deg[to_blk] + node_total_in + self_edges_in;
  int32_t to_out_deg_new = blk_out_deg[to_blk] + node_total_out + self_edges_in;
  int32_t from_in_deg_new = blk_in_deg[from_blk] - node_total_in - self_edges_in;
  int32_t from_out_deg_new = blk_out_deg[from_blk] - node_total_out - self_edges_in;

  double log_to_in_old = log( (double) blk_in_deg[to_blk] );
  double log_to_out_old = log( (double) blk_out_deg[to_blk] );
  double log_from_in_old = log( (double) blk_in_deg[from_blk] );
  double log_from_out_old = log( (double) blk_out_deg[from_blk] );
  double log_to_in_new = log( (double) to_in_deg_new );
  double log_to_out_new = log( (double) to_out_deg_new );
  double log_from_in_new = log( (double) from_in_deg_new);
  double log_from_out_new = log( (double) from_out_deg_new);

  double del_ent = 0;

  int32_t to_to=0, to_from=0, from_to=0, from_from=0;

	/* Had to add several statements similar to 
	 * if (to_row_pos != iblk_ec->nnz)
	 * to prevent trying to read from an empty row or column at the end of the matrix (whose values would then
	 * be stored beyond the allocated array. This can probably be avoided, or at least be less intrusive with
	 * rewriting. Reorganizing should at least help the problem.
	 */

  /* Calculate entropy change from node's outgoing edges (rows in adjacency matrix) */
	/* These terms are all of the form M_rs log(M_rs). 
	 * M_rs log(1/(d_r,out * d_s,in)) terms will be handled later */
  int32_t to_row_pos = iblk_ec->row_ptr[to_blk];
  int32_t from_row_pos = iblk_ec->row_ptr[from_blk];
  for (int32_t i=0; i<nblks; i++) {
		if (to_row_pos != iblk_ec->nnz) {
	    while ( iblk_ec->row_ind[to_row_pos] < i && to_row_pos < iblk_ec->row_ptr[to_blk+1]-1) {
				to_row_pos++;
	    }
		}
		if (from_row_pos != iblk_ec->nnz) {
	    while ( iblk_ec->row_ind[from_row_pos] < i && from_row_pos < iblk_ec->row_ptr[from_blk+1]-1) {
				from_row_pos++;
	    }
		}
    if (node_blk_out[i] != 0 && i != from_blk && i != to_blk) {    /* Node has outgoing edges that move */
      int32_t iblk_edges_orig;
      /* First consider entropy change from block node originated from */
			if (from_row_pos != iblk_ec->nnz) {
	      if ( iblk_ec->row_ind[from_row_pos] == i ) {    /* from_blk originally had edges to block being considered */
	        iblk_edges_orig = iblk_ec->row_val[from_row_pos];
	      }
	      else {
	        iblk_edges_orig = 0;
	      }
	
	      if (iblk_edges_orig != 0) {
	        del_ent += iblk_edges_orig * log( (double) iblk_edges_orig );
	      }
	      if ( iblk_edges_orig - node_blk_out[i] != 0 ) {    /* Not all inter-block edges from given node */
	        del_ent -= (iblk_edges_orig - node_blk_out[i]) * log( (double) iblk_edges_orig - node_blk_out[i] );
	      }
			}

      /* Second, consider entropy change from block node is moving to */
			if (to_row_pos != iblk_ec->nnz) {
	      if ( iblk_ec->row_ind[to_row_pos] == i ) {    /* to_blk originally had edges to block being considered */
	        iblk_edges_orig = iblk_ec->row_val[to_row_pos];
	      }
	      else {
	        iblk_edges_orig = 0;
	      }
	
	      if (iblk_edges_orig != 0) {
	        del_ent += iblk_edges_orig * log( (double) iblk_edges_orig );
	      }
	      del_ent -= (iblk_edges_orig + node_blk_out[i]) * log( (double) iblk_edges_orig + node_blk_out[i] );
			}
    }
  }

  /* Calculate entropy change from node's incoming edges (columns in adjacency matrix) */
  int32_t to_col_pos = iblk_ec->col_ptr[to_blk];
  int32_t from_col_pos = iblk_ec->col_ptr[from_blk];
  for (int32_t i=0; i<nblks; i++) {
		if (to_col_pos != iblk_ec->nnz) {
	    while (iblk_ec->col_ind[to_col_pos] < i && to_col_pos < iblk_ec->col_ptr[to_blk+1]-1) {
				to_col_pos++;
	    }
		}
		if (from_col_pos != iblk_ec->nnz) {
	    while (iblk_ec->col_ind[from_col_pos] < i && from_col_pos < iblk_ec->col_ptr[from_blk+1]-1) {
				from_col_pos++;
	    }
		}
    if (node_blk_in[i] != 0 && i != from_blk && i != to_blk) {    /* Node has outgoing edges that move */
      int32_t iblk_edges_orig;
      /* First consider entropy change from block node originated from */
			if (from_col_pos != iblk_ec->nnz) {
	      if ( iblk_ec->col_ind[from_col_pos] == i ) {    /* from_blk originally had edges to block being considered */
	        iblk_edges_orig = iblk_ec->col_val[from_col_pos];
	      }
	      else {
	        iblk_edges_orig = 0;
	      }

      	if (iblk_edges_orig != 0) {
      	  del_ent += iblk_edges_orig * log( (double) iblk_edges_orig );
      	}
      	if ( iblk_edges_orig - node_blk_in[i] != 0 ) {    /* Not all inter-block edges from given node */
      	  del_ent -= (iblk_edges_orig - node_blk_in[i]) * log( (double) iblk_edges_orig - node_blk_in[i] );
      	}
			}

      /* Second, consider entropy change from block node is moving to */
			if (to_col_pos != iblk_ec->nnz) {
	      if ( iblk_ec->col_ind[to_col_pos] == i ) {    /* to_blk originally had edges to block being considered */
	        iblk_edges_orig = iblk_ec->col_val[to_col_pos];
	      }
	      else {
	        iblk_edges_orig = 0;
	      }
	
	      if (iblk_edges_orig != 0) {
	        del_ent += iblk_edges_orig * log( (double) iblk_edges_orig );
	      }
	      del_ent -= (iblk_edges_orig + node_blk_in[i]) * log( (double) iblk_edges_orig + node_blk_in[i] );
			}
    }

    /* Store intersections of from rows and columns for later */
    if (i == from_blk) {
			if (to_col_pos != iblk_ec->nnz) {
	      if (iblk_ec->col_ind[to_col_pos] == i) {
	        from_to = iblk_ec->col_val[to_col_pos];
	      }
			}

			if (from_col_pos != iblk_ec->nnz) {
	      if (iblk_ec->col_ind[from_col_pos] == i) {
	        from_from = iblk_ec->col_val[from_col_pos];
	      }
			}
    }

    if (i == to_blk) {
			if (to_col_pos != iblk_ec->nnz) {
	      if (iblk_ec->col_ind[to_col_pos] == i) {
	        to_to = iblk_ec->col_val[to_col_pos];
	      }
			}

			if (from_col_pos != iblk_ec->nnz) {
	      if (iblk_ec->col_ind[from_col_pos] == i) {
	        to_from = iblk_ec->col_val[from_col_pos];
	      }
			}
    }
  }

	/* Handle all M_rs log( 1/(d_r,out * d_s,in) ) terms
	 * Notice that summing M_rs over rows or columns gives out and in degrees for blocks, respectively */
  if (blk_out_deg[from_blk] != 0) {
    del_ent -= blk_out_deg[from_blk] * log_from_out_old;
  }
  if (from_out_deg_new != 0) {
    del_ent += from_out_deg_new * log_from_out_new;
  }
  if (blk_out_deg[to_blk] != 0) {
    del_ent -= blk_out_deg[to_blk] * log_to_out_old;
  }
  if (to_out_deg_new != 0) {
    del_ent += to_out_deg_new * log_to_out_new;
  }
  if (blk_in_deg[from_blk] != 0) {
    del_ent -= blk_in_deg[from_blk] * log_from_in_old;
  }
  if (from_in_deg_new != 0) {
    del_ent += from_in_deg_new * log_from_in_new;
  }
  if (blk_in_deg[to_blk] != 0) {
    del_ent -= blk_in_deg[to_blk] * log_to_in_old;
  }
  if (to_in_deg_new != 0) {
    del_ent += to_in_deg_new * log_to_in_new;
  }

  /* Handle intersections of from and to rows and columns separately */
  if (from_from != 0) {
    del_ent += from_from * log( (double) from_from);
    int32_t new_from_from = from_from - node_blk_in[from_blk] - node_blk_out[from_blk] - self_edges_out;
    assert( new_from_from >= 0 );
    if (new_from_from != 0) {
      del_ent -= new_from_from * log( (double) new_from_from );
    }
  }

  if (from_to != 0) {
    del_ent += from_to * log( (double) from_to);
  }
  int32_t new_from_to = from_to + node_blk_in[from_blk] - node_blk_out[to_blk];
  assert( new_from_to >= 0 );
  if (new_from_to != 0) {
    del_ent -= new_from_to * log( (double) new_from_to);
  }

  if (to_from != 0) {
    del_ent += to_from * log( (double) to_from);
  }
  int32_t new_to_from = to_from + node_blk_out[from_blk] - node_blk_in[to_blk];
  assert( new_to_from >= 0 );
  if (new_to_from != 0) {
    del_ent -= new_to_from * log( (double) new_to_from);
  }

  if (to_to != 0) {
    del_ent += to_to * log( (double) to_to);
  }
  int32_t new_to_to = to_to + node_blk_in[to_blk] + node_blk_out[to_blk] + self_edges_out;
  assert( new_to_to >= 0 );
  if (new_to_to != 0) {
    del_ent -= new_to_to * log( (double) new_to_to);
  }

  del_ent_time += monotonic_seconds() - start;

  free(node_blk_in);
  free(node_blk_out);

  return del_ent;

}

/* Method for computing Metropolis-Hastings given in "Efficient Monte Carlo and greedy heuristic for the inference of stochastic block models" by Peixoto. */
double sbm_met_has(
    const sbm_dyn_csr * const iblk_ec,
    const int32_t * const blk_in_deg,
    const int32_t * const blk_out_deg,
    const int32_t * const blk_assn,
    const int32_t * const in_ngbrs,
    const int32_t * const out_ngbrs,
    const int32_t nblks,
    const int32_t from_blk,
    const int32_t to_blk)
{
  double num = 0;
  double den = 0;

  /* Compute denominator (terms are before proposed move) */
  int32_t * terms = malloc(nblks * sizeof(*terms));    /* Array to store numerator of transition probability for each potential neighbor block */
  /* Initialize terms to be 1 (the chosen value for epsilon) */
  for (int32_t i=0; i<nblks; i++) {
    terms[i] = 1;
  }

  for (int32_t i=iblk_ec->row_ptr[to_blk]; i<iblk_ec->row_ptr[to_blk+1]; i++) {
    terms[ iblk_ec->row_ind[i] ] += iblk_ec->row_val[i];
  }
  for (int32_t i=iblk_ec->col_ptr[to_blk]; i<iblk_ec->col_ptr[to_blk+1]; i++) {
    terms[ iblk_ec->col_ind[i] ] += iblk_ec->col_val[i];
  }

  for (int32_t i=0; i<nblks; i++) {
    den += (in_ngbrs[i] + out_ngbrs[i]) * terms[i] / ( (double) blk_in_deg[i] + blk_out_deg[i] + nblks );
  }

  /* Compute numerator (terms are after proposed move) */
  for (int32_t i=0; i<nblks; i++) {
    terms[i] = 1;
  }

  for (int32_t i=iblk_ec->row_ptr[from_blk]; i<iblk_ec->row_ptr[from_blk+1]; i++) {
    terms[ iblk_ec->row_ind[i] ] += iblk_ec->row_val[i];
  }

  for (int32_t i=iblk_ec->col_ptr[from_blk]; i<iblk_ec->col_ptr[from_blk+1]; i++) {
    terms[ iblk_ec->col_ind[i] ] += iblk_ec->col_val[i];
  }

  for (int32_t i=0; i<nblks; i++) {
    terms[i] -= in_ngbrs[i] + out_ngbrs[i];
  }

  terms[to_blk] += in_ngbrs[from_blk] + out_ngbrs[from_blk];
  terms[from_blk] -= in_ngbrs[from_blk] + out_ngbrs[from_blk];

  for (int32_t i=0; i<nblks; i++) {
    if (i != to_blk && i != from_blk) {
      num += (in_ngbrs[i] + out_ngbrs[i]) * terms[i] / ( (double) blk_in_deg[i] + blk_out_deg[i] + nblks );
    }
  }

  /* Need to handle contribution from blocks involved in move separately */
  num += (in_ngbrs[to_blk] + out_ngbrs[to_blk]) * terms[to_blk] / ( (double) blk_in_deg[to_blk] + in_ngbrs[nblks] + blk_out_deg[to_blk] + out_ngbrs[nblks] + nblks );
  num += (in_ngbrs[from_blk] + out_ngbrs[from_blk]) * terms[from_blk] / ( (double) blk_in_deg[from_blk] - in_ngbrs[nblks] + blk_out_deg[from_blk] - out_ngbrs[nblks] + nblks );

  free(terms);

  return num/den;
}

void sbm_stream_partition(
		char * ifile_prefix,
		int nfiles,
		double beta)
{
  //srand(1);
  pcg32_srandom(40, 50);

	size_t append_len = 1 + ( (int) log10( (double) nfiles ) ) + 6;
	size_t len = strlen(ifile_prefix) + append_len;
	size_t truth_len = strlen(ifile_prefix) + 19;

	char * suffix = malloc(append_len * sizeof(*suffix));
	char * fname = malloc(len * sizeof(*fname));
	char * truth_fname = malloc(truth_len * sizeof(*truth_fname));

	int32_t * blk_truth;

	sbm_dyn_csr * adj_mat = sbm_dyn_csr_alloc();

	double * partition_times = malloc(nfiles * sizeof(*partition_times));

	/* Find number of nodes in first graph to initialize block assignments */
	sprintf(suffix, "_1.tsv");
	strcpy(fname, ifile_prefix);
	strcat(fname, suffix);
	int32_t nnodes = sbm_count_nodes_tsv(fname, 1);
	int32_t nblks = 0;

	int32_t * blk_assn = malloc(sizeof(*blk_assn));
	/*
	int32_t * blk_assn = malloc(nnodes * sizeof(*blk_assn));
	for (int32_t i=0; i<nnodes; i++) {
		blk_assn[i] = i;
	}
	*/

	for (int i=0; i<nfiles; i++) {
		printf("\n\n**********************************************************************\n");
		printf("Processing file %d of %d\n", i+1, nfiles);
		/* Create name of file to open */
		sprintf(suffix, "_%d.tsv", i+1);
		strcpy(fname, ifile_prefix);
		strcat(fname, suffix);

		/* If nodes are added, add them to their own block */
		nnodes = sbm_count_nodes_tsv(fname, 1);
		if (nnodes > adj_mat->nrows) {    /* New nodes added */
			blk_assn = (int32_t *) realloc(blk_assn, nnodes * sizeof(*blk_assn));
			for (int32_t i=adj_mat->nrows; i<nnodes; i++) {
				blk_assn[i] = nblks++;
			}
		}

		sbm_read_tsv_append(fname, adj_mat, 1);

		double start = monotonic_seconds();
		nblks = sbm_partition(adj_mat, beta, blk_assn, nblks);
		partition_times[i] = monotonic_seconds() - start;

		printf("File %d finished with %d blocks\n", i, nblks);
	}

	/* Average partition times */
	double avg_partition_time = 0;
	for (int32_t i=1; i<nfiles; i++) {    /* Exclude first partitioning */
		avg_partition_time += partition_times[i];
	}
	if (nfiles > 1) {
		avg_partition_time /= (nfiles-1);
	}

	strcpy(truth_fname, ifile_prefix);
	strcat(truth_fname, "_truePartition.tsv");
	blk_truth = sbm_read_truth_tsv(truth_fname, adj_mat->nrows, 1);

	sbm_partition_metrics * prec_recall = sbm_pairwise_prec_recall(blk_assn, blk_truth, adj_mat->nrows);
  printf("\nNumber of vertices: %d\nNumber of edges: %d\n", adj_mat->nrows, adj_mat->nnz);
  printf("\nPrecision: %0.04f\n", prec_recall->prec);
  printf("Recall: %0.04f\n", prec_recall->recall);

	printf("Initial partition time: %0.04f\n", partition_times[0]);
	printf("Average partition time: %0.04f\n", avg_partition_time);

  free(suffix);
  free(fname);
  free(truth_fname);
  free(blk_truth);
  free(partition_times);
  free(blk_assn);
  free(prec_recall);

  sbm_dyn_csr_free(adj_mat);

}

int32_t sbm_partition(
    sbm_dyn_csr * adj_mat,
    const double beta,
		int32_t * blk_assn,
		int32_t nblks)
{
	double start; /* Used for timing */
	int32_t * opt_blk_assn = malloc(adj_mat->nrows * sizeof(*blk_assn));    /* Store best block assignment so far */

  int32_t * blk_in_deg = malloc(nblks * sizeof(*blk_in_deg));
  int32_t * blk_out_deg = malloc(nblks * sizeof(*blk_out_deg));

  int32_t * blk_moves = malloc(nblks * sizeof(*blk_moves));    /* Store best block merge for each block */
  double * blk_ent_change = malloc(nblks * sizeof(*blk_ent_change));    /* Store entropy changes for associated block merges */
  int32_t * blk_map = malloc(nblks * sizeof(*blk_map));    /* Gives map for assigning nodes to blocks after merge step */
  int32_t new_nblks;

  //int32_t * blk_ngbrs_in = malloc((nblks+1) * sizeof(*blk_ngbrs_in));
  //int32_t * blk_ngbrs_out = malloc((nblks+1) * sizeof(*blk_ngbrs_out));
  int32_t * blk_ngbrs_in, * blk_ngbrs_out;

  double moving_delta_ent[3];

  double delta_ent = 0;
  double avg_delta_ent;

	/* Overkill!!! Setting blk_lbl to have the largest possible number of blocks
	 * ***************************************************************************************/
  int32_t * blk_lbl = malloc(adj_mat->nrows * sizeof(*blk_lbl));
  int32_t * blk_sizes = malloc(adj_mat->nrows * sizeof(*blk_sizes));
  for (int32_t i=0; i<adj_mat->nrows; i++) {
    blk_lbl[i] = i;
    blk_sizes[i] = 1;
  }

	/* Initialize variables for Golden Ratio Search */
  double best_desc_lens[3] = {-1, -1, -1};
  int32_t best_nblks[3] = {-1, -1, -1};
  int32_t * best_blk_assn[3] = {NULL, NULL, NULL};
  sbm_dyn_csr * best_iblk_ec[3] = {NULL, NULL, NULL};

  sbm_dyn_csr * iblk_ec = sbm_compute_iblk_ec(adj_mat, blk_assn, nblks);
  sbm_compute_blk_deg(blk_in_deg, blk_out_deg, iblk_ec, nblks);

  double desc_len = sbm_desc_len(iblk_ec, blk_in_deg, blk_out_deg, nblks, adj_mat->nrows);
  printf("Description Length: %0.04f\n", desc_len);
  printf("\tNumber of blocks: %d\n", nblks);

	new_nblks = sbm_new_nblks(desc_len, nblks, blk_assn, iblk_ec, best_desc_lens, best_nblks, best_blk_assn, best_iblk_ec, adj_mat->nrows, CONTRACTION_FACTOR);

  while (new_nblks > 0) {

		if (new_nblks < nblks) {
	    /* Block Merges */
	    printf("\nMerging from %d blocks to %d blocks...\n", nblks, new_nblks);
	
	    start = monotonic_seconds();
	    int32_t * iblk_diag = malloc(nblks * sizeof(*iblk_diag));
			//blk_ngbrs_in = (int32_t *) realloc(blk_ngbrs_in, (nblks+1)*sizeof(*blk_ngbrs_in));    /* nblks may have increased since array was allocated */
			//blk_ngbrs_out = (int32_t *) realloc(blk_ngbrs_out, (nblks+1)*sizeof(*blk_ngbrs_out));
			blk_moves = (int32_t *) realloc(blk_moves, nblks * sizeof(*blk_moves));
			blk_ent_change = (double *) realloc(blk_ent_change, nblks * sizeof(*blk_ent_change));
			blk_map = (int32_t *) realloc(blk_map, nblks * sizeof(*blk_map));
	    for (int32_t i=0; i<nblks; i++) {
	      iblk_diag[i] = 0;
	      for (int32_t j=iblk_ec->row_ptr[i]; j<iblk_ec->row_ptr[i+1]; j++) {
	        if (iblk_ec->row_ind[j] == i) {
	          iblk_diag[i] = iblk_ec->row_val[j];
	        }
	      }
	    }
	
      blk_ngbrs_in = malloc( (nblks+1) * sizeof(*blk_ngbrs_in) );
      blk_ngbrs_out = malloc( (nblks+1) * sizeof(*blk_ngbrs_out) );
      for (int32_t i=0; i<nblks; i++) {
        blk_moves[i] = -1;
        blk_ent_change[i] = 0;
        memset(blk_ngbrs_in, 0, (nblks+1) * sizeof(*blk_ngbrs_in));
        memset(blk_ngbrs_out, 0, (nblks+1) * sizeof(*blk_ngbrs_out));
        sbm_node_ngbr_blk_count(iblk_ec, i, blk_lbl, nblks, blk_ngbrs_in, blk_ngbrs_out);
        for (int32_t trial=0; trial<MAX_TRIALS; trial++) {
	
          int32_t merge_blk = sbm_propose_new_blk(i, iblk_ec, blk_lbl, iblk_ec, blk_in_deg, blk_out_deg, nblks, 1);

          delta_ent = sbm_delta_entropy(iblk_ec, iblk_ec, blk_in_deg, blk_out_deg, blk_lbl, i, nblks, blk_lbl[i], merge_blk);
	
          if (delta_ent < blk_ent_change[i] || trial == 0) {
            blk_ent_change[i] = delta_ent;
            blk_moves[i] = merge_blk;
          }
        }
      }
      free(blk_ngbrs_in);
      free(blk_ngbrs_out);
	
	    new_nblks = sbm_merge_blks(blk_moves, blk_ent_change, nblks - new_nblks, nblks, blk_map);    /* Returns actual number of blocks left */

      int32_t * new_blk_sizes = malloc(nblks * sizeof(*new_blk_sizes));
	
	    sbm_relabel_blks(blk_assn, blk_map, adj_mat->nrows);
      sbm_merged_blk_sizes(new_blk_sizes, blk_sizes, blk_map, nblks);

      nblks = new_nblks;
	
	    sbm_dyn_csr_free(iblk_ec);
	    iblk_ec = sbm_compute_iblk_ec(adj_mat, blk_assn, nblks);
	    sbm_compute_blk_deg(blk_in_deg, blk_out_deg, iblk_ec, nblks);
	    block_merge_time += monotonic_seconds() - start;
	
	    free(iblk_diag);
      free(blk_sizes);
      blk_sizes = new_blk_sizes;
		}
		else {
			/* Block bisections */
			printf("\nSplitting from %d blocks to %d blocks...\n", nblks, new_nblks);

      blk_sizes = realloc(blk_sizes, new_nblks * sizeof(*blk_sizes));
      for (int32_t i=nblks; i<new_nblks; i++) {
        blk_sizes[i] = 0;
      }

			sbm_blk_size * blk_sizes_sort = malloc(nblks * sizeof(*blk_sizes_sort));
			for (int32_t i=0; i<nblks; i++) {
				blk_sizes_sort[i].blk_lbl = i;
        blk_sizes_sort[i].blk_size = 0;
				//blk_sizes_sort[i].blk_size = blk_sizes[i];
			}
			/* Compute size of each block. */
		  /* SHOULD PROBABLY MOVE THIS ELSEWHERE *********************************/
			for (int32_t i=0; i<adj_mat->nrows; i++) {
				blk_sizes_sort[ blk_assn[i] ].blk_size++;
			}

			qsort(blk_sizes_sort, nblks, sizeof(*blk_sizes_sort), sbm_lt_sbm_blk_size);

			/* Split each of the largest blocks */
			int32_t i=0;
			while(nblks < new_nblks) {
				sbm_random_bisect(blk_sizes_sort[i].blk_lbl, blk_assn, blk_sizes, nblks, adj_mat->nrows);
				nblks++;
				i++;
			}
			sbm_dyn_csr_free(iblk_ec);
	    iblk_ec = sbm_compute_iblk_ec(adj_mat, blk_assn, nblks);
			blk_in_deg = (int32_t *) realloc(blk_in_deg, nblks * sizeof(*blk_in_deg));    /* Number of blocks has grown */
			blk_out_deg = (int32_t *) realloc(blk_out_deg, nblks * sizeof(*blk_out_deg));
	    sbm_compute_blk_deg(blk_in_deg, blk_out_deg, iblk_ec, nblks);

      free(blk_sizes_sort);
		}

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

        new_blk = sbm_propose_new_blk(i, adj_mat, blk_assn, iblk_ec, blk_in_deg, blk_out_deg, nblks, 0);  /* Need to add 1 to first argument when nodes are 1-indexed */

        if (new_blk == blk_assn[i]) {    /* Trying to move node to same block */
          continue;
        }

        int32_t * node_ngbrs_in = calloc( (nblks+1), sizeof(*node_ngbrs_in) );
        int32_t * node_ngbrs_out = calloc( (nblks+1), sizeof(*node_ngbrs_out) );
        sbm_node_ngbr_blk_count(adj_mat, i, blk_assn, nblks, node_ngbrs_in, node_ngbrs_out);

        double met_has = sbm_met_has(iblk_ec, blk_in_deg, blk_out_deg, blk_assn, node_ngbrs_in, node_ngbrs_out, nblks, blk_assn[i], new_blk);

        double d_ent = sbm_delta_entropy(adj_mat, iblk_ec, blk_in_deg, blk_out_deg, blk_assn, i, nblks, blk_assn[i], new_blk);    /* Assumes no self-edges (0 in argument 4) */

        double accept_prob = met_has * exp(-1*beta*d_ent);

        double rand_dbl = sbm_rand_dbl(0, 1);

        if (rand_dbl <= accept_prob) {     /* Move is accepted */
          sbm_update_iblk_ec(adj_mat, iblk_ec, blk_assn, nblks, i, blk_assn[i], new_blk);
          sbm_update_blk_deg(blk_in_deg, blk_out_deg, adj_mat, nblks, i, blk_assn[i], new_blk);

          blk_sizes[ blk_assn[i] ]--;
          blk_sizes[ new_blk]++;
          blk_assn[i] = new_blk;
          delta_ent += d_ent;
          moves++;
        }

        free(node_ngbrs_in);
        free(node_ngbrs_out);

      }

      if (delta_ent < opt_del_ent) {    /* Current assignment is better than previous best */
        opt_del_ent = delta_ent;
        memcpy(opt_blk_assn, blk_assn, adj_mat->nrows * sizeof(*blk_assn));
      }

      moving_delta_ent[iter%3] = delta_ent;
      desc_len += delta_ent;
      iter++;

      avg_delta_ent = 0;
      if (iter >=3) {
        for (int k=0; k<3; k++) {
          avg_delta_ent += moving_delta_ent[k];
        }
        avg_delta_ent = avg_delta_ent / 3;
      }
    } while ( ( (avg_delta_ent < -0.0001 * desc_len) || (iter < 3) ) && iter<20 );
    nodal_update_time += monotonic_seconds() - start;

    double desc_len = sbm_desc_len(iblk_ec, blk_in_deg, blk_out_deg, nblks, adj_mat->nrows);
    printf("Total nodal moves: %d\nDescription Length: %0.04f\n", moves, desc_len);

    new_nblks = sbm_new_nblks(desc_len, nblks, opt_blk_assn, iblk_ec, best_desc_lens, best_nblks, best_blk_assn, best_iblk_ec, adj_mat->nrows, CONTRACTION_FACTOR);

		if (new_nblks == nblks) {    /* Trying same number of edges again */
			printf("Attempting to use same number of blocks, %d, again\n", nblks);
			break;
		}

    printf("[ %d, %d, %d ]\n[ %0.02f, %0.02f, %0.02f ]\n", best_nblks[0], best_nblks[1], best_nblks[2], best_desc_lens[0], best_desc_lens[1], best_desc_lens[2]);

    /* Start with best assignments with larger number of blocks */
    if (new_nblks <= best_nblks[1] || best_nblks[2]==-1) {
      nblks = best_nblks[1];
      memcpy(blk_assn, best_blk_assn[1], adj_mat->nrows*sizeof(*blk_assn));
      if (iblk_ec != NULL) {
        sbm_dyn_csr_free(iblk_ec);
      }
      iblk_ec = sbm_dyn_csr_deep_cpy(best_iblk_ec[1]);
      sbm_compute_blk_deg(blk_in_deg, blk_out_deg, iblk_ec, nblks);    /* Could store instead of recomputing */
    }
    else {
      nblks = best_nblks[2];
      memcpy(blk_assn, best_blk_assn[2], adj_mat->nrows*sizeof(*blk_assn));
      if (iblk_ec != NULL) {
        sbm_dyn_csr_free(iblk_ec);
      }
      iblk_ec = sbm_dyn_csr_deep_cpy(best_iblk_ec[2]);
      sbm_compute_blk_deg(blk_in_deg, blk_out_deg, iblk_ec, nblks);
    }
  }

  memcpy(opt_blk_assn, best_blk_assn[1], adj_mat->nrows*sizeof(*opt_blk_assn));

  free(blk_in_deg);
  free(blk_out_deg);
  free(blk_moves);
  free(blk_ent_change);
  free(blk_map);
  free(blk_lbl);
  //free(blk_ngbrs_in);
  //free(blk_ngbrs_out);

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
  printf("\nEntropy Change time: %0.04f\n", del_ent_time);

  memcpy(opt_blk_assn, blk_assn, adj_mat->nrows*sizeof(*opt_blk_assn));
	free(opt_blk_assn);

	return best_nblks[1];
}
