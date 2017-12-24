
#include <stdint.h>

/*
 * @brief Struct for holding edges
 */
typedef struct {
  /* ID of node at beginning of edge */
  int32_t from;

  /* ID of node at end of edge */
  int32_t to;

  /* Weight of edge */
  double wgt;
} edge;

/*
 * @brief Compare from IDs of edges
 *
 * @param a Pointer to first edge
 * @param b Pointer to second edge
 *
 * @return b->from - a->from
 */
int sbm_lt_from(
    const void *a,
    const void *b);

/*
 * @brief Compare to IDs of edges
 *
 * @param a Pointer to first edge
 * @param b Pointer to second edge
 *
 * @return b->to - a->to
 */
int sbm_gt_to(
    const void *a,
    const void *b);

/*
 * @brief Compare weights of edges
 *
 * @param a POinter to first edge
 * @param b POinter to second edge
 *
 * @return b->wgt - a->wgt
 */
int sbm_gt_wgt(
    const void *a,
    const void *b);
