
#include "graph.h"

int sbm_lt_from(
    const void *a,
    const void *b)
{
  edge * a_ptr = (edge *) a;
  edge * b_ptr = (edge *) b;

  if (b_ptr->from < a_ptr->from) {
    return 1;
  }
  else if (b_ptr->from > a_ptr->from) {
    return -1;
  }
  else {
    return a_ptr->to - b_ptr->to;
  }
}

int sbm_gt_to(
    const void *a,
    const void *b)
{
  edge * a_ptr = (edge *) a;
  edge * b_ptr = (edge *) b;

  return b_ptr->to - a_ptr->to;
}

int sbm_gt_wgt(
    const void *a,
    const void *b)
{
  edge * a_ptr = (edge *) a;
  edge * b_ptr = (edge *) b;

  if (b_ptr->wgt > a_ptr->wgt)
    return 1;
  else if (b_ptr->wgt < a_ptr->wgt)
    return -1;
  else
    return 0;
}
