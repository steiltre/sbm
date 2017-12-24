
#include <stdlib.h>
#include "sbm_util.h"

double sbm_rand_dbl(
    double low,
    double high)
{
  /* Create random double between 0 and 1 */
  double res = ( (double) rand() ) / RAND_MAX;
  /* Scale and shift value before returning */
  return res * (high - low) + low;
}

int sbm_gt_dbl(
    const void * a_ptr,
    const void * b_ptr)
{
  double a = *(double *) a_ptr;
  double b = *(double *) b_ptr;

  if (a < b) {
    return -1;
  }
  else if (a > b) {
    return 1;
  }
  else {
    return 0;
  }
}
