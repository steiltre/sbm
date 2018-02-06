
#include <stdlib.h>
#include <assert.h>
#include "sbm_util.h"

double sbm_rand_dbl(
    double low,
    double high)
{
  /* Create random double between 0 and 1 */
  double res = ( (double) rand() ) / (RAND_MAX+1.0);
  /* Scale and shift value before returning */
  return res * (high - low) + low;
}

/* Returns an integer in the range [0, n).
 *
 * Uses rand(), and so is affected-by/affects the same seed.
*/
int32_t randint(int32_t n) {

  int32_t res = 0;
  int32_t n_rand_max = RAND_MAX;
  if (n % (n_rand_max+1) == 0) {
    for (int32_t i=0; i<n/(n_rand_max+1) + 1; i++) {
      res += rand() * pow_int(n_rand_max+1,i);
    }
    return res;
  } else {
    for (int32_t i=0; i<n/(n_rand_max+1); i++) {
      res += rand() * pow_int(n_rand_max+1,i);
    }
    // Chop off all of the values that would cause skew...
    long end = n_rand_max / (n % (n_rand_max)); // truncate skew
    assert (end > 0L);
    end *= n;

    // ... and ignore results from rand() that fall above that limit.
    // (Worst case the loop condition should succeed 50% of the time,
    // so we can expect to bail out of this loop pretty quickly.)
    int r;
    while ((r = rand()) >= end);

    res += (r%n) * pow_int(n_rand_max+1, n/(n_rand_max+1));

    return res;
  }
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

int32_t pow_int(
    const int32_t base,
    const int32_t exp)
{
  if (exp == 0) {
    return 1;
  }
  else if (exp % 2 == 0) {
    int32_t tmp = pow_int(base, exp/2);
    return tmp * tmp;
  }
  else {
    return base * pow_int(base, exp-1);
  }
}
