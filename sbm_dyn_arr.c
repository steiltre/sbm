
#include <stdlib.h>
#include <string.h>

#include "sbm_dyn_arr.h"

/*
 * @brief Find smallest power of 2 greater than value
 *
 * @param x Value to find power of 2 from
 *
 * @return Power of 2
 */
static inline int32_t smallest_pow2(
    int32_t x)
{
  if (x < 0) {
    return 0;
  }

  x--;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return x+1;
}

sbm_dyn_arr * sbm_dyn_arr_alloc()
{
  sbm_dyn_arr * arr = malloc( sizeof(*arr) );

  arr->cap = DYN_ARR_INIT_CAP;
  arr->n_elmnts = 0;

  arr->val = malloc(DYN_ARR_INIT_CAP * sizeof(*arr->val));

  return arr;
}

sbm_dyn_arr * sbm_dyn_arr_init(
    int32_t * vals,
    int32_t n_elmnts)
{
  sbm_dyn_arr * arr = malloc( sizeof(*arr) );

  arr->cap = smallest_pow2(n_elmnts);
  arr->n_elmnts = n_elmnts;

  arr->val = malloc(arr->cap * sizeof(*arr->val));

  sbm_dyn_arr_add_values(arr, vals, n_elmnts);

  return arr;
}

void sbm_dyn_arr_free(
    sbm_dyn_arr * arr)
{
  free(arr->val);
  free(arr);
}

void sbm_dyn_arr_double_capacity(
    sbm_dyn_arr * arr)
{
  arr->cap = arr->cap << 1;

  arr->val = (int32_t *) realloc(arr->val, arr->cap * sizeof(arr->val));
}

void sbm_dyn_arr_add_values(
    sbm_dyn_arr * arr,
    int32_t *vals,
    int32_t size)
{
  while( arr->n_elmnts + size > arr->cap ) {
    sbm_dyn_arr_double_capacity(arr);
  }

  memcpy( &arr->val[arr->n_elmnts], vals, size * sizeof(*vals) );

  arr->n_elmnts += size;
}
