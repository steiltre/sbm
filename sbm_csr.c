
#include <stdlib.h>

#include "sbm_csr.h"

sbm_csr * sbm_csr_alloc(
    int32_t nrows,
    int32_t nnz)
{
  sbm_csr * arr = malloc(sizeof(*arr));

  arr->nrows = nrows;
  arr->nnz = nnz;

  arr->row_ptr = malloc((nrows+1) * sizeof(*arr->row_ptr));
  arr->row_ind = malloc(nnz * sizeof(*arr->row_ind));
  arr->val = malloc(nnz * sizeof(*arr->val));

  return arr;
}

void sbm_csr_free(
    sbm_csr * arr)
{
  free(arr->row_ptr);
  free(arr->row_ind);
  free(arr->val);
  free(arr);
}

