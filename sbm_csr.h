
#ifndef SBM_CSR_H
#define SBM_CSR_H

#include <stdint.h>

/*
 * Sparse CSR format array
 */
typedef struct {
  /* Number of rows in array */
  int32_t nrows;

  /* Number of nonzeroes in array */
  int32_t nnz;

  /* Pointer to beginning of rows */
  int32_t * row_ptr;

  /* Indices for row entries */
  int32_t * row_ind;

  /* Nonzero values in array */
  int32_t * val;
} sbm_csr;


/*
 * @brief Allocate a CSR array
 *
 * @param nrows Number of rows
 * @param nnz Number of nonzeroes
 *
 * @return Allocated CSR structure
 */
sbm_csr * sbm_csr_alloc(
    int32_t nrows,
    int32_t nnz);

/*
 * @brief Free memory allocated to CSR array
 *
 * @param arr CSR array to free
 */
void sbm_csr_free(
    sbm_csr * arr);

#endif
