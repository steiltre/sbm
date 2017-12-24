
#ifndef SBM_DYN_ARR_H
#define SBM_DYN_ARR_H

#include <stdint.h>

static int32_t const DYN_ARR_INIT_CAP = 1;

/*
 * @brief Structure for holding a dynamic array
 */
typedef struct sbm_dyn_arr {
  /* Capacity of current array */
  int32_t cap;

  /* Number of elements in array */
  int32_t n_elmnts;

  /* Values in array */
  int32_t * val;
} sbm_dyn_arr;

/*
 * @brief Allocate space for a dynamic array
 *
 * @return Empty dynamic array
 */
sbm_dyn_arr * sbm_dyn_arr_alloc();

/*
 * @brief Initialize dynamic array to have space for specific number of elements
 *
 * @param vals Array of values to put in array
 * @param n_elmnts Number of elements to allocate space for
 *
 * @return Dynamic array
 */
sbm_dyn_arr * sbm_dyn_arr_init(
    int32_t * vals,
    int32_t n_elmnts);

/*
 * @brief Free memory for dynamic array
 *
 * @param arr Array to free
 */
void sbm_dyn_arr_free(
    sbm_dyn_arr * arr);

/*
 * @brief Double capacity of array
 *
 * @param arr Dynamic array
 */
void sbm_dyn_arr_double_capacity(
    sbm_dyn_arr * arr);

/*
 * @brief Add values to dynamic array
 *
 * @param arr Dynamic array to add values to
 * @param vals Array of values to add
 * @param size Number of elements to add
 */
void sbm_dyn_arr_add_values(
    sbm_dyn_arr * arr,
    int32_t *vals,
    int32_t size);

#endif
