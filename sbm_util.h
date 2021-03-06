
#ifndef SBM_UTIL_H
#define SBM_UTIL_H

/* Gives us high-resolution timers. */
#define _POSIX_C_SOURCE 200809L
#include <time.h>

/*
 * @brief Generate a random double
 *
 * @param low Minimum allowable value
 * @param high Maximum allowable value
 *
 * @return Random float between low and high
 */
double sbm_rand_dbl(
    double low,
    double high);

/* @brief Generate a random int
 *
 * @param n Maximum value generated (excludes max)
 *
 * @return Randomly generated int
 */
int randint(
    int n);

/*
 * @brief Compares two doubles
 *
 * @param a First double
 * @param b Second double
 *
 * @return Value indicating if a > b
 */
int sbm_gt_dbl(
    const void * a,
    const void * b);
/**
 * * @brief Return the number of seconds since an unspecified time (e.g., Unix
 * *        epoch). This is accomplished with a high-resolution monotonic timer,
 * *        suitable for performance timing.
 * *
 * * @return The number of seconds.
 * */
static inline double monotonic_seconds()
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/*
 * @brief Compute an exponential involving integers
 *
 * @param base Base of expression
 * @param exp Exponent of expression
 *
 * @return base**exp
 */
int32_t pow_int(
    const int32_t base,
    const int32_t exp);

#endif
