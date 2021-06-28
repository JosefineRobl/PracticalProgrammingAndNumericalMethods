#include <gsl/gsl_vector.h> // Gnu Scientific Library with vectors
#include <assert.h> // For assertions throughtout the code

/*
 * Locates the interval for z by binary searching through the entire interval. Runs O(ln(n)) time due to halfing of interval every times (traversing a binary search tree).
 *
 * n:
 * x:
 * z:
 *
 * Returns the integer i for which x[i] < z < x[i+1].
 */
int binarySearch(int n, gsl_vector* x, double z);
