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
int binarySearch(int n, gsl_vector* x, double z){
	// Checks that z is in the interval given. Otherwise results in error
	assert(gsl_vector_get(x, 0) <= z && z <= gsl_vector_get(x, n-1));
	// Initialisation
	int i = 0, j = n - 1;
	// Binary search part
	while (j - i > 1) {
		// Finds center of interval
		int mid = (i + j) / 2;
		// Change the interval to upper og lower half depending on which part 
		(z > gsl_vector_get(x, mid)) ? i = mid : j = mid;
	}
	return i;
}
