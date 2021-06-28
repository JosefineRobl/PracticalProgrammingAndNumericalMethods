#include <gsl/gsl_vector.h> // Gnu Scientific Library with vectors
#include <assert.h> // For assertions throughtout the code
#include "binarySearch.h"
#include "linearSpline.h"

/*
 * ...
 *
 * x:
 * y:
 * z:
 *
 * Returns a double containing the spline.
 */
double linearInterpolation(gsl_vector x, gsl_vector y, double z){
	// Checks that the lengths of the two vectors are equal
	assert(x->size == y->size);
	// Finding the interval around z for the splice
	int i = binarySearch(x->size, x, z);
	// Splice (eq. 5 and 6 in lecture notes)
	double yi = gsl_vector_get(y, i);
	double yiPlus1 = gsl_vector_get(y, i+1);
	double xi = gsl_vector_get(x, i);
	double xiPlus1 = gsl_vector_get(x, i+1);
	return spline = yi + (yiPlus1 - yi) / (xiPlus1 - xi) * (z - xi);
}

/*
 * ...
 *
 * x:
 * y:
 * z:
 *
 * Returns a double containing the result of the integration of the spline on x[0] to z.
 */
double linearInterpolationIntegration(gsl_vector x, gsl_vector y, double z){
	// Checks that the lengths of the two vectors are equal
	assert(x->size == y->size);
	// Checks that z is inside the range of x
	assert(z >= gsl_vector_get(x, 0) && z <= gsl_vector_get(x, x->size - 1));
	//
	int i = 0;
	double integral = 0;
	double deltaX = 0;
	double deltaY = 0;
	while (z > gsl_vector_get(x, i + 1)) {
		deltaX = gsl_vector_get(x, i + 1) - gsl_vector_get(x, i);
		deltaY = gsl_vector_get(y, i + 1) - gsl_vector_get(y, i);
		integral += gsl_vector_get(y, i) * deltaX + 0.5 * deltaY * deltaX;
		i++;
	}
	double dX = z - gsl_vector_get(x, i);
	deltaX = gsl_vector_get(x, i + 1) - gsl_vector_get(x, i);
	deltaY = gsl_vector_get(y, i + 1) - gsl_vector_get(y, i);
	integral += gsl_vector_get(y, i) * dX + 0.5 * deltaY / deltaX * dX * dX;
	return integral;
}
