#include <gsl/gsl_vector.h> // Gnu Scientific Library with vectors
#include <assert.h> // For assertions throughtout the code
#include "binarySearch.h"

/*
 * ...
 *
 * x:
 * y:
 * z:
 *
 * Returns a double containing the spline.
 */
double linearInterpolation(gsl_vector x, gsl_vector y, double z);

/*
 * ...
 *
 * x:
 * y:
 * z:
 *
 * Returns a double containing the result of the integration of the spline on x[0] to z.
 */
double linearInterpolationIntegration(gsl_vector x, gsl_vector y, double z);
