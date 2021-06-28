#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

/*
 * Calculates the numeric gradient.
 */
void numericGradient(double f(gsl_vector* x), gsl_vector* x, gsl_vector* gradient);

/*
 * ...
 *
 * F: Objective function.
 * x: Pointer to vector; on input is the starting point, while on exit is the approximation to the root.
 * epsilon: Double containing the accuracy goal; on exit the absolute value of the gradient should be less than epsilon.
 */
void quasiNewton(double F(gsl_vector* x), gsl_vector* x, double epsilon);
