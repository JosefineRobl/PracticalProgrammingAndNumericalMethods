#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>

/*
 * ...
 *
 * PER DMITRI THE PARAMETERS ARE
 * function f: takes the input vector x and fills fx=f(x)
 * vector x: on entry -- the stating point; upon exit -- the approximation to the root
 * double epsilon: the accuracy goal; on exit the condition ||f(x)|| < epsilon should be satisfied
 */
void newton(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* x, double epsilon);
