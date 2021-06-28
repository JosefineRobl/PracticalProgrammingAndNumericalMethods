#include <stdio.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "quasiNewton.h"
double delta = sqrt (DBL_EPSILON); 

/*
 * Calculates the numeric gradient.
 */
void numericGradient(double f(gsl_vector* x), gsl_vector* x, gsl_vector* gradient){
	double fx = f(x);
	for (int i = 0; i < x->size; i++){
		double dx;
		double x_i = gsl_vector_get (x, i);
		if (fabs (x_i) < sqrt (delta)) {
			dx = delta;
		} else {
			dx = fabs (x_i)*delta;
		}
		gsl_vector_set(x, i, x_i + dx); 
		gsl_vector_set(gradient, i, (f(x) - fx)/dx);
		gsl_vector_set(x, i, x_i); 
	}
}

/*
 * ...
 *
 * F: Objective function.
 * x: Pointer to vector; on input is the starting point, while on exit is the approximation to the root.
 * epsilon: Double containing the accuracy goal; on exit the absolute value of the gradient should be less than epsilon.
 */
void quasiNewton(double F(gsl_vector* x), gsl_vector* x, double epsilon){
	// Initialization
	int n = x->size; // Dimension of vector x
	int steps = 0; // Steps taken all in all for finding the minimum
	int goodSteps = 0; // Good steps taken in search of the minimum
	int badSteps = 0; // Bad steps taken in search of the minimum
	//
	gsl_matrix* B = gsl_matrix_alloc(n, n);
	gsl_vector* gradient = gsl_vector_alloc(n);

	gsl_vector* dx = gsl_vector_alloc(n);
	gsl_vector* xNew = gsl_vector_alloc(n);
	gsl_vector* gradientNew = gsl_vector_alloc(n);

	gsl_vector* dGradient = gsl_vector_alloc(n);
	gsl_vector* u = gsl_vector_alloc(n);

	gsl_matrix_set_identity(B);
	numericGradient(f, x, gradient);

	double fx = f(x);
	double fxNew;

	while (steps < 10000) {
		steps++; // Update the steps variable
		// Perform matrix multiplication
		gsl_blas_dgemv(CblasNoTrans, -1, B, gradient, 0, dx);
		// Error handling
		if (gsl_blas_dnrm2(dx) < delta*gsl_blas_dnrm2(x)) {
			fprintf(stderr, "Error in function quasiNewton: |dx < delta*|x|.\n");
			break;
		}
		if (gsl_blas_dnrm2(gradient) < epsilon) {
			fprintf(stderr, "Error in function quasiNewton: |gradient| < epsilon.\n");
			break;
		}
		
		double lambda = 1;
		while (1) { // true til we break the loop
			gsl_vector_memcpy(xNew, x);
			gsl_vector_add(xNew, dx);
			fxNew = f(xNew);

			double sTGradient;
			gsl_blas_ddot(dx, gradient, &sTGradient);
			// Update good and bad step variable
			if (fxNew < fx + 0.01*sTGradient) {
				goodSteps++;
				break;
			}
			if (lambda < delta) {
				badSteps++;
				break;
			}
			// Update lambda
			lambda *= 0.5;
			gsl_vector_scale(dx, 0.5);
		}
		// Update the gradient vectors
		numericGradient(F, xNew, gradientNew);
		gsl_vector_memcpy(dGradient, gradientNew);
		gsl_blas_daxpy(-1, gradient , dgradient);
		// Copy the derivative due to use of vector in function
		gsl_vector_memcpy(u, dx);
		gsl_blas_dgemv(CblasNoTrans, -1, B, dGradient, 1, u);
		// Define new variables to keep dot product of gradients
		double sTDGradient, uTDGradient;
		if (fabs(sTDGradient) > 1e-12) {
			gsl_blas_ddot(u, dGradient, &uTDGradient);
			double gamma = uTDGradient/(2*sTDGradient);
			gsl_blas_daxpy (-gamma, dx, u);
			gsl_blas_dger (1./sTdgradient, u, dx, B);
			gsl_blas_dger (1./sTdgradient, dx, u, B);
		}
		// Copy the value from xNew into the original x variable, and likewise for the gradient
		gsl_vector_memcpy(x, xNew);
		gsl_vector_memcpy(gradient, gradientNew);
	}
	// Free all matrices and vectors
	gsl_matrix_free(B);
	gsl_vector_free(gradient);
	gsl_vector_free(dx);
	gsl_vector_free(xNew);
	gsl_vector_free(gradientNew);
	gsl_vector_free(dGradient);
	gsl_vector_free(u);
	// Print the number of steps, the good steps, the bad steps and the fucntion value at the minimum
	FILE* quasiNewtonFile = fopen("quasiNewtonPrints.txt", "w");
	fprintf(quasiNewtonFile, "quasi_newton : steps = %i, good steps = %i, bad steps = %i\n", steps, good_steps, bad_steps); 
	fprintf(quasiNewtonFile, "quasi_newton : function value at min %g \n", fx);
	fclose(quasiNewtonFile);
}
