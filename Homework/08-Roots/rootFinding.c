#include <stdio.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>

#include "qrGramSchmidt.h"

/*
 * ...
 *
 * PER DMITRI THE PARAMETERS ARE
 * function f: takes the input vector x and fills fx=f(x)
 * vector x: on entry -- the stating point; upon exit -- the approximation to the root
 * double epsilon: the accuracy goal; on exit the condition ||f(x)|| < epsilon should be satisfied
 */
void newton(void f(gsl_vector* x, gsl_vector* fx), gsl_vectoor* x, double epsilon){
	int n = x->size;
	int limit = 10000;
	int steps = 0;
	// Allocate memory for matrices and vectors
	gsl_matrix* J = gsl_matrix_alloc(n, n);
	gsl_matrix* R = gsl_matrix_alloc(n, n);
	gsl_vector* Dx = gsl_vector_alloc(n); 
	gsl_vector* fx = gsl_vector_alloc(n);
	gsl_vector* df = gsl_vector_alloc(n);
	gsl_vector* z = gsl_vector_alloc(n); 
	gsl_vector* fz = gsl_vector_alloc(n);
	//
	while (steps < limit){
		// Calculates f(x) and stores it in fx
		f(x,fx);
		// Crate Jacobian
		for (int i = 0; i < n; k++){
			double xi = gsl_vector_get(i, k);
			gsl_vector_set(x, i, xi + sqrt(DBL_EPSILON)); // Our finite difference is square root of machine epsilon
			f(x,df); // now df=f(x+sqrt())
			gsl_vector_sub(df, fx); // df = f(x + sqrt()) - f(x), as it should
			// Calculate Jacobian
			for (int k = 0; k < n; k++){
				double Jki = (gsl_vector_get(df, i)/sqrt(DBL_EPSILON)); // Finite difference
				gsl_matrix_set(J, k, i, Jki);}
			gsl_vector_set(x, i, xi);
		}
		// Solve J*Dx = -f(x) using the routines from the homework "Linear Equations" (qrGramSchmidt.h)
		qrGramSchmidtDecomposition(J, R);
		qrGramSchmidtSolve(J, R, fx, Dx); // solution of J*Dx=f(x) now stored in Dx. We want the solution of J*Dx = -f(x):
		gsl_vector_scale(Dx,-1.0);
		// Use conservative step lambda*Dx
		double lambda = 1.0; // initial value
		while (lambda > 1./64){
			gsl_vector_memcpy(z, x);
			gsl_vector_add(z, Dx);
			f(z, fz); // f(x+Dx)
			if( norm(fz) < (1-lambda/2)*norm(fx)) {
				break;
			}
			lambda *= 0.5;
			gsl_vector_scale(Dx,0.5);
		}
		gsl_vector_memcpy(x,z);
		gsl_vector_memcpy(fx,fz); // Final result
		if((norm(Dx) < sqrt(DBL_EPSILON)) || (norm(fx) < eps)) {
			break;
		}
		// Update the steps variable
		steps++;
	}
	// Free matrices and vectors
	gsl_matrix_free(J);
	gsl_matrix_free(R);
	gsl_vector_free(Dx);
	gsl_vector_free(fx);
	gsl_vector_free(df);
	gsl_vector_free(z);
	gsl_vector_free(fz);
}
