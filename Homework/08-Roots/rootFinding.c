#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

/*
 * ...
 *
 * PER DMITRI THE PARAMETERS ARE
 * function f: takes the input vector x and fills fx=f(x)
 * vector x: on entry -- the stating point; upon exit -- the approximation to the root
 * double epsilon: the accuracy goal; on exit the condition ||f(x)|| < epsilon should be satisfied
 */
void newton(void f(gsl_vector* x, gsl_vector* fx), gsl_vectoor* x, double epsilon = 1e-3){
	// ...
	
	while (gsl_blas_dnrm2(fx) > epsilon && max(abs(delta_x)) > dx){
		for (int i = 0; i < N; i++) {
			gsl_vector_memcpy(xTemporary, x);
			gsl_vector_set(xTemporary, i, gsl_vector_get(xTemporary, i) + dx);
			dfdx = (f(xTemporary) - fx)/dx;
			for (int j = 0; j < N; j++) {
				gsl_matrix_set(J, j, i, gsl_vector_get(dfdx, j));
			}
		}
		R = new matrix(n, n);
		qr_decomp(J, R);
		delta_x = qr_solve(J, R; -1.0*fx);

		lam = 1.0;
		//...
	}
}
