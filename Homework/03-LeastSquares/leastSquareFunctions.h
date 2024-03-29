// General
#include <stdio.h>
#include <math.h>
#include <assert.h>
// GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

/*
 * The dot product of two GSL vectors.
 */
double cdot(gsl_vector* A, gsl_vector* B);


/*
 * Gram-Schmidt backsubstitution.
 */
void gramScmidtBackSub(gsl_matrix* R, gsl_vector* x);

/*
 * Gram-Schmidt solver.
 */
void gramSchmidtSolve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

/*
 * Gram-Schmidt decomposition
 */
void gramSchmidtDecomp(gsl_matrix* A,gsl_matrix* R);

/*
 * Gram-Schmidt inversion of matrix.
 */
void gramSchmidtInverse(gsl_matrix* A,gsl_matrix* Inv);

void matrixMultiplication(gsl_matrix* A, gsl_matrix* B,gsl_matrix* result);

void transposeMatrix(gsl_matrix* A, gsl_matrix* result);

void leastSquareFit(int m, double f(int i, double x), gsl_vector* x, gsl_vector* y, gsl_vector* dy, gsl_vector* c, gsl_matrix* S);
