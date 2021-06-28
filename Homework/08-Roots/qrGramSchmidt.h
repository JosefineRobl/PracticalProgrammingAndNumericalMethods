#include <gsl/gsl_vector.h> // Enables Gnu Scientific Library vectors and some of their functions to be used
#include <gsl/gsl_matrix.h> // Enables Gnu Scientific Library matrices and some of their functions to be used
#include <gsl/gsl_blas.h> // Enables Gnu Scientific Library functions for calculations with gsl vectors and matrices.

double norm(gsl_vector* x)

void qrGramSchmidtDecomposition(gsl_matrix* A, gsl_matrix* R);

void qrGramSchmidtSolve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

void qrGramSchmidtInverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B);
