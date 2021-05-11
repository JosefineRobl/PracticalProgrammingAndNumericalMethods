#include <stdio.h> // Standard Input/Output, contains printf()
#include <gsl/gsl_matrix.h> // Enables Gnu Scientific Library matrices and some of their functions to be used
#include <gsl/gsl_vector.h> // Enables Gnu Scientific Library vectors and some of their functions to be used

extern void printMatrix(gsl_matrix* M, char* matrixName, char* matrixType);

extern void printVector(gsl_vector* V, char* vectorName);
