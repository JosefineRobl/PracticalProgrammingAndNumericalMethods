#include <stdio.h> // Standard Input/Output, contains printf()
#include <string.h> // Contains string comparison
#include <gsl/gsl_matrix.h> // Enables Gnu Scientific Library matrices and some of their functions to be used
#include <gsl/gsl_vector.h> // Enables Gnu Scientific Library vectors and some of their functions to be used

/*
 * Prints the matrix and takes into acount if only some of the matrix is stored due to it being symmetric or antisymmetric.
 *
 * M: Pointer to gsl_matrix containing the matrix, which shall be printed.
 * matrixName: Pointer to string containing the name of the matrix, which shall be printed before the equality sign.
 * matrixType: Pointer to string containing either "symmetric upper" for an upper symmetric matrix, "symmetric lower" for a lower symmetric matrix, "antisymmetric upper" for an upper antisymmetric matrix, "antisymmetric lower" for a lower antisymmetric matrix, and another string, i.e. "normal" for a matrix matching none of the above.
 */
void printMatrix(gsl_matrix* M, char* matrixName, char* matrixType);

/*
 * Prints the vector.
 *
 * V: Pointer to gsl_vector containing the vector, which shall be printed.
 * vectorName: String containing the name of the vector.
 */
void printVector(gsl_vector* V, char* vectorName);
