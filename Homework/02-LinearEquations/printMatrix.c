#include <stdio.h> // Standard Input/Output, contains printf()
#include <string.h> // Contains string comparison
#include <gsl/gsl_matrix.h> // Enables Gnu Scientific Library matrices and some of their functions to be used
#include "printMatrix.h" // File containing the heads of the function for printing a matrix.

/*
 * Switch cases with string comparison
 */
#ifndef SWITCH_CASE_INIT
#define SWITCH_CASE_INIT
	#define SWITCH(X) for (char* __switch_p__ = X, __switch_next__= 1 ; __switch_p__ ; __switch_p__= 0, __switch_next__= 1) {{
	#define CASE(X)	} if (!__switch_next__ || !(__switch_next__ = strcmp(__switch_p__, X))) {
	#define DEFAULT } {
	#define END }}
#endif

/*
 * Prints the matrix and takes into acount if only some of the matrix is stored due to it being symmetric or antisymmetric.
 *
 * M: Pointer to gsl_matrix containing the matrix, which shall be printed.
 * matrixName: Pointer to string containing the name of the matrix, which shall be printed before the equality sign.
 * matrixType: Pointer to string containing either "symmetric upper" for an upper symmetric matrix, "symmetric lower" for a lower symmetric matrix, "antisymmetric upper" for an upper antisymmetric matrix, "antisymmetric lower" for a lower antisymmetric matrix, and another string, i.e. "normal" for a matrix matching none of the above.
 */
void printMatrix(gsl_matrix* M, char* matrixName, char* matrixType){
	// Initializing a double for the matrix element to be kept in.
	double matrixElement;
	// Printing the matrix name before equality sign
	prinf("%s =\n", matrixName)
	// Printing the start bracket of the matrix.
	printf("[");
	// Running through the rows of the matrix.
	for (int i = 0; i < M->size1; i++) {
		// Running through the columns of the matrix.
		for (int j = 0; j < M->size2; j++) {
			// The following cases ensures to print the entire matrix since some functions may only store the upper or lower part of the matrix if it is symmetric or antisymmetric.
			SWITCH (matrixType)
				CASE ("symmetric upper")
					if (i > j) {
						matrixElement = gsl_matrix_get(M, j, i);
					}
					// Otherwise go to DEFAULT. This is done by giving no break for this case, thus it rolls over to subsequent CASEs through DEFAULT, but no other cases than DEFAULT will apply.
				CASE ("symmetric lower")
					if (i < j) {
						matrixElement = gsl_matrix_get(M, j, i);
					}
					// Otherwise go to DEFAULT. This is done by giving no break for this case, thus it rolls over to subsequent CASEs through DEFAULT, but no other cases than DEFAULT will apply.
				CASE ("antisymmetric upper")
					if (i > j) {
						matrixElement = - gsl_matrix_get(M, j, i);
					}
					// Otherwise go to DEFAULT. This is done by giving no break for this case, thus it rolls over to subsequent CASEs through DEFAULT, but no other cases than DEFAULT will apply.
				CASE ("antisymmetric lower")
					if (i < j) {
						matrixElement = - gsl_matrix_get(M, j, i);
					}
					// Otherwise go to DEFAULT. This is done by giving no break for this case, thus it rolls over to subsequent CASEs through DEFAULT, but no other cases than DEFAULT will apply.
				DEFAULT
					// Just print all matrix elements as they are writting in the matrix.
					matrixElement = gsl_matrix_get(M, i, j);
			END
			matrixElement = gsl_matrix_get(M, i, j);
			// Adding an extra horizontal space before first element (j == 0) of every row but the first (i > 0) for aligment of the matrix elements when printed.
			if ((j == 0) && (i > 0)) {
				printf(" ");
			}
			// Printing the actual matrix element with 10 digits. Space added in the beginning for space between the matrix elements.
			printf(" %10f", matrixElement);
			// Comma added only in-between matrix elements, thus giving better distinguishability.
			if (j < M->size2 - 1) {
				printf(",");
			}
		}
		// After each row but the last (i < M->size1 - 1) a newline is printed for the new line of the matrix to start on a new line.
		if (i < M->size1 - 1) {
			printf("\n");
		}
	}
	// Printing the end bracket one horizontal space from the last matrix element.
	printf(" ]\n");
}

/*
 * Prints the vector.
 *
 * V: Pointer to gsl_vector containing the vector, which shall be printed.
 */
void printVector(gsl_vector* V){
	// Printing the start bracket of the matrix.
	printf("[");
	// Running through the rows of the matrix.
	for (int i = 0; i < V->size; i++) {
		// Adding an extra horizontal space before the element of every row but the first (i > 0) for aligment of the matrix elements when printed.
		if (i > 0) {
			printf(" ");
		}
		// Printing the actual element with 10 digits. Space added in the beginning for extra spacing.
		printf(" %10f", gsl_vector_get(V, i));
		// After each row but the last (i < V->size - 1) a newline is printed for the new line of the vector to start on a new line.
		if (i < V->size - 1) {
			printf("\n");
		}
	}
	// Printing the end bracket one horizontal space from the last vector element.
	printf(" ]\n");
}
