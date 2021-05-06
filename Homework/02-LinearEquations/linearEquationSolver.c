#include <stdio.h> // Standard Input/Output, contains printf()
#include <stdlib.h> // Contains rand()
// #include <string.h> // Contains string comparison
// #include <gsl/gsl_vector.h> // Enables Gnu Scientific Library vectors and some of their functions to be used
#include <gsl/gsl_matrix.h> // Enables Gnu Scientific Library matrices and some of their functions to be used
#include <gsl/gsl_blas.h> // Enables Gnu Scientific Library functions for calculations with gsl vectors and matrices.
#include "qrGramSchmidt.h" // File containing the heads of the functions for the QR-decomposition, the solver and likewise.
#include "printMatrix.h" // File containing the heads of the function for printing a matrix.

/*
 * Printing the exercise number, letter or name as "== Exercise <exercise> ==".
 *
 * exercise: String containing the number, letter or name of the exercise
 */
void printExercise(char* exercise){
	printf("========== Exercise %s ==========\n", exercise);
}

/*
 * To get a random double between -1 and +1.
 *
 * returns a random double in the interval [1, +1].
 */
double randomBetweenPlusMinus1(){
	// Initialize
	double randomDouble;
	// For a random seed (random due to dependence on computer's internal clock).
	//srand(time(NULL));
	// Generating random double value from random int
	randomDouble = (double)rand()/RAND_MAX*2.0-1.0; // Double in range -1 to 1
	// Returning the random double value between -1 and 1
	return randomDouble;
}

/*
 * Performing the QR-decomposition and testing if the QR-decomposition algorithm using the Gram-Schmidt orthogonality is implemented correct. Here the tests performed are to check that R is upper triangular, that Q^TQ = 1, and that QR = A.
 *
 * A:
 * R:
 */
void performAndTestQRGramSchmidtDecomposition(gsl_matrix* A, gsl_matrix* R){
	// Copy A before QR-decomposition for comparison
	gsl_matrix* ABeforeQRDecomposition = gsl_matrix_alloc(A->size1, A->size2); // Allocating space for the copied matrix.
	gsl_matrix_memcpy(ABeforeQRDecomposition, A); // Copy the matrix A before QR-decomposition.
	// Performing the QR-decomposition
	qrGramSchmidtDecomposition(A, R);
	// Printing the matrices A (in the theory called Q) and R after the QR-decomposition.
	printf("----- Matrix A after QR-decomposition -----\n");
	printMatrix(A, "normal");
	printf("----- Matrix R after QR-decomposition -----\n");
	printMatrix(R, "normal");
	// Check that R is upper triangular
	// ...
	// Check that Q^TQ = 1
	gsl_matrix* resultOfQTransposedTimesQ = gsl_matrix_alloc(A->size2, A->size2);
	gsl_blas_dsyrk(CblasUpper, CblasTrans, 1, A, 0, resultOfQTransposedTimesQ);
		// For filling out lower part of result-matrix, since it only stores either upper (CblasUpper) og lower (CblasLower) part of the matrix due to it being symmetric.
	
	printf("----- Calculatoin of Q^TQ -----\n");
	printMatrix(resultOfQTransposedTimesQ, "symmetric upper");
	gsl_matrix_free(resultOfQTransposedTimesQ);
	// Check that QR = A
	// ...
	// Freeing all left 
}

/*
 * The main function of the document.
 *
 * returns 0 for succesfull execution, and non-zero for error.
 */
int main(void){
	// EXERCISE A
	printExercise("A");
	// Generating the A matrix
	int n = 5; // 1st dimension
	int m = 3; // 2nd dimension
	gsl_matrix* A = gsl_matrix_alloc(n, m); // Allocating space for the matrix
	// Assigning entries in A with floats of random values ranging from -1 to 1.
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			gsl_matrix_set(A, i, j, randomBetweenPlusMinus1());
		}
	}
	// Generating the R matrix
	gsl_matrix* R = gsl_matrix_alloc(m, m);
	// Printing the 
	printf("----- Matrix A before QR-decomposition -----\n");
	printMatrix(A, "normal");
	// Performing the QR-decomposition and checks that it is correct.
	performAndTestQRGramSchmidtDecomposition(A, R);
	// Freeing the allocated space for the matrices
	gsl_matrix_free(A);
	gsl_matrix_free(R);

	// EXERCISE B
	printExercise("B");
	// <SOME MORE HERE>

	// EXERCISE C
	printExercise("C");
	// <SOME MORE HERE>
	return 0;
}
