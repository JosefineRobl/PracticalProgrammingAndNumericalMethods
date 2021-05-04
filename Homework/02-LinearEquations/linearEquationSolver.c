#include <stdio.h> // Standard Input/Output, contains printf()
#include <stdlib.h> // Contains rand()
#include <string.h> // Contains string comparison
#include <gsl/gsl_vector.h> // Enables Gnu Scientific Library vectors and some of their functions to be used
#include <gsl/gsl_matrix.h> // Enables Gnu Scientific Library matrices and some of their functions to be used
#include <gsl/gsl_blas.h> // Enables Gnu Scientific Library functions for calculations with gsl vectors and matrices.

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
 * ...
 *
 * A:
 * R:
 */
void qrGramSchmidtDecomposition(gsl_matrix* A, gsl_matrix* R){
	// Initializing m, the size of matrix R is (m, m)
	int m = A->size2;
	// Defining a double for later use for the dot product value
	double dotProduct;
	// Gram-Schmidt decomposition
	for (int i = 0; i < m; i++) {
		// Create a view of the i-th column of the matrix
		gsl_vector_view col = gsl_matrix_column(A, i);
		// Compute the norm of this column (sum of entries)
		gsl_matrix_set(R, i, i, gsl_blas_dnrm2(&col.vector));
		// Divide all entries of this column by its norm
		gsl_vector_scale(&col.vector, 1 / gsl_matrix_get(R, i, i)); // Updates &col.vector which in turn updates the i'th column of A, since gsl_vector_view points to the location of the column in the matrix
		// gsl_matrix_set_col(A, i, &col.vector); // Not needed due to above comment ???
		//
		for (int j = i + 1; j < m; j++) {
			// Create a view of the j-th colum of the matrix
			gsl_vector_view colJ = gsl_matrix_column(A, j);
			// Compute the dot product between A_i and A_j
			gsl_blas_ddot(&col.vector, &colJ.vector, &dotProduct);
			gsl_matrix_set(R, i, j, dotProduct);
			// Set A_j to be A_j - q_i * R_{ij}
			gsl_blas_daxpy(-gsl_matrix_get(R, i, j), &col.vector, &colJ.vector);
			//gsl_vector_scale(&col.vector, gsl_matrix_get(R, i, j)); // Updates &col.vector
			//gsl_vector_sub(&colJ.vector, &col.vector); // Updates &colJ.vector <(write more here, like that in the i for-loop)>
			//gsl_matrix_set_col(A, j, colJ.vector); // Same as for same line in i for-loop ???
		}
	}
}

// Solve given Ax = QRx = b
//     static public vector qr_gs_solve(matrix Q, matrix R, vector b){
//             vector x = Q.T*b;
//             return backsub(R, x);
//     }

/*
 * Prints the matrix and takes into acoount if only some of the matrix is stored due to it being symmetric or antisymmetric.
 *
 * M: Pointer to gsl_matrix containing the matrix, which shall be printed.
 * matrixType: Pointer to string containing either "symmetric upper" for an upper symmetric matrix, "symmetric lower" for a lower symmetric matrix, "antisymmetric upper" for an upper antisymmetric matrix, "antisymmetric lower" for a lower antisymmetric matrix, and another string, i.e. "normal" for a matrix matching none of the above.
 */
void printMatrix(gsl_matrix* M, char* matrixType){
	// Initializing a double for the matrix element to be kept in.
	double matrixElement;
	// Printing the start bracket of the matrix.
	printf("[");
	// Running through the rows of the matrix.
	for (int i = 0; i < M->size1; i++) {
		// Running through the columns of the matrix.
		for (int j = 0; j < M->size2; j++) {
			/*
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
					// Otherwise go to DEFAULT. This is done by giving no break for this case, thus it rolls over to subsequent CASEs through DEFAULT, but no other cases than DEFAULT will apply.
				DEFAULT
					// Just print all matrix elements as they are writting in the matrix.
					matrixElement = gsl_matrix_get(M, i, j);
			END
			*/
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
		// After each row but the last (i < M->size - 1) a newline is printed for the new line of the matrix to start on a new line.
		if (i < M->size1 - 1) {
			printf("\n");
		}
	}
	// Printing the end bracket one horizontal space from the last matrix element.
	printf(" ]\n");
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
 * ...
 *
 * Q:
 * R:
 * b:
 * x:
 */
void qrGramSchmidtSolve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

/*
 * Printing the exercise number, letter or name as "== Exercise <exercise> ==".
 *
 * exercise: String containing the number, letter or name of the exercise
 */
void printExercise(char* exercise){
	printf("========== Exercise %s ==========\n", exercise);
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
