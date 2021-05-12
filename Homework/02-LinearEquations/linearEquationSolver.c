#include <stdio.h> // Standard Input/Output, contains printf()
#include <stdlib.h> // Contains rand()
#include <math.h> // Contains mathematical expressions and functions
// #include <stdbool.h> // Contains the type 'bool'
// #include <string.h> // Contains string comparison
#include <gsl/gsl_vector.h> // Enables Gnu Scientific Library vectors and some of their functions to be used
#include <gsl/gsl_matrix.h> // Enables Gnu Scientific Library matrices and some of their functions to be used
#include <gsl/gsl_blas.h> // Enables Gnu Scientific Library functions for calculations with gsl vectors and matrices
#include <gsl/gsl_linalg.h> // Enables linear algebra functions from Gnu Scientific Library
#include <sys/time.h> // Measuring time taken between points in the execution
#include "qrGramSchmidt.h" // File containing the heads of the functions for the QR-decomposition, the solver and likewise
#include "printMatrix.h" // File containing the heads of the function for printing a matrix

/*
 * Printing the exercise number, letter or name as "== Exercise <exercise> ==".
 *
 * exercise: String containing the number, letter or name of the exercise
 */
void printExercise(char* exercise){
	printf("========== Exercise %s ==========\n", exercise);
}

/*
 * Printing the exercise subtitle as "----- <subtitle> -----".
 *
 * subtitle: String containing the subtitle of the exercise.
 */
void printExerciseSubtitle(char* subtitle){
	printf("----- %s -----\n", subtitle);
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
 * Fills a matrix M of dimensions (dim1, dim2) with randomly generated double values between -1 and 1.
 *
 * dim1: Integer containing the first dimension.
 * dim2: Integer containing the second dimension.
 * M: Pointer to gsl_matrix which shall be filled with randomly generated double values between -1 and 1.
 */
void fillRandomMatrix(int dim1, int dim2, gsl_matrix* M){
	for (int i = 0; i < dim1; i++) {
		for (int j = 0; j < dim2; j++) {
			// Sets the (i, j)'th element to a random double between -1 and 1
			gsl_matrix_set(M, i, j, randomBetweenPlusMinus1());
		}
	}
}

/*
 * Tests if two matrices are equal and prints if they are equal or not.
 * 
 * M1: Pointer to gsl_matrix containing the first of two matrices to be compared.
 * M1Name: String containing the name of the first matrix to be compared.
 * M2: Pointer to gsl_matrix containing the second of two matrices to be compared.
 * M2Name: String containing the name of the second matrix to be compared.
 * tolerance: Double containing the tolerance for the equality.
 */
void areMatricesEqual(gsl_matrix* M1, char* M1Name, gsl_matrix* M2, char* M2Name, double tolerance){
	// For the two matrices to be equal they must have the same dimensions
	if ((M1->size1 == M2->size1) && (M1->size2 == M2->size2)) {
		// Comparing each element of M1 with the corresponding element of M2
		for (int col = 0; col < M1->size2; col++) {
			for (int row = 0; row < M1->size1; row++) {
				// If the elements differ by more than the tolerance the loops are broken and the faliure is printed
				if (fabs(gsl_matrix_get(M1, row, col) - gsl_matrix_get(M2, row, col)) > tolerance) {
					printf("%s is not equal to %s with a tolerance of %g.\n", M1Name, M2Name, tolerance);
					// Invalidates the 'col' loop thus breaking the outer loop before next iteration
					col = M1->size2;
					// The 'row' loop is broken if a non-diagonal element is non-zero taking into acount the tolerance
					break;
				}
				// Write out the succes if not stopped when last element is processed
				if ((row == M1->size1 - 1) && (col == M1->size2 - 1)) {
					printf("%s is equal to %s with a tolerance of %g.\n", M1Name, M2Name, tolerance);
				}
			}
		}
	} else {
		printf("The two matrices have different dimensions and thus are not equal.\n");
	}
}

/*
 * Tests if two vectors are equal and prints if they are equal or not.
 *
 * V1: Pointer to gsl_vector containing the first of two vectors to be compared.
 * V1Name: String containing the name of the first vector to be compared.
 * V2: Pointer to gsl_matrix containing the second of two vectors to be compared.
 * V2Name: String containing the name of the second vector to be compared.
 * tolerance: Double containing the tolerance for the equality.
 */
void areVectorsEqual(gsl_vector* V1, char* V1Name, gsl_vector* V2, char* V2Name, double tolerance){
	// For the two vectors to be equal they must have the same dimensions
	if (V1->size == V2->size) {
		// Comparing each element of V1 with the corresponding element of V2
		for (int row = 0; row < V1->size; row++) {
			// If the elements differ by more than the tolerance the loop is broken and the faliure is printed
			if (fabs(gsl_vector_get(V1, row) - gsl_vector_get(V2, row)) > tolerance) {
				printf("%s is not equal to %s with a tolerance of %g.\n", V1Name, V2Name, tolerance);
				break;
			}
			// Write out the succes if not stopped when last element is processed
			if (row == V1->size - 1) {
				printf("%s is equal to %s with a tolerance of %g.\n", V1Name, V2Name, tolerance);
			}
		}
	} else {
		printf("The two vectors have different dimensions and this are not equal.\n");
	}
}

/*
 * Performing the QR-decomposition and testing if the QR-decomposition algorithm using the Gram-Schmidt orthogonality is implemented correct. Here the tests performed are to check that R is upper triangular, that Q^TQ = 1, and that QR = A.
 *
 * A:
 * R:
 */
void performAndTestQRGramSchmidtDecomposition(gsl_matrix* A, gsl_matrix* R){
	// Copy A before QR-decomposition for comparison
	gsl_matrix* ABeforeQRDecomposition = gsl_matrix_alloc(A->size1, A->size2); // Allocating space for the copied matrix
	gsl_matrix_memcpy(ABeforeQRDecomposition, A); // Copy the matrix A before QR-decomposition
	
	// Performing the QR-decomposition
	qrGramSchmidtDecomposition(A, R);
	
	// Printing the matrices A (in the theory called Q) and R after the QR-decomposition
	printExerciseSubtitle("Matrix A (in theory Q) after QR-decomposition");
	printMatrix(A, "Q", "normal");
	printExerciseSubtitle("Matrix R after QR-decomposition");
	printMatrix(R, "R", "normal");
	
	// Check that R is upper triangular (everything below diagonal shall be zero)
	double tolerance = 1e-4;
	for (int col = 1; col < R->size2 - 1; col++) {
		for (int row = col + 1; row < R->size1; row++) {
			// If an element below the diagonal is non-zero (more than the tolerance) the loops are broken and the faliure is printed
			if (fabs(gsl_matrix_get(R, row, col)) > tolerance) {
				printf("R is not upper triangular with a tolerance of %g.\n", tolerance);
				// Invalidates the 'col' loop thus breaking the outer loop before next iteration
				col = R->size2 - 1;
				// The 'row' loop is broken if en element in the lower triangle of matrix is non-zero taking into account the tolerance
				break;
			}
		}
	}
	printf("R is upper triangular with a tolerance of %g.\n", tolerance);
	
	// Check that Q^TQ = 1
	gsl_matrix* resultOfQTransposedTimesQ = gsl_matrix_alloc(A->size2, A->size2);
	gsl_blas_dsyrk(CblasUpper, CblasTrans, 1, A, 0, resultOfQTransposedTimesQ);
		// Only stores either upper (CblasUpper) og lower (CblasLower) part of the matrix due to it being symmetric.
	printExerciseSubtitle("Checking for Q^TQ = 1");
	printMatrix(resultOfQTransposedTimesQ, "Q^TQ", "symmetric upper");
	/* // The below is instead calculated in the function 'areMatricesEqual'
	bool conditionDiagonal, conditionNonDiagonal;
	for (int col = 0; col < R->size2; col++) {
		for (int row = 0; row < R->size1; row++) {
			// If a diagonal element is different from 1 (more than the tolerance) or a non-diagonal element is non-zero (more than the tolerance) the loops are broken an the faliure is printed
			conditionDiagonal = ((row == col) && (fabs(gsl_matrix_get(resultOfQTransposedTimesQ, row, col) - 1) > tolerance));
			conditionNonDiagonal = ((row != col) && (fabs(gsl_matrix_get(resultOfQTransposedTimesQ, row, col)) > tolerance));
			if (conditionDiagonal || conditionNonDiagonal) {
				printf("Q^TQ is not the identity matrix with a tolerance of %g.\n", tolerance);
				// Invalidates the 'col' loop thus breaking the outer loop before next iteration
				col = R->size2;
				// The 'row' loop is broken if a non-diagonal element is non-zero taking into acount the tolerance
				break;
			}
		}
	}
	*/
	gsl_matrix* identityMatrix = gsl_matrix_alloc(resultOfQTransposedTimesQ->size1, resultOfQTransposedTimesQ->size2); // Allocates for the identity matrix for comparison
	gsl_matrix_set_identity(identityMatrix); // Sets the matrix equal to the identity matrix
	areMatricesEqual(resultOfQTransposedTimesQ, "Q^TQ", identityMatrix, "the identity matrix", 1e-3);
	gsl_matrix_free(resultOfQTransposedTimesQ);
	gsl_matrix_free(identityMatrix);
	/* printf("Q^TQ is the indentity matrix with a tolerance of %g.\n", tolerance); */
	
	// Check that QR = A
	gsl_matrix* resultOfQTimesR = gsl_matrix_alloc(A->size1, R->size2);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, R, 0, resultOfQTimesR); // Multiplying Q (in the code called A) and R, while resultOfQTimesR holds the result
	printExerciseSubtitle("Checking for QR = A");
	printMatrix(resultOfQTimesR, "QR", "normal");
	printMatrix(ABeforeQRDecomposition, "A", "normal");
	areMatricesEqual(resultOfQTimesR, "QR", ABeforeQRDecomposition, "A", 1e-3);
	
	// Freeing all left matrices
	gsl_matrix_free(ABeforeQRDecomposition);
	gsl_matrix_free(resultOfQTimesR);
}

/*
 * Performin the QR-decomposition and the solves for QRx = b alongside testing if Ax = b.
 *
 * A:
 * R:
 * b:
 */
void performAndTestQRGramSchmidtSovler(gsl_matrix* A, gsl_matrix* R, gsl_vector* b){
	// Copy A before QR-decomposition for comparison
	gsl_matrix* ABeforeQRDecomposition = gsl_matrix_alloc(A->size1, A->size2); // Allocating space for the copied matrix
	gsl_matrix_memcpy(ABeforeQRDecomposition, A); // Copy the matrix A before QR-decomposition

	// Performing the QR-decomposition
	qrGramSchmidtDecomposition(A, R);
	
	// Printing the matrixes A (in the theory called Q) and R after the QR-decomposition
	printExerciseSubtitle("Matrix A (in theory Q) and R after QR-decomposition");
	printMatrix(A, "Q", "normal");
	printMatrix(R, "R", "normal");
	
	// Allocating space for the x vector
	gsl_vector* x = gsl_vector_alloc(b->size);
	
	// Solving QRx = b
	qrGramSchmidtSolve(A, R, b, x);
	printExerciseSubtitle("Vector x after using QR-solver on QRx = b");
	printVector(x, "x");
	
	// Cheking that Ax = b
	printExerciseSubtitle("Checking that Ax = b");
	gsl_vector* resultATimesX = gsl_vector_alloc(x->size); // allocating space for the vector Ax
	printf("Recalling, the random vector b generated is:\n");
	printVector(b, "b");
	printf("Multiplying our A matrix with the found x from the solver yields:\n");
	gsl_blas_dgemv(CblasNoTrans, 1, ABeforeQRDecomposition, x, 0, resultATimesX); // Multiplying A and x and saving the result in resultATimesX
	printVector(resultATimesX, "Ax");
	areVectorsEqual(resultATimesX, "Ax", b, "b", 1e-3);
	
	// Freeing all left matrices and vectors
	gsl_vector_free(x);
	gsl_vector_free(resultATimesX);
}

/*
 * ...
 *
 * A:
 * R:
 */
void performAndTestQRGramSchmidtInverse(gsl_matrix* A, gsl_matrix* R) {
	// Copy A before QR-decomposition for later calculations
	gsl_matrix* ABeforeQRDecomposition = gsl_matrix_alloc(A->size1, A->size2); // Allocating space for the copied matrix
	gsl_matrix_memcpy(ABeforeQRDecomposition, A); // Copy the matrix A before QR-decomposition
	
	// Performing the QR-decomposition
	qrGramSchmidtDecomposition(A, R);
	
	// Allocating space for the matrix B
	gsl_matrix* B = gsl_matrix_alloc(A->size1, A->size1);
	
	// Performing the inverse
	qrGramSchmidtInverse(A, R, B);
	printExerciseSubtitle("Inverse matrix B found from QR Gram-Schmidt inverse");
	printMatrix(B, "B", "normal");

	// Checking that A*B = Q*R*B = I
	gsl_matrix* resultATimesB = gsl_matrix_alloc(ABeforeQRDecomposition->size1, B->size2); // Allocates memory for the result of A*B (here A is ABeforeQRDecomposition)
	gsl_matrix* resultBTimesA = gsl_matrix_alloc(B->size1, ABeforeQRDecomposition->size2); // Allocates memory for the result of B*A (here A is ABeforeQRDecomposition)
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, ABeforeQRDecomposition, B, 0, resultATimesB); // Calculating AB = A*B
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, B, ABeforeQRDecomposition, 0, resultBTimesA); // Calculating BA = B*A
	gsl_matrix* identityMatrix = gsl_matrix_alloc(resultATimesB->size1, resultATimesB->size2); // Allocates memory for the identity matrix
	gsl_matrix_set_identity(identityMatrix); // Sets the matrix equal to the identity matrix
	printExerciseSubtitle("Checking that AB = BA = I");
	printMatrix(resultATimesB, "AB", "normal"); // Prints the resulting matrix for A*B
	printMatrix(resultBTimesA, "BA", "normal"); // Prints the resulting matrix for B*A
	areMatricesEqual(resultATimesB, "AB", identityMatrix, "the identity matrix", 1e-3);
	areMatricesEqual(resultBTimesA, "BA", identityMatrix, "the identity matrix", 1e-3);
	areMatricesEqual(resultATimesB, "AB", resultBTimesA, "BA", 1e-3);
	
	// Freeing all left matrices
	gsl_matrix_free(ABeforeQRDecomposition);
	gsl_matrix_free(B);
	gsl_matrix_free(resultATimesB);
	gsl_matrix_free(resultBTimesA);
	gsl_matrix_free(identityMatrix);
}

/*
 * The answers to exercise A.
 */
void exerciseA(void){
	// EXERCISE A PART 1
	printExercise("A, part 1");
	// Generating the A matrix
	int n = 5; // 1st dimension - (A random dimension at max 5 could be generated as n = rand()%5 + 1)
	int m = 3; // 2nd dimension - (random dimension less than 1st dimension: m = rand()%(n - 1) + 1)
	gsl_matrix* A = gsl_matrix_alloc(n, m); // Allocating space for the matrix
	// Assigning entries in A with doubles of random values ranging from -1 to 1
	fillRandomMatrix(n, m, A);
	// Aloocating space for the R matrix
	gsl_matrix* R = gsl_matrix_alloc(m, m);
	// Printing the generated matrix A
	printExerciseSubtitle("Matrix A before QR-decomposition");
	printMatrix(A, "A", "normal");
	// Performing the QR-decomposition and checks that it is correct
	performAndTestQRGramSchmidtDecomposition(A, R);
	// Freeing the allocated space for the matrices
	gsl_matrix_free(A);
	gsl_matrix_free(R);

	// EXERCISE A PART 2
	printf("\n");
	printExercise("A, part 2");
	// Generating the new square matrix A and vector b
	n = 5; // Dimension of square matrix
	A = gsl_matrix_alloc(n, n); // Allocating space for the matrix
	gsl_vector* b = gsl_vector_alloc(n); // Allocating space for the vector (same dimension as the matrix)
	// Assigning entries in A and b with doubles of random values ranging from -1 to 1
	for (int i = 0; i < n; i++) {
		gsl_vector_set(b, i, randomBetweenPlusMinus1()); // Generates random values for vector b
		for (int j = 0; j < n; j++) {
			gsl_matrix_set(A, i, j, randomBetweenPlusMinus1()); // Generates random values for matrix A
		}
	}
	// Allocating spaec for the R matrix
	R = gsl_matrix_alloc(n, n);
	// Printing the generated matrix A
	printExerciseSubtitle("Matrix A before QR-decomposition");
	printMatrix(A, "A", "normal");
	// Printing the generated vector b
	printExerciseSubtitle("Vector b before QR-decomposition");
	printVector(b, "b");
	// Performing the QR-solve and checks that it is correct
	performAndTestQRGramSchmidtSovler(A, R, b);
	// Freeing the allocated space for the matrices and vector
	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_vector_free(b);
}

/*
 * The answers to Exercise B.
 */
void exerciseB(void){
	// Printing the exercise title
	printf("\n");
	printExercise("B");
	// Generating a square matrix A
	int n = 5; // Dimension of square matrix
	gsl_matrix* A = gsl_matrix_alloc(n, n); // Allocating space for the matrix
	// Assigning entries in A with doubles of random values ranging from -1 to 1
	fillRandomMatrix(n, n, A);
	// Allocates space for the the R matrix
	gsl_matrix* R = gsl_matrix_alloc(n, n);
	// Printing the generated matrix A
	printExerciseSubtitle("Matrix A before QR-decomposition");
	printMatrix(A, "A", "normal");
	// Performing the QR-inverse and checks that it is correct
	performAndTestQRGramSchmidtInverse(A, R);
	// Freeing the allocated space for the matrices (and vector)
	gsl_matrix_free(A);
	gsl_matrix_free(R);
}

/*
 * The answers to Exercise C.
 */
void exerciseC(void){
	// Printing the exercise title
	printf("\n");
	printExercise("C");
	// Creating objects
	struct timeval start, end;
	double runtimeOwn, runtimeGSL; // Runtime for own and GSL QR-decomposition implementation respectively
	// Create new file in writing mode for inserting the running times in
	FILE* fileTemporary = fopen("temporary.txt", "w");
	// Calculating runtimes and inserting them into data.txt
	for (int dim = 1; dim < 100; dim++) {
		// Allocating space in memory for the three matrices and the vector for use in the calculations
		gsl_matrix* AOwn = gsl_matrix_alloc(dim, dim); // A for own implementation
		gsl_matrix* Agsl = gsl_matrix_alloc(dim, dim); // A for GSL implementation
		gsl_matrix* R = gsl_matrix_alloc(dim, dim); // R for own implementation
		gsl_vector* tau = gsl_vector_alloc(dim); // tau for GSL implementation
		// Filling the matrix A with random double values between -1 and 1
		fillRandomMatrix(dim, dim, AOwn);
		gsl_matrix_memcpy(Agsl, AOwn);
		// Timing own QR-decomposition
		gettimeofday(&start, NULL);
		qrGramSchmidtDecomposition(AOwn, R);
		gettimeofday(&end, NULL);
		runtimeOwn = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)*1e-6;
		// Timing GSL QR-decomposition
		gettimeofday(&start, NULL);
		gsl_linalg_QR_decomp(Agsl, tau);
		gettimeofday(&end, NULL);
		runtimeGSL = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)*1e-6;
		// Constructing a string of the runtimes using format specifiers and inserts the string into 'data.txt'
		fprintf(fileTemporary, "%d \t %f \t %f\n", dim, runtimeOwn, runtimeGSL);
		// Freeing the matrices and the vector
		gsl_matrix_free(AOwn);
		gsl_matrix_free(Agsl);
		gsl_matrix_free(R);
		gsl_vector_free(tau);
	}
	// Checks if the file 'data.txt' already exists: If it does then delete it and rename 'temporary.tmp' to 'data.txt', otherwise just rename 'temporary.tmp'. (One cannot easily overwrite files in C and this script is run multiple times, thus I do not want the values just appended the original file)
	char* filePath = "data.txt";
	FILE* file;
	if ((file = fopen(filePath, "r")) != NULL) {
		// 'data.txt' exitst, thus close and delete it
		fclose(file);
		remove(filePath);
	}
	rename("temporary.txt", filePath);
	// Printing the text for exercise C
	printf("The dimension value and the correspinding times can be seen in 'data.txt'.\n");
	printf("The plot of the runtimes can be seen in 'timeDependence.svg'.\n");
}

/*
 * The main function of the document.
 *
 * returns 0 for successful execution, and non-zero for error.
 */
int main(void){
	exerciseA();
	exerciseB();
	exerciseC();
	printf("====================\n");
	return 0;
}
