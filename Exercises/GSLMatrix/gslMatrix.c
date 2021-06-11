#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <math.h>

#include "printMatrix.h"

void printExercise(char* exercise){
	printf("=============== Exercise %s ===============\n", exercise);
}

void printSubExercise(char* subExercise){
	printf("--------------- %s ---------------\n", subExercise);
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

int main(void){
	// EXERCISE A
	printExercise("A");
	// Allocate space for the matrices and vectors
	gsl_matrix* A = gsl_matrix_alloc(3, 3); // For the given A matrix
	gsl_matrix* ACopy = gsl_matrix_alloc(A->size1, A->size2); // For the copy of A
	gsl_vector* x = gsl_vector_alloc(A->size1); // For the vector when system A*x=b solved
	gsl_vector* b = gsl_vector_alloc(A->size1); // For the given b vector
	gsl_vector* bCalculated = gsl_vector_alloc(b->size); // For the calculation A*x when cheking correct solution
	// Adding values to matrix A
	double AValues[3][3] = {{6.13, -2.90, 5.86}, {8.08, -6.31, -3.89}, {-4.36, 1.00, 0.19}};
	for (int row = 0; row < A->size1; row++) {
		for (int col = 0; col < A->size2; col++) {
			gsl_matrix_set(A, row, col, AValues[row][col]);
		}
	}
	// Adding values to vector b
	double bValues[3] = {6.23, 5.37, 2.29};
	for (int row = 0; row < b->size; row++) {
		gsl_vector_set(b, row, bValues[row]);
	}
	// Copy A to ACopy
	gsl_matrix_memcpy(ACopy, A);
	
	printSubExercise("Subexercise I");
	// Calculates x from A and b (A*x=b)
	gsl_linalg_HH_solve(ACopy, b, x);
	printf("Using GSL Householder sovler for liniar systems yields\n");
	printVector(x, "x");

	printSubExercise("Subexercise II");
	// Checks if A*x is equal to b
	gsl_blas_dgemv(CblasNoTrans, 1, A, x, 0, bCalculated);
	printf("Now for checking if A*x (for calculated x) is equal to the given b, we now calculate c = A*x, and we will compare b and c.\n");
	printVector(bCalculated, "c");
	printVector(b, "b");
	areVectorsEqual(bCalculated, "c", b, "b", 0.001);

	// Freeing allocated space
	gsl_matrix_free(A);
	gsl_matrix_free(ACopy);
	gsl_vector_free(x);
	gsl_vector_free(b);
	gsl_vector_free(bCalculated);
	
	// EXERCISE B
	//printExercise("B");
		
	return 0;
}
