#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <math.h>

/*
 * Modifying a given matrix to be the matrix itself multiplied by the Jacobi matrix from the right (A <- AJ).
 *
 * A:
 * p:
 * q:
 * theta:
 */
void timesJ(gsl_matrix* A, int p, int q, double theta){
	// Initialize the sin and cos values
	double c = cos(theta);
	double s = sin(theta);
	// Loop through the rows of matrix A
	for(int i = 0; i < A->size1; i++){
		double newAip = c*gsl_matrix_get(A, i, p) - s*gsl_matrix_get(A, i, q);
		double newAiq = s*gsl_matrix_get(A, i, p) + c*gsl_matrix_get(A, i, q);
		// Seys the ith matrix row to the updated matrix row
		gsl_matrix_set(A, i, p, newAip);
		gsl_matrix_set(A, i, q, newAiq);
	}
}

/*
 * Modifying a given matrix to be the matrix itself multiplied by the Jacobi matrix from the left (A <- JA).
 * A:
 * p:
 * q:
 * theta:
 */
void Jtimes(gsl_matrix* A, int p, int q, double theta){
	// Initialize the sin and cos values
	double c = cos(theta);
	double s = sin(theta);
	// Loop through the columns of matrix A
	for(int j=0 ; j < A->size2; j++){
		double newApj = c*gsl_matrix_get(A, p, j) + s*gsl_matrix_get(A, q, j);
		double newAqj = -s*gsl_matrix_get(A, p, i) + c*gsl_matrix_get(A, q, j);
		// Sets the jth matrix column to the updated matrix column
		gsl_matrix_set(A,p,j,newApj);
		gsl_matrix_set(A,q,j,newAqj);
	}
}
