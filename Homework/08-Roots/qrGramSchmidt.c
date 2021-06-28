#include <stdio.h> // Standard Input/Output, contains printf()
#include <math.h>
#include <gsl/gsl_vector.h> // Enables Gnu Scientific Library vectors and some of their functions to be used
#include <gsl/gsl_matrix.h> // Enables Gnu Scientific Library matrices and some of their functions to be used
#include <gsl/gsl_blas.h> // Enables Gnu Scientific Library functions for calculations with gsl vectors and matrices.
#include "qrGramSchmidt.h" // File containing the heads of the functions for the QR-decomposition, the solver and likewise.

// Defining dot product and norm to ease the notation
double dot(gsl_vector* x, gsl_vector* y){
	double xy;
	gsl_blas_ddot(x, y, &xy); //GSL function that calculates dot product
	return xy;
}
double norm(gsl_vector* x){
	return sqrt(dot(x, x));
}

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

/*
 * ...
 *
 * M:
 * v:
 */
void backSubstitution(gsl_matrix* M, gsl_vector* v){
	for (int i = v->size - 1; i >= 0; i--) {
		double sum = gsl_vector_get(v, i);
		for (int j = i + 1; j < v->size; j++) {
			sum -= gsl_matrix_get(M, i, j)*gsl_vector_get(v, j);
		}
		gsl_vector_set(v, i, sum/gsl_matrix_get(M, i, i));
	}
}

/*
 * ...
 *
 * Q:
 * R:
 * b:
 * x:
 */
void qrGramSchmidtSolve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){
	// Calculates Q*b and stores the solution in x
	gsl_blas_dgemv(CblasTrans, 1, Q, b, 0, x);
	// Calls the recursive function doing the backsubstitution
	backSubstitution(R, x);
}

/*
 * ...
 *
 * Q:
 * R:
 * B:
 */
void qrGramSchmidtInverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B){
	// Initializes a vector to be the basis vector
	gsl_vector* basisVector = gsl_vector_alloc(Q->size2);
	for (int i = 0; i < Q->size2; i++) {
		// Sets the basis vector
		gsl_vector_set_basis(basisVector, i);
		// Generates a view (temp. object on stack) containing the i'th column of matrix B
		gsl_vector_view col = gsl_matrix_column(B, i);
		// The i'th column in the inverse matrix B is found using the Gram-Schmidt solver (Q*R*B_i = I_i).
		qrGramSchmidtSolve(Q, R, basisVector, &col.vector);
	}
	// Freeing the allocated space for the basis vector
	gsl_vector_free(basisVector);
}

