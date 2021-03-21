#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

void gramSchmidtDecomposition(gsl_matrix* A, gsl_matrix* R){
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
			gsl_vector_scale(&col.vector, gsl_matrix_get(R, i, j)); // Updates &col.vector
			gsl_vector_sub(&colJ.vector, &col.vector); // Updates &colJ.vector <(write more here, like that in the i for-loop)>
			//gsl_matrix_set_col(A, j, colJ.vector); // Same as for same line in i for-loop ???
		}
	}
}

// Solve given Ax = QRx = b
//     static public vector qr_gs_solve(matrix Q, matrix R, vector b){
//             vector x = Q.T*b;
//             return backsub(R, x);
//     }

void printMatrix(gsl_matrix* M){
	for (int i = 0; i < M->size1; i++) {
		for (int j = 0; j < M->size2; j++) {
			printf("%f	", gsl_matrix_get(M, i, j));
		}
		printf("\n");
	}
}

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

void gramSchmidtSolve(gsl_matrix* Q, gsl_matrix R, gsl_vector* b, gsl_vector* x);

int main(void){
	// Making matrices
	int n = 5; // dimension
	int m = 3; // dimension
	gsl_matrix* A = gsl_matrix_alloc(n, m);
	// Assigning entries in A with floats of random values ranging from -1 to 1.
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			gsl_matrix_set(A, i, j, randomBetweenPlusMinus1());
		}
	}
	printf("Matrix A before QR-decomposition\n");
	printMatrix(A);
	gsl_matrix* R = gsl_matrix_alloc(m, m);
	gramSchmidtDecomposition(A, R);
	printf("Matrix A after QR-decomposition\n");
	printMatrix(A);
	printf("Matrix R after QR-decomposition\n");
	printMatrix(R);
	gsl_matrix_free(A);
	gsl_matrix_free(R);
	return 0;
}
