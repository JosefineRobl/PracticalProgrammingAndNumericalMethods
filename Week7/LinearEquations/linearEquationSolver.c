#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

void gramSchmidtDecomposition(gsl_matrix* A, gsl_matrix* R){
	// Initializing m, the size of matrix R is (m, m)
	int m = A->size2;
	// Gram-Schmidt decomposition
	for (int i = 0; i < m; i++) {
		// Create a view of the i-th column of the matrix
		gsl_vector_view col = gsl_matrix_column(A, i);
		// Compute the norm of this column (sum of entries)
		gsl_matrix_set(R, i, i, gsl_blas_dnrm2(&col.vector));
		// Divide all entries of this column by its norm
		gsl_matrix_set_col(A, i, &col.vector / gsl_matrix_get(R, i, i));
		//
		for (int j = i + 1; j < m; j++) {
			// Create a view of the j-th colum of the matrix
			gsl_vector view colJ = gsl_matrix_column(A, j);
			// Compute the dot product between A_i and A_j
			gsl_matrix_set(R, i, j, gsl_blas_ddot(&col.vector, &colJ.vector));
			// Set A_j to be A_j - q_i * R_{ij}
			gsl_matrix_set_col(A, j, &colJ.vector - &col.vector * gsl_matrix_get(R, i, j));
		}
	}
}

int main(void){
	return 0;
}
