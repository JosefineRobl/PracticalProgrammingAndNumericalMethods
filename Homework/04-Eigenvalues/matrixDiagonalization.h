#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

/*
 * Modifying a given matrix to be the matrix itself multiplied by the Jacobi matrix from the right (A <- AJ).
 *
 * A:
 * p:
 * q:
 * theta:
 */
void timesJ(gsl_matrix* A, int p, int q, double theta);

/*
 * Modifying a given matrix to be the matrix itself multiplied by the Jacobi matrix from the left (A <- JA).
 * A:
 * p:
 * q:
 * theta:
 */
void Jtimes(gsl_matrix* A, int p, int q, double theta);

/*
 * Considering only upper half part (Exercise C). We write a function that do the total rotation J^T*A*J in one step, only considering the upper half of the matrix A.
 */
void JTAJ_up(gsl_matrix* A, int p, int q, double theta);

/*
 * Jacobi eigenvalue algorithm for matrix diagonalization.
 */
void jacobiDiagonalisation(gsl_matrix* A, gsl_matrix* V);

/*
 * The optimized diagonalization algorithm (Exercise C). Same approach as above but under rotation we only change the upper half of A.
 */
void jacobiDiagonalisationOptimised(gsl_matrix* A, gsl_matrix* V);
