#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "matrixDiagonalization.h"

/*
 * Modifying a given matrix to be the matrix itself multiplied by the Jacobi matrix from the right (A <- AJ).
 *
 * A:
 * p:
 * q:
 * theta:
 */
void timesJ(gsl_matrix* A, int p, int q, double theta){
	double c = cos(theta);
	double s = sin(theta);
	int n = A->size1;
	for(int i = 0; i < n; i++){
		double aip = gsl_matrix_get(A,i,p);
		double aiq = gsl_matrix_get(A,i,q);
		double newAip = c*aip - s*aiq; 
		double newAiq = s*aip + c*aiq;
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
	double c = cos(theta);
	double s = sin(theta);
	int m = A->size2;
	for(int j = 0; j < m; j++){
		double apj = gsl_matrix_get(A,p,j);
		double aqj = gsl_matrix_get(A,q,j);
		double newApj = c*apj + s*aqj;
		double newAqj = -s*apj + c*aqj;
		gsl_matrix_set(A, p, j, newApj);
		gsl_matrix_set(A, q, j, newAqj);
	}
}

/*
 * Considering only upper half part (Exercise C). We write a function that do the total rotation J^T*A*J in one step, only considering the upper half of the matrix A.
 */
void JTAJ_up(gsl_matrix* A, int p, int q, double theta){
	double c = cos(theta);
	double s = sin(theta);
	int n = A -> size1;
	double app = gsl_matrix_get(A,p,p);
	double aqq = gsl_matrix_get(A,q,q);
	double apq = gsl_matrix_get(A,p,q);
	// Updating app, aqq and apq according to Eq. 10, i.e. by construction the angle zeroes apq after rotation
	gsl_matrix_set(A, p, q, 0);
	gsl_matrix_set(A, p, p, c*c*app - 2*s*c*apq + s*s*aqq);
	gsl_matrix_set(A, q, q, s*s*app + 2*s*c*apq + c*c*aqq);
	// Updating the rest of upper half elements. Because of the structure of J, we only have to consider the upper half of column p and q and the right part of row p and row q
	for (int i = 0; i < p; i++){
		double aip = gsl_matrix_get(A, i, p);
		double aiq = gsl_matrix_get(A, i, q);
		gsl_matrix_set(A, i, p, c*aip - s*aiq); // Here we ake use of the symmetry of A (A is symmetric)
		gsl_matrix_set(A, i, q, s*aip + c*aiq);
	}
	for (int i = p + 1; i < q; i++){
		double api = gsl_matrix_get(A,p,i);
		double aiq = gsl_matrix_get(A,i,q);
		gsl_matrix_set(A,p,i,c*api-s*aiq);
		gsl_matrix_set(A,i,q,s*api+c*aiq);
	}
	for (int i = q + 1; i < n; i++){
		double api = gsl_matrix_get(A,p,i);
		double aqi = gsl_matrix_get(A,q,i);
		gsl_matrix_set(A,p,i,c*api-s*aqi);
		gsl_matrix_set(A,q,i,s*api+c*aqi);
	}
}

/*
 * Jacobi eigenvalue algorithm for matrix diagonalization.
 */
void jacobiDiagonalisation(gsl_matrix* A, gsl_matrix* V){
	int n = A->size1;
	int changed;
	do {
		changed = 0;
		for (int p = 0; p < n-1; p++){
			for (int q = p + 1; q < n; q++){
				// The relevant elements in the matrix are
				double apq = gsl_matrix_get(A,p,q);
				double app = gsl_matrix_get(A,p,p);
				double aqq = gsl_matrix_get(A,q,q);
				// Rotation angle
				double theta = 1./2 * atan2(2*apq, aqq - app);
				double c = cos(theta);
				double s = sin(theta);
				// Calculation new elements 
				double newApp = c*c*app - 2*s*c*apq + s*s*aqq;
				double newAqq = s*s*app + 2*s*c*apq + c*c*aqq;
				// Check if the new element differ from the old one
				if ((newApp != app) || (newAqq != aqq)) {
				// If they differ we have not reached convergence and we must rotate again
					changed = 1;
					timesJ(A, p, q, theta);
					Jtimes(A, p, q, -theta); // Note that J^T(theta) = J(-theta)
					timesJ(V, p, q, theta);
				}
			}
		}
	} while (changed != 0);
}

/*
 * The optimized diagonalization algorithm (Exercise C). Same approach as above but under rotation we only change the upper half of A.
 */
void jacobiDiagonalisationOptimised(gsl_matrix* A, gsl_matrix* V){
	int n = A -> size1;
	int changed;
	do {
		changed = 0;
		for(int p = 0; p < n-1; p++){
			for(int q = p + 1; q < n; q++){
				double apq = gsl_matrix_get(A, p, q);
				double app = gsl_matrix_get(A, p, p);
				double aqq = gsl_matrix_get(A, q, q);
				double theta = 1./2 * atan2(2*apq, aqq - app);
				double c = cos(theta);
				double s = sin(theta);
				double newApp = c*c*app - 2*s*c*apq + s*s*aqq;
				double newAqq = s*s*app + 2*s*c*apq + c*c*aqq;
				if((newApp != app) || (newAqq != aqq)) {
					changed = 1;
					JTAJ_up(A, p, q, theta);
					timesJ(V, p, q, theta);
				}
			}
		}
	} while (changed != 0);
}
