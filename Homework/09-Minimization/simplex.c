#include <stdio.h>
#include <math.h>

void simplexContraction(double* highest, double* centroid, double* contracted, int dim){
	for (int i = 0; i < dim; i++) {
		contracted[i] = 1./2 * (centroid[i] + highest[i]);
}

void simplexReflection(double* highest, double* centroid, double* reflected, int dim){
	for (int i = 0; i < dim; i++) {
		reflected[i] = 2*centroid[i] - highest[i];
}

void simplexExpansion(double* highest, double* centroid, double* expanded, int dim){
	for (int i = 0; i < dim; i++) {
		expanded[i] = 3*centroid[i] - 2*highest[i];
}

void simplexReduction(double** simplex, int dim, int lo){
	for (int i = 0; i < dim + 1; i++) {
		if (i != lowValIndex) {
			for(int j = 0; j < dim; j++) {
				simplex[i][j] = 1./2 * (simplex[i][j] + simplex[lowValIndex][j]);
			}
		}
	}
}

double simplexDistance(double* a, double* b, int dim){
	// Initialize variable and update it in for-loop
	double s = 0;
	for (int i = 0; i < dim; i++) {
		s += pow(b[i] - a[i], 2);
	return sqrt(s);
}

/*
 * Calculates the size of simplex
 */
double simplexSize(double** simplex, int dim){
	// Initialize and 
	double s = 0;
	for(int i = 1; i < dim + 1; j++){
		double d = distance(simplex[0], simplex[i], dim);
		if(d > s) {
			s = d;
		}
	}
	return s;
}

/*
 * Updates simplex.
 */
void simplexUpdate(int dim, double** simplex, double* f_val, int* highValIndex, int* lowValIndex, double* centroid){
	*highValIndex = 0;
	*lowValIndex = 0;
	double fHigh = fVal[highValIndex];
	double fLow = fVal[lowValIndex];

	// Find highest and lowest value
	for (int i = 1; i < dim + 1; i++){
		double fx = fVal[i];
		// Update index for high or low val if such is found
		if (fx > fHigh) {
			fHigh = fx;
			*highValIndex = i;
		}
		if (fx < fLow) {
			fLow = fx;
			*lowValIndx = i;
		}
	}
	//Finding centroid
	for (int i = 0; i < dim; i++){
		double sum = 0;
		for (int j = 0; j < dim + 1; j++) {
			if (j != *highValIndex) {
				sum += simplex[j][i];
			}
			centroid[i] = sum/dim;
		}
	}
}
	
/*
 * Initialize simplex.
 */
void simplexInitialize(int dim, double fun(double*), double** simplex, double* fVal, int* highValIndex, int* lowValIndex, double* centroid){
	// Calculate the function value at vertices
	for(int i = 0; i < dim + 1; i++) {
		fVal[i] = fun(simplex[i]);
	}
	simplexUpdate(dim, simplex, fVal, highValIndex, lowValIndex, centroid);
}


int amoeba(int dim, double f(double*), double** simplex, double simplexSizeGoal){
	int steps = 0;
	int highValIndex, lowValIndex;
	double centroid[dim];
	double fVal[dim + 1];
	double p1[dim], p2[dim];
	
	// Initialization of simplex
	simplexInitialize(dim, f, simplex, fVal, &highValIndex, &lowValIndex, centroid);
	
	// Loop until it is broken
	while(1){
		if (simplexSize(simplex, dim) < simplexSizeGoal) {
			break;
		}
		// Perform reflection
		simplexUpdate(dim, simplex, fVal, &highValIndex, &lowValIndex, centroid);
		reflection(simplex[highValIndex], centroid, p1, dim);
		double fReflection = f(p1); // Reflection variable
		if (fReflection < fVal[lowValIndex]) {
			// If the reflection was performed okay try then the expansion
			simplexExpansion(simplex[highValIndex], centroid, p2, dim);
			double fExpansion = f(p2); // Expansion variable
			// Check if expansion okay
			if (fExpansion < fReflection) {
				// Accept of the expansion
				for (int i = 0; i < dim; i++) {
					simplex[highValIndex][i] = p2[i];
					fVal[highValIndex] = fExpansion;
				}
			} else { 
				// Reject expansion and accept reflection
				for (int i = 0; i < dim; i++) {
					simplex[highValIndex][i] = p1[i];
					fVal[highValIndex] = fReflection;
				}
			}
		} else {
			// If the reflection is not good, the following is tried
			if (fReflection < fVal[highValIndex]) {
				// We accept reflection anyway
				for (int i = 0; i < dim; i++) {
					simplex[highValIndex][i] = p1[i];
					fVal[highValIndex] = fReflection;
				}
			} else {
				// If the reflection is not accepted, the contraction is tried
				simplexContraction(simplex[highValIndex], centroid, p1, dim); 
				double fContraction = f(p1); // Contraction variable
				if (fContraction < fVal[highValIndex]){
					// Accepting contraction
					for (int i = 0; i < dim; i++) {
						simplex[highValIndex][i] = p1[i];
						fVal[highValIndex] = fReflection;
					}
				} else {
					// try reduction
					simplexReduction(simplex, dim, lo);
					simplexInitialize(dim, f, simplex, f_val, &hi, &lo, centroid);
				}
			}
		}
		// Updating the steps variable
		steps++;
	}
	return steps;
}
