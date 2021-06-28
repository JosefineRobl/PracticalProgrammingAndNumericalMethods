#include <stdio.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "quasiNewton.h"
#include "simplex.h"

/*
 * For exercise A.
 */
double Rosenbrock(gsl_vector* r){
	double x = gsl_vector_get(r, 0);
	double y = gsl_vector_get(r, 1);
	return pow(1 - x, 2) + 100*pow(y - pow(x, 2), 2);
}

/*
 * For exercise C. Same as for in exercise A, but without using GSL vectors.
 */
double RosenbrockC(double* r){
	return pow(1 - r[0], 2) + 100*pow(r[1] - pow(r[0], 2), 2);
}

/*
 * For exercise A.
 */
double Himmelblau(gsl_vector* r){
	double x = gsl_vector_get(r, 0);
	double y = gsl_vector_get(r, 1);
	return pow((pow(x, 2) + y - 11), 2) + pow(x + pow(y, 2) - 7, 2);
}

/*
 * For exercise C. Same as for in exercise A, but without using GSL vectors.
 */
double HimmelblauC(double* r){
	return pow(pow(r[0], 2) + r[1] - 11, 2) + pow(r[0] + pow(r[1], 2) - 7, 2);
}

/*
 * Breit-Wigner function. For exercise B.
 */
double BreitWigner(double energy, gsl_vector* p){ 
	double m = gsl_vector_get(p, 0);
	double Gamma = gsl_vector_get(p, 1);
	double A = gsl_vector_get(p, 2);
	return A/(pow(energy - m, 2) + pow(Gamma, 2)/4);
}

/*
 * Deviation function for exercise B.
 */
static int data;
static gsl_vector* E;
static gsl_vector* sigma;
static gsl_vector* dsigma;
double deviation(gsl_vector* p){
	double sum =0;
	for(int i = 0; i < data; i++){
		double Ei = gsl_vector_get(E, i);
		double sigmai = gsl_vector_get(sigma, i); // Function value for i'th point
		double dsigmai = gsl_vector_get(dsigma, i); // Derivative of function with value for i'th point
		sum += pow(BreitWigner(Ei, p) - sigmai, 2)/pow(dsigmai, 2);
	}
	return sum;
}

/*
 * Exercise A.
 */
void exerciseA(void){
	// Open file for writing
	FILE* ExcA = fopen("exerciseA.txt", "w");
	
	// Minimum of Rosenbrock
	int n = 2;
	gsl_vector* x0 = gsl_vector_alloc(n);
	gsl_vector* x = gsl_vector_alloc(n);
	// The initial guess
	gsl_vector_set(x0, 0, -2);
	gsl_vector_set(x0, 1, 8);
	// Accuracy
	double acc = 1e-4;
	// Copy the initial guess for one of the vectors to be used for update to found minimum
	gsl_vector_memcpy(x,x0);
	// Perform
	int steps = quasiNewton(Rosenbrock, x, acc);
	fprintf(ExcA, "=============== Rosenbrock's valley function: ===============\n");
	fprintf(ExcA, "Tolerance       = %g\n", acc);
	fprintf(ExcA, "Initial guess   = (%g,%g)\n", gsl_vector_get(x0, 0), gsl_vector_get(x0, 1));
	fprintf(ExcA, "Found minimum   = (%g,%g)\n", gsl_vector_get(x, 0), gsl_vector_get(x, 1));
	fprintf(ExcA, "Value at min    = %g\n", Rosenbrock(x));
	fprintf(ExcA, "Number of steps = %d\n", steps);
	
	// Set steps = 0 for use in new function also
	steps = 0;
	
	// Minimum of Himmelblau
	gsl_vector* y0 = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_alloc(n);
	// Initial guess on minimum
	gsl_vector_set(y0, 0, 4);
	gsl_vector_set(y0, 1, 2);
	// Accuracy
	double eps = 1e-5;
	// Copy the initial guess for one of the vectors to be used for update to found minimum
	gsl_vector_memcpy(y, y0);
	// Perform
	int steps = quasiNewton(Himmelblau, y, eps);
	fprintf(ExcA, "\n=============== Himmelblau's function: ===============\n");
	fprintf(ExcA, "Tolerance       = %g\n", eps);
	fprintf(ExcA, "Initial guess   = (%g,%g)\n", gsl_vector_get(y0, 0), gsl_vector_get(y0, 1) );
	fprintf(ExcA, "Found minimum   = (%g,%g)\n", gsl_vector_get(y, 0), gsl_vector_get(y, 1) );
	fprintf(ExcA, "Value at min    = %g\n", Himmelblau(y));
	fprintf(ExcA, "Number of steps = %i\n", steps);
	
	// Free GSL vectors
	gsl_vector_free(x0);
	gsl_vector_free(x);
	gsl_vector_free(y0);
	gsl_vector_free(y);
	// Close file for writing
	fclose(ExcA);
}

/*
 * Exercise B.
 */
void exerciseB(void){
	FILE* ExcB = fopen("exerciseB.txt", "w");
	
	// Getting the Higgs data
	data = 30;
	E = gsl_vector_alloc(data);
	sigma = gsl_vector_alloc(data);
	dsigma = gsl_vector_alloc(data);
	// Open files for Higss data
	FILE* energyFile = fopen("higgsEnergy.txt", "r");
	FILE* sigmaFile = fopen("higgsSigma.txt", "r");
	FILE* dsigmaFile = fopen("higgsDsigma.txt", "r");
	// Read files for data 
	gsl_vector_fscanf(energyFile, E);
	gsl_vector_fscanf(sigmaFile, sigma);
	gsl_vector_fscanf(dsigmaFile, dsigma);
	// Close Higgs data files
	fclose(energyFile);
	fclose(sigmaFile);
	fclose(dsigmaFile);
	
	// Number of parameters (m, gamma, A)
	int N = 3;
	// Initialize vectors for the parameters (guess and found values)
	gsl_vector* param0 = gsl_vector_alloc(N);
	gsl_vector* param = gsl_vector_alloc(N);	
	// Initial guess
	gsl_vector_set(param0, 0, 125);
	gsl_vector_set(param0, 1, 1);
	gsl_vector_set(param0, 2, 1);
	// Copy the initial guess for one of the vectors to be used for update to found minimum
	gsl_vector_memcpy(param, param0);
	// Initialize a tolerance
	double tolerance = 1e-3;
	// Perform
	int nsteps = quasiNewton(D, param, tolerance);
	fprintf(ExcB, "=============== Higgs Boson: ===============\n");
	fprintf(ExcB, "Tolerance                     = %g\n", tolerance);
	fprintf(ExcB, "Initial guess (m, Gamma, A)   = (%g,%g,%g)\n", gsl_vector_get(param0, 0), gsl_vector_get(param0, 1), gsl_vector_get(param0, 2));
	fprintf(ExcB, "Estimate values (m, Gamma, A) = (%g,%g,%g)\n", gsl_vector_get(param, 0), gsl_vector_get(param, 1), gsl_vector_get(param, 2));
	fprintf(ExcB, "Number of steps               = %i\n", nsteps);
	
	// Free GSL vectors
	gsl_vector_free(param0);
	gsl_vector_free(param);
	// Close file for writing
	fclose(ExcB);
}

/*
 * Exercise C.
 */
void exerciseC(void){
	// Open file for writing
	FILE* ExcC = fopen("exerciseC.txt", "w");
	
	//Minimum of Rosenbrock's valley function
	int dim = 2;
	// Initial value of simplex. How hard it is for convergence depends on this value
	double* simplex[] = {
		(double[]) {1.0, 2.0},
		(double[]) {0.0, 2.0},
		(double[]) {8.0, 0.0}
	};
	// Initialize size goal for simplex
	double sizeGoal = 1e-4;
	// Initialize index for the lowest value
	int lowValIndex = 0;
	// Perform
	int steps = amoeba(dim, RosenbrockC, simplex, sizeGoal);
	fprintf(ExcC, "=============== Rosenbrock's valley function: ===============\n");
	fprintf(ExcC, "Simplex size goal = %g\n", sizeGoal);
	fprintf(ExcC, "Found minimum at  = (%g, %g)\n", simplex[lowValIndex][0], simplex[lowValIndex][1]);
	fprintf(ExcC, "Number of steps   = %i\n", steps);
	
	// Set steps = 0 for use in new function also
	steps = 0;
	
	// Initial value of simplex
	double* Simplex[] = {
		(double[]) {1.0, 2.0},
		(double[]) {0.0, 2.0},
		(double[]) {8.0, 0.0}
	};
	// Initialize size goal for simplex
	sizeGoal = 1e-4;
	// Initialize index for the lowest value
	lowValIndex = 0;
	// Perform
	int steps = amoeba(dim, HimmelblauC, Simplex, sizeGoal);
	fprintf(ExcC, "=============== Himmelblau's function: ===============\n");
	fprintf(ExcC, "Simplex size goal = %g\n", sizeGoal);
	fprintf(ExcC, "Found minimum at  = (%g,%g)\n", Simplex[lowValIndex][0], Simplex[lowValIndex][1]);
	fprintf(ExcC, "Number of steps   = %i\n", steps);
	
	// Close file for writing
	fclose(ExcC);
}

/*
 * Main function.
 */
int main(void){
	exerciseA();
	exerciseB();
	exerciseC();
	return 0;
}
