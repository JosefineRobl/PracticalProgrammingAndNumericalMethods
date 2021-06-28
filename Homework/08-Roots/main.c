#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<float.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>

#include "qrGramSchmidt.h"
#include "rungeKutta23.h"
#include "rootFinding.h"

int ncalls;

/*
 * Gradient of Rosenbrock's valley function for exercise A.
 */
void RosenbrockGradient(gsl_vector* r, gsl_vector* Rr){
	ncalls++;
	double x = gsl_vector_get(r, 0);
	double y = gsl_vector_get(r, 1);
	double Rx = 2*(x - 1) + 400*(pow(x, 2) - y)*x;
	double Ry = 200*(y - pow(x, 2));
	gsl_vector_set(Rr, 0, Rx);
	gsl_vector_set(Rr, 1, Ry);
}

/*
 * Schrödinger ODE to be solved, exercise B.
 */
static double e; // energy
void schroedingerEquation(int n, double x, double* y, double* dydx) {
	dydx[0] = y[1];
	dydx[1] = -2*e*y[0] - 2./x * y[0];
}

/*
 * Function giving M(e)=F_e(rmax) for exercise B.
 */
static double rMax=8;
char* pathToFile;
void solverForSchrodinger(gsl_vector* x, gsl_vector* M){
	ncalls++;
	// The energy is the variable of M
	E = gsl_vector_get(x,0);
	// Second order ODE
	int n = 2;
	// For avoid to dividing by zero a is non-zero but only close
	double a = 1e-3;
	double b = rMax;
	// Initialize lists for the initial values
	double ya[n];
	double yb[n];
	// Constants to be declared
	double h = 0.001;
	double acc = 1e-4;
	double epsilon = 1e-4;
	// Settign the initial values. These come from the equations f=r-r^2 and dfdr = 1-2r for a small r.
	ya[0] = a - pow(a, 2);
	ya[1] = 1 - 2*a;
	// Solve
	rungeKuttaDrive23(n, schroedingerEquation, a, b, ya, yb, h, acc, epsilon, pathToFile);
	// M is found as the solution evaluated at the endpoint
	gsl_vector_set(M, 0, yb[0]);
}

/*
 * Exercise A.
 */
void exerciseA(void){
	FILE* excA = fopen("exerciseA.txt", "w");
	fprintf(excA, "=============== Extremum of Rosenbrock valley function ===============\n");
	// Initialize accuracy
	double epsilon = 0.01;
	int d = 2;
	gsl_vector* r = gsl_vector_alloc(d);
	// The initial guess for roots
	double x0 = -2;
	double y0 = 8;
	gsl_vector_set(r, 0, x0);
	gsl_vector_set(r, 1, y0);
	// Initialize the calls variable
	ncalls = 0;
	// Calculate
	newton(RosenbrockGradient, r, epsilon);
	// Print the result
	fprintf(excA, "We start searching at (x, y) = (%g, %g)\n", x0, y0);
	fprintf(excA, "The function is called %i times\n", calls);
	fprintf(excA, "The found extremum is:\n");
	for(int i=0; i<d; i++){
	       fprintf(excA, "%g\n", gsl_vector_get(r,i));
	}
	// Close the file for writing
	fclose(Exc_A);
	// Free allocated memory
	gsl_vector_free(r);
}

/*
 * Exercise B.
 */
void exerciseB(void){
	// Open file for writing
	FILE* excB = fopen("exerciseB.txt", "w");
	pathToFile = "boundStatesHydrogen.txt";
	// Initialize accuracy
	double epsilon = 0.01;
	// Allocate space for the solution (starting from the guess)
	int d = 1;
	gsl_vector* x = gsl_vector_alloc(d);
	double x0 = -2; // Initial guess
	gsl_vector_set(x, 0, x0);
	// Set the calls variable to 0
	ncalls = 0;
	// Perform
	newton(solverForSchrodinger, x, epsilon);
	// Print results
	fprintf(excB, "=============== Lowest root of M(E) = 0: ===============\n");
	fprintf(excB, "Chosen rmax           = %g\n", rMax);
	fprintf(excB, "We start searching at = %g\n", x0);
	fprintf(excB, "Lowest found root     = %g\n", gsl_vector_get(x,0));
	fprintf(excB, "Exact result would be = -1/2\n");
	// Close file
	fclose(excB);
	// Free allocated vector
	gsl_vector_free(x);
}

int main(){
	exerciseA();
	exerciseB();
	
	return 0;
}
