#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "integratorTridivision.h"

/*
 * DUE TO WSL NOT BEING ABLE TO WORK WITH NESTED FUNCTIONS AND INSTEAD RESULT IN SEGMENTATION FAULT. THE ACTUAL NESTED CODE IS COMMENTED OUT, SUCH THAT IT MIGHT STILL BE TESTED ON ANOTHER SYSTEM, BUT DUE TO THIS, THE FUNCTIONS ARE MOVED OUT HERE.
 */

// Not accesible from other files, and only here due to moved nested functions needing a and b
static double A, B;

// If only a is INFINITY: Using converted integral from table in "Numerical Integration" PDF (eq. 62)
static double gOnlyLowerLimitInf(double f(double), double t){
	return f( B - (1 - t)/t ) / pow(t, 2);
}

// If only b is INFINITY: Using converted integral from table in "Numerical Integration" PDF (eq. 60)
static double gOnlyUpperLimitInf(double f(double),double t){// variable transformation formula
	return f( A + (1 - t)/t ) / pow(t, 2);
}

// If both a and b are INFINITY: Using converted integral from table in "Numerical Integration" PDF (eq. 58)
static double gBothLimitsInf(double f(double), double t){
	return ( f((1-t)/t) + f(-(1-t)/t) ) / pow(t, 2);
}

/*
 * END OF FUNCTIONS MOVED FROM NESTED DUE TO WSL.
 */

/*
 * Adaptively and recursively integrates a function.
 * 
 * f: Function to be integrated. Takes a double as input and returns a double (the function value at that point).
 * a: Double containing the lower integration limit.
 * b: Double containing the upper integration limit.
 * delta: Double containing the absolute accuracy goal.
 * epsilon: Double containing the relative accuracy goal.
 * f2: Double containing the evaluated function for a specified x value.
 * recursionLimit: Integer determining the maximum number of recusions before exiting the program with an error.
 * variableTransformationFormula: Integer determining the variable transformation formula to use. 0 if no integration limit is infinity, 1 if lower integration limit is infinity, 2 if upper integration limit is infinity, and 3 if both lower and upper integration limit is infinity.
 * 
 * returns a double containing the recursively and adaptively integrated value of the function in the integration limit.
 */
double recursiveIntegrateTridivision(double f(double), double a, double b, double delta, double epsilon, double f2, int nrec, int variableTransformationFormula){
	// Initialization of the upper and lower x-values used in the surrounding functions below
	double xValLowerLimit = a + (b - a)/6,
	       xValUpperLimit = a + 5*(b - a)/6;
	// Initialization of functions surrounding f2 (like in 'Numerical Integration' PDF eq. 51)
	double f1, f3;
	switch (variableTransformationFormula) {
		case 1:
			// Lower limit infinity transformation
			f1 = gOnlyLowerLimitInf(f, xValLowerLimit);
			f3 = gOnlyLowerLimitInf(f, xValUpperLimit);
			break;
		case 2:
			// Upper limit infinity transformation
			f1 = gOnlyUpperLimitInf(f, xValLowerLimit);
			f3 = gOnlyUpperLimitInf(f, xValUpperLimit);
			break;
		case 3:
			// Both limits infinity transformation
			f1 = gBothLimitsInf(f, xValLowerLimit);
			f3 = gBothLimitsInf(f, xValUpperLimit);
			break;
		default:
			// No variable transformation since variableTransformationFormula == 0
			f1 = f(xValLowerLimit);
			f3 = f(xValUpperLimit);
	}
	// Initialization of higher and lower order rules (Q and q respectively)
	double Q = (3*f1 + 2*f2 + 3*f3) / 8 * (b - a);
	double q = (f1 + f2 + f3) / 3 * (b - a);
	// Initialization of the error and the tolerance
	double tolerance = delta + epsilon*fabs(Q);
	double error = fabs(Q - q);
	// Assert
	assert(nrec > 0);
	// Return value: Depending on the found error vs. the tolerance keep doing recursions or return the found integral value
	if (error < tolerance) { //((error < tolerance) || (nrec == 0)) {
		return Q;
	} else {
		double deltaNew = delta/sqrt(3);
		double Q1 = recursiveIntegrateTridivision(f, a, a+(a+b)/3, deltaNew, epsilon, f1, nrec - 1, variableTransformationFormula);
		double Q2 = recursiveIntegrateTridivision(f, a+(a+b)/3, a+2*(a+b)/3, deltaNew, epsilon, f2, nrec - 1, variableTransformationFormula);
		double Q3 = recursiveIntegrateTridivision(f, a+2*(a+b)/3, b, deltaNew, epsilon, f3, nrec - 1, variableTransformationFormula);
		return Q1 + Q2 + Q3;
	}
}

/*
 * Integrates a function with limits possibly being infinity, using three sub-divisions of the integration interval instead of two and an adaptive and recursive technique.
 * 
 * f: Function to be integrated. Taks a double as input and returns a double (the function value at that point).
 * a: Double containing the lower integration limit.
 * b: Double containing the upper integration limit.
 * delta: Double containing the absolute accuracy goal.
 * epsilon: Double containing the relative accuracy goal.
 * variableTransformationFormula: Integer determining the variable transformation formula to use. 0 if no integration limit is infinity, 1 if lower integration limit is infinity, 2 if upper integration limit is infinity, and 3 if both lower and upper integration limit is infinity.
 * 
 * returns a double containing the recursively and adaptively integrated value of the function in the integration limit.
 */
double integrateTridivision(double f(double), double a, double b, double delta, double epsilon, int variableTransformationFormula){
	// Initialize field-scope variables with the integration limits
	A = a; B = b;
	// Set new limits for use in the transformation formulas
	switch (variableTransformationFormula) {
		case 1:
			a = 0; b = 1;
			break;
		case 2:
			a = 0; b = 1;
			break;
		case 3:
			a = 0; b = 1;
	}
	// Initialized f2 (as from eq. 51 and 48 in the 'Numerical Integration' PDF)
	double xVal = a + 3*(b - a)/6;
	double f2;
	switch (variableTransformationFormula) {
		case 1:
			// Lower limit infinity transformation
			f2 = gOnlyLowerLimitInf(f, xVal);
			break;
		case 2:
			// Upper limit infinity transformation
			f2 = gOnlyUpperLimitInf(f, xVal);
			break;
		case 3:
			// Both limits infinity transformation
			f2 = gBothLimitsInf(f, xVal);
			break;
		default:
			// No variable transformation since variableTransformationFormula == 0
			f2 = f(xVal);
	}
	// Initialize the number of recursions
	int nrec = (int) 99;
	// Begin recursion
	return recursiveIntegrateTridivision(f, a, b, delta, epsilon, f2, nrec, variableTransformationFormula);
}

/*
 * Integrates a function with one or more limits possibly being infinity, using three sub-divisions of the integration interval instead of two, and an adaptive and recursive technique.
 * 
 * f: Function to be integrated. Takes a double as input and returns a double (the function value at the point).
 * a: Double containing the lower integration limit.
 * b: Double containing the upper integration limit.
 * delta: Double containing the absolute accuracy goal.
 * epsilon: Double containing the relative accuracy goal.
 * 
 * returns a double containing the recursively and adaptively integrated value of the function in the integration limit.
 */
double generalisedIntegratorTridivision(double f(double), double a, double b, double delta, double epsilon){
	int variableTransformationFormula;
	if (isinf(a) != 0) {
		if (isinf(b) != 0) {
			variableTransformationFormula = 3;
		} else {
			variableTransformationFormula = 1;
		}
	} else if (isinf(b) != 0) {
		variableTransformationFormula = 2;
	} else {
		variableTransformationFormula = 0;
	}
	return integrateTridivision(f, a, b, delta, epsilon, variableTransformationFormula);
}
