#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "integratorBidivision.h"

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
 * f3: Double containing the evaluated function for a specified x value (other than the above).
 * recursionLimit: Integer determining the maximum number of recusions before exiting the program with an error.
 * variableTransformationFormula: Integer determining the variable transformation formula to use. 0 if no integration limit is infinity, 1 if lower integration limit is infinity, 2 if upper integration limit is infinity, and 3 if both lower and upper integration limit is infinity.
 * 
 * returns a double containing the recursively and adaptively integrated value of the function in the integration limit.
 */
double recursiveIntegrateBidivision( double f(double),double a, double b, double delta, double epsilon, double f2, double f3, int nrec, int variableTransformationFormula){
	// Initialization of the upper and lower x-values used in the surrounding functions below
	double xValLowerLimit = a + (b - a)/6,
	       xValUpperLimit = a + 5*(b - a)/6;
	// Initialization of functions surrounding func2 and func3 using 'Numerical Integration' PDF eq. 51 and 48
	double f1, f4;
	switch (variableTransformationFormula) {
		case 1:
			// Lower limit infinity transformation
			f1 = gOnlyLowerLimitInf(f, xValLowerLimit);
			f4 = gOnlyLowerLimitInf(f, xValUpperLimit);
			break;
		case 2:
			// Upper limit infinity transformation
			f1 = gOnlyUpperLimitInf(f, xValLowerLimit);
			f4 = gOnlyUpperLimitInf(f, xValUpperLimit);
			break;
		case 3:
			// Both limits infinity transformation
			f1 = gBothLimitsInf(f, xValLowerLimit);
			f4 = gBothLimitsInf(f, xValUpperLimit);
			break;
		default:
			// No variable transformation since variableTransformationFormula == 0
			f1 = f(xValLowerLimit);
			f4 = f(xValUpperLimit);
	}
	// Initialization of higher and lower order rules (Q and q respectively)
	double Q = (2*f1 + f2 + f3 + 2*f4) / 6 * (b - a);
	double q = (f1 + f2 + f3 + f4) / 4 * (b - a);
	// Initialization of the error and the tolerance
	double tolerance = delta + epsilon*fabs(Q);
	double error = fabs(Q - q);
	// Assert the number of recursion
	assert(nrec < 99);
	// Return value: Depending on the found error vs. the tolerance keep doing recursions or return the found integral value
	if (error < tolerance) {
		return Q;
	} else {
		double Q1 = recursiveIntegrateBidivision(f, a, (a+b)/2, delta/sqrt(2), epsilon, f1, f2, nrec + 1, variableTransformationFormula);
		double Q2 = recursiveIntegrateBidivision(f, (a+b)/2, b, delta/sqrt(2), epsilon, f3, f4, nrec + 1, variableTransformationFormula);
		return Q1 + Q2;
	}
}

/*
 * Integrates a function with limits possibly being infinity, using two sub-divisions of the integration interval, and an adaptive and recursive technique.
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
double integrateBidivision(double f(double), double a, double b, double delta, double epsilon, int variableTransformationFormula){
	// Initialize field-scope variables with the integration limits
	A=a; B=b;
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
			break;
		case 4:
			a = 0; b = M_PI;
	}
	// Initialized func2 and func3 from eq. 51 and 48 in the 'Numerical Integration' PDF
	double xValLowerLimit = a + 2*(b - a)/6,
	       xValUpperLimit = a + 4*(b - a)/6;
	double f2, f3;
	switch (variableTransformationFormula) {
		case 1:
			// Lower limit infinity transformation
			f2 = gOnlyLowerLimitInf(f, xValLowerLimit);
			f3 = gOnlyLowerLimitInf(f, xValUpperLimit);
			break;
		case 2:
			// Upper limit infinity transformation
			f2 = gOnlyUpperLimitInf(f, xValLowerLimit);
			f3 = gOnlyUpperLimitInf(f, xValUpperLimit);
			break;
		case 3:
			// Both limits infinity transformation
			f2 = gBothLimitsInf(f, xValLowerLimit);
			f3 = gBothLimitsInf(f, xValUpperLimit);
			break;
		default:
			// No variable transformation since variableTransformationFormula == 0
			f2 = f(xValLowerLimit);
			f3 = f(xValUpperLimit);
	}
	// Initialize the number of recursions
	int nrec = 0;
	// Begin recursion
	return recursiveIntegrateBidivision(f, a, b, delta, epsilon, f2, f3, nrec, variableTransformationFormula);
}

/*
 * Integrates a function with one or more limits possibly being infinity, using two sub-divisions of the integration interval, and an adaptive and recursive technique.
 * 
 * f: Function to be integrated. Takes a double as input and returns a double (the function value at the point).
 * a: Double containing the lower integration limit.
 * b: Double containing the upper integration limit.
 * delta: Double containing the absolute accuracy goal.
 * epsilon: Double containing the relative accuracy goal.
 * 
 * returns a double containing the recursively and adaptively integrated value of the function in the integration limit.
 */
double generalisedIntegratorBidivision(double f(double), double a, double b, double delta, double epsilon){
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
	return integrateBidivision(f, a, b, delta, epsilon, variableTransformationFormula);
}
