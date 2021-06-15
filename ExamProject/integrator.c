#include <stdio.h>
#include <math.h>
#include <assert.h>

double adaptiveRecursiveIntegrate(double f(double), double a, double b, double delta, double epsilon, double func2, int recursionLimit){
	// Initialization of functions surrounding func2 as sain in README.md (like eq. 51 and 48 in the 'Numerical Integration' PDF)
	double func1 = f(a + (b - a)/6);
	double func3 = f(a + 5*(b - a)/6);
	// Initialization of higher and lower order rules (Q and q respectively)
	double Q = (3*func1 + 2*func2 + 3*func3)*(b - a)/6;
	double q = (func1 + func2 + func3)*(b - a)/3;
	// Initialization of the error and the tolerance
	double error = fabs(Q - q);
	double tolerance = delta + epsilon*fabs(Q);
	if (recursionLimit == 0) {
		fprintf(stderr, "Function 'adaptiveRecursiveIntegrate' have reached the recursion limit.\n");
		return Q;
	}
	if (error < tolerance) {
		return Q;
	} else {
		double rescaledAbsoluteAccuracyGoal = delta/sqrt(3);
		return adaptiveRecursiveIntegrate(f, a, (a + b)/3, rescaledAbsoluteAccuracyGoal, epsilon, func1, recursionLimit - 1)
			+ adaptiveRecursiveIntegrate(f, (a + b)/3, 2*(a + b)/3, rescaledAbsoluteAccuracyGoal, epsilon, func2, recursionLimit - 1)
			+ adaptiveRecursiveIntegrate(f, 2*(a + b)/3, b, rescaledAbsoluteAccuracyGoal, epsilon, func3, recursionLimit - 1);
	}
}

double integrateTridivision(double f(double), double a, double b, double delta, double epsilon){
	// Initialized func2 as said in README.md (like eq. 51 and 48 in the 'Numerical Integration' PDF)
	double func2 = f(a + 3*(b - a)/6);
	// Initialize the limit of recursions to 99 times - the 100th time should resolve in an error
	int recursionLimit = 99;
	// Begin recursion
	return adaptiveRecursiveIntegrate(f, a, b, delta, epsilon, func2, recursionLimit);
}

/*
 * DUE TO WSL NOT BEING ABLE TO WORK WITH NESTED FUNCTIONS AND INSTEAD RESULT IN SEGMENTATION FAULT. THE ACTUAL NESTED CODE IS COMMENTED OUT, SUCH THAT IT MIGHT STILL BE TESTED ON ANOTHER SYSTEM, BUT DUE TO THIS, THE FUNCTIONS ARE MOVED OUT HERE.
 */

/*
static double A, B; // Not accesible from other files, and only here due to moved nested functions needing a and b

// If both a and b are INFINITY: Using converted integral from table in "Numerical Integration" PDF (eq. 58)
double gBothLimitsInf(double f(double), double t){
	return ( f((1 - t)/t) + f(- (1 - t)/t ))/pow(t, 2);
}
// If only a is INFINITY: Using converted integral from table in "Numerical Integration" PDF (eq. 62)
double gOnlyLowerLimitInf(double f(double), double t){
	return ( f(B - (1 - t)/t) )/pow(t, 2);
}
// If only b is INFINITY: Using converted integral from table in "Numerical Integration" PDF (eq. 60)
double gOnlyUpperLimitInf(double f(double), double t){
	return ( f(A + (1 - t)/t) )/pow(t, 2);
}
*/

/*
 * END OF FUNCTIONS MOVED FROM NESTED DUE TO WSL.
 */

double generalisedIntegrator(double f(double), double a, double b, double delta, double epsilon){
	//A = a; B = b; // Instead of nested functions
	// Check if a and/or b is INFINITY (isinf returns 0 for not infinity and nonzero for argument being infinity)
	if (isinf(a) != 0) {
		if (isinf(b) != 0) {
			// If both a and b are INFINITY: Using converted integral from table in notes (eq. 58)
			double g(double t){return ( f((1 - t)/t) + f(- (1 - t)/t ))/pow(t, 2);} // Nested function
			return integrateTridivision(g, 0, 1, delta, epsilon);
		} else {
			// If only a is INFINITY: Using converted integral from table in notes (eq. 62)
			double g(double t){return ( f(b - (1 - t)/t) )/pow(t, 2);} // Nested function
			return integrateTridivision(g, 0, 1, delta, epsilon);
		}
	} else if (isinf(b) != 0) {
		// If only b is INFINITY: Using converted integral from table in notes (eq. 60)
		double g(double t){return ( f(a + (1 - t)/t) )/pow(t, 2);} // Nested function
		return integrateTridivision(g, 0, 1, delta, epsilon);
	} else {
		// If neither a nor b are INFINITY: Regular integrate function.
		return integrateTridivision(f, a, b, delta, epsilon);
	}
}
