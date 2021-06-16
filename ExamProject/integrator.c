#include <stdio.h>
#include <math.h>

/*
 * DUE TO WSL NOT BEING ABLE TO WORK WITH NESTED FUNCTIONS AND INSTEAD RESULT IN SEGMENTATION FAULT. THE ACTUAL NESTED CODE IS COMMENTED OUT, SUCH THAT IT MIGHT STILL BE TESTED ON ANOTHER SYSTEM, BUT DUE TO THIS, THE FUNCTIONS ARE MOVED OUT HERE.
 */

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
 * func2: Double containing the evaluated function for a specified x value.
 * recursionLimit: Integer determining the maximum number of recusions before exiting the program with an error.
 * variableTransformationFormula: Integer determining the variable transformation formula to use. 0 if no integration limit is infinity, 1 if lower integration limit is infinity, 2 if upper integration limit is infinity, and 3 if both lower and upper integration limit is infinity.
 * 
 * returns a double containing the recursively and adaptively integrated value of the function in the integration limit.
 */
double adaptiveRecursiveIntegrate(double f(double), double a, double b, double delta, double epsilon, double func2, int recursionLimit, int variableTransformationFormula){
	// Initialize field-scope variables with the integration limits
	A = a; B = b;
	// Initialization of functions surrounding func2 as sain in README.md (like eq. 51 and 48 in the 'Numerical Integration' PDF)
	double xValLower = a + (b - a)/6,
	       xValUpper = a + 5*(b - a)/6;
	double func1, func3;
	/*
	if (variableTransformationFormula == 0) {
		// No variable transformation
		func1 = f(xValLower);
		func3 = f(xValUpper);
	}
	*/
	if (variableTransformationFormula == 1) {
		// Lower limit infinity transformation
		func1 = gOnlyLowerLimitInf(f, xValLower);
		func3 = gOnlyLowerLimitInf(f, xValUpper);
	} else if (variableTransformationFormula == 2) {
		// Upper limit infinity transformation
		func1 = gOnlyUpperLimitInf(f, xValLower);
		func3 = gOnlyUpperLimitInf(f, xValUpper);
	} else if (variableTransformationFormula == 3) {
		// Both limits infinity transformation
		func1 = gBothLimitsInf(f, xValLower);
		func3 = gBothLimitsInf(f, xValUpper);
	} else {
		// No variable transformation since variableTransformationFormula == 0.
		func1 = f(xValLower);
		func3 = f(xValUpper);
	}
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
		return adaptiveRecursiveIntegrate(f, a, (a + b)/3, rescaledAbsoluteAccuracyGoal, epsilon, func1, recursionLimit - 1, variableTransformationFormula)
			+ adaptiveRecursiveIntegrate(f, (a + b)/3, 2*(a + b)/3, rescaledAbsoluteAccuracyGoal, epsilon, func2, recursionLimit - 1, variableTransformationFormula)
			+ adaptiveRecursiveIntegrate(f, 2*(a + b)/3, b, rescaledAbsoluteAccuracyGoal, epsilon, func3, recursionLimit - 1, variableTransformationFormula);
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
	// Initialized func2 as said in README.md (like eq. 51 and 48 in the 'Numerical Integration' PDF)
	double xVal = a + 3*(b - a)/6;
	double func2;
	/*
	if (variableTransformationFormula == 0) {
		// No variable transformation
		func2 = f(xVal);
	}
	*/
	if (variableTransformationFormula == 1) {
		// Lower limit infinity transformation
		func2 = gOnlyLowerLimitInf(f, xVal);
	} else if (variableTransformationFormula == 2) {
		// Upper limit infinity transformation
		func2 = gOnlyUpperLimitInf(f, xVal);
	} else if (variableTransformationFormula == 3) {
		// Both limits infinity transformation
		func2 = gBothLimitsInf(f, xVal);
	} else {
		// No variable transformation since variableTransformationFormula == 0
		func2 = f(xVal);
	}
	/*
	else {
		fprintf(stderr, "Error in integratTridivision: The argument passed for the type of variable transformation is not an integer in [0, 3].\n");
	}
	*/
	// Initialize the limit of recursions to 99 times - the 100th time should resolve in an error
	int recursionLimit = 99;
	// Begin recursion
	return adaptiveRecursiveIntegrate(f, a, b, delta, epsilon, func2, recursionLimit, variableTransformationFormula);
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
double generalisedIntegrator(double f(double), double a, double b, double delta, double epsilon){
	// Initializes the variable for knowing which variable function formula to use
	int variableTransformationFormula;
	// Check if a and/or b is INFINITY (isinf returns 0 for not infinity and nonzero for argument being infinity)
	if (isinf(a) != 0) {
		if (isinf(b) != 0) {
			// If both a and b are INFINITY: Using converted integral from table in notes (eq. 58)
			variableTransformationFormula = 3;
			a = 0; b = 1;
			//double g(double t){return ( f((1 - t)/t) + f(- (1 - t)/t ))/pow(t, 2);} // Nested function
			//return integrateTridivision(f, 0, 1, delta, epsilon); // When nested function
		} else {
			// If only a is INFINITY: Using converted integral from table in notes (eq. 62)
			variableTransformationFormula = 1;
			a = 0; b = 1;
			//double g(double t){return ( f(b - (1 - t)/t) )/pow(t, 2);} // Nested function
			//return integrateTridivision(g, 0, 1, delta, epsilon); // When nested function
		}
	} else if (isinf(b) != 0) {
		// If only b is INFINITY: Using converted integral from table in notes (eq. 60)
		variableTransformationFormula = 2;
		a = 0; b = 1;
		//double g(double t){return ( f(a + (1 - t)/t) )/pow(t, 2);} // Nested function
		//return integrateTridivision(f, 0, 1, delta, epsilon); // When nested function
	} else {
		// If neither a nor b are INFINITY: Regular integrate function.
		variableTransformationFormula = 0;
	}
	// Return value
	return integrateTridivision(f, a, b, delta, epsilon, variableTransformationFormula);
}