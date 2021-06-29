#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "integratorTridivision.h"
#include "integratorBidivision.h"

void printExercise(char* exercise){
	printf("==================== Exercise %s ====================\n", exercise);
}

void printSubtext(char* subtext){
	printf("---------- %s ----------\n", subtext);
}

// Define functions

int calls; // Variable for storing the number of calls the integration function makes

double fun1(double x){ // fun1(x) = sqrt(x)
	calls++;
	return sqrt(x);
}

double fun2(double x){ // fun2(x) = 4*sqrt(1-x^2)
	calls++;
	return 4*sqrt(1 - pow(x, 2));
}

double fun3(double x){ // fun3(x) = exp(-x)
	calls++;
	return exp(-x);
}

double fun4(double x){ // fun4(x) = exp(x)
	calls++;
	return exp(x);
}

double fun5(double x){ // fun5(x) = 1/(1+x^2)
	calls++;
	return 1./(1 + pow(x, 2));
}

void testAndPrintIntegrate(double f(double), char* functag, double a, double b, double delta, double epsilon, double trueVal, int integrator){
	assert((integrator == 1) || (integrator == 2));
	// Initialize the call count variable
	calls = 0;
	// Perform integration
	double Q;
	char* integrationType;
	switch (integrator) {
		case 1:
			Q = generalisedIntegratorTridivision(f, a, b, delta, epsilon);
			integrationType = "adaptive and recursive integration with bi-division";
		case 2:
			Q = generalisedIntegratorBidivision(f, a, b, delta, epsilon);
			integrationType = "adaptive and recursive integration with tri-division";
	}
	// Calculate the estimated and actual error
	double estimatedError = delta + fabs(Q)*epsilon;
	double actualError = fabs(Q - trueVal);
	// Print the found values
	printf("The integration of %s from %g to %g using %s:\n", functag, a, b, integrationType);
	printf("\t delta = %g, and epsilon = %g.\n", delta, epsilon);
	printf("\t Found value     = %g.\n", Q);
	printf("\t Exact value     = %g.\n", trueVal);
	printf("\t Estimated error = %g.\n", estimatedError);
	printf("\t Actual error    = %g.\n", actualError);
	printf("\t Number of calls = %d.\n", calls);
	if (estimatedError > actualError) {
		printf("Test passed: Estimated error is larger than the actual error.\n");
	} else {
		printf("Test failed: Actual error is larger than the estimated error.\n");
	}
}

int main(void){
	printf("\n\nMy student number is 201706760, and mod(60, 22) = 16.\n");

	printf("\n###################################");
	printf("\n# ------------------------------- #");
	printf("\n# ------- EXAM PROJECT 16 ------- #");
	printf("\n# ------------------------------- #");
	printf("\n###################################\n");
	printf("\nTitle: Adaptive recursive integrator with subdivision into three subintervals.\n");
	printf("\nDescription:\n");
	printf("Implement a (one-dimensional) adaptive recursive integrator (open or closed quadrature, at your choice) which at each iteration subdivides the interval not into two, but into three sub-intervals. Reuse points. Compare with the adaptive integrator from your homework.\n\n");

	// Absolute and relative accuracy goals
	double delta = 1e-4, epsilon = 1e-4;

	// fun1
	printExercise("Cauclating integral of √(x) from 0 to 1");
	testAndPrintIntegrate(fun1, "√(x)", 0, 1, delta, epsilon, 2./3, 1);
	testAndPrintIntegrate(fun1, "√(x)", 0, 1, delta, epsilon, 2./3, 2);
	// fun2
	printExercise("Cauclating integral of 4√(1-x²) from 0 to 1");
	testAndPrintIntegrate(fun2, "4√(1-x²)", 0, 1, delta, epsilon, M_PI, 1);
	testAndPrintIntegrate(fun2, "4√(1-x²)", 0, 1, delta, epsilon, M_PI, 2);
	// fun3
	printExercise("Calculating integral of exp(-x) from 0 to infinity");
	testAndPrintIntegrate(fun3, "exp(-x)", 0, INFINITY, delta, epsilon, 1, 1);
	testAndPrintIntegrate(fun3, "exp(-x)", 0, INFINITY, delta, epsilon, 1, 2);
	// fun4
	printExercise("Calculating integral of exp(x) from -infinity to 0");
	testAndPrintIntegrate(fun4, "exp(x)", -INFINITY, 0, delta, epsilon, 1, 1);
	testAndPrintIntegrate(fun4, "exp(x)", -INFINITY, 0, delta, epsilon, 1, 2);
	// fun5
	printExercise("Calculating integral of 1/(1+x²) from -infinity to infinity");
	testAndPrintIntegrate(fun5, "1/(1+x²)", -INFINITY, INFINITY, delta, epsilon, M_PI, 1);
	testAndPrintIntegrate(fun5, "1/(1+x²)", -INFINITY, INFINITY, delta, epsilon, M_PI, 2);

	printf("========================================\n");
	printf("\nConclusion:\n");
	printf("\nThe integration using tri-division uses less recursive calls in general compared to the integration with bi-division, thus it performs better overall.\n");
	
	return 0;
}
