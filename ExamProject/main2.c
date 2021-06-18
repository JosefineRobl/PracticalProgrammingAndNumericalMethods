#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "integrator2.h"

int calls;
double f(double x){
	calls++;
	return sqrt(x);
}
double h(double x){
	calls++;
	return 4*sqrt(1 - pow(x, 2));
}

void printIntegrate(double f(double), char* functag, double a, double b, double delta, double epsilon, double trueVal){
	// Initialize the call count variable
	int calls = 0;
	// Perform integration
	double Q = integrate(f, a, b, delta, epsilon);
	// Calculate the estimated and actual error
	double estimatedError = delta + fabs(Q)*epsilon;
	double actualError = fabs(Q - trueVal);
	// Print results
	printf("====================\n");
	printf("The integration of %s from %g to %g using recursive adaptive integration with tri-division:\n", functag, a, b);
	printf("\t delta = %g, and epsilon = %g.\n", delta, epsilon);
	printf("\t Found value     = %g.\n", Q);
	printf("\t Exact value     = %g.\n", trueVal);
	printf("\t Estimated error = %g.\n", estimatedError);
	printf("\t Actual error    = %g.\n", actualError);
	printf("\t Number of calls = %d.\n", calls);
	// Print if test passed, i.e. if actual calculated error is less than the estimated error
	if (estimatedError > actualError) {
		printf("Test passed: Estimated error is larger than the actual error.\n");
	} else {
		printf("Test failed: Actual error is larger than the estimated error.\n");
	}
}

int main(){
	double delta = 1e-6, epsilon = 1e-6;
	printIntegrate(f, "sqrt(t)", 0, 1, delta, epsilon, 2./3);
	printIntegrate(h, "4*sqrt(1-x^2)", 0, 1, delta, epsilon, M_PI);
	/*
	double delta = 0.001, epsilon = 0.001;
	calls = 0;
	double exact = 2./3;
	// sqrtXIntegral = THE BELOW
	double Q = integrate(f, 0, 1, delta, epsilon);
	printf("----------\n");
	printf("∫01 dx √(x)\n");
	printf("              Q = %g\n", Q);
	printf("          exact = %g\n", exact);
	printf("          calls = %d\n", calls);
	printf("estimated error = %g\n", delta + fabs(Q)*epsilon);
	printf("   actual error = %g\n", fabs(Q - exact));
	//
	calls = 0;
	exact = M_PI;
	// fourTimessqrt1MinusXSquaredIntegral = THE BELOW
	Q = integrate(h, 0, 1, delta, epsilon);
	printf("----------\n");
	printf("∫01 dx 4√(1-x²) = π\n");
	printf("              Q = %g\n", Q);
	printf("          exact = %g\n", exact);
	printf("          calls = %d\n", calls);
	printf("estimated error = %g\n", delta + fabs(Q)*epsilon);
	printf("   actual error = %g\n", fabs(Q - exact));
	*/
	//
	return 0;
}
