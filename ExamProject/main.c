#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "integrator.h"

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
/*
double fun3(double x){ // fun3(x) = exp(-x)
	calls++;
	return exp(-x);
}
*/

void printIntegrate(double f(double), char* functag, double a, double b, double delta, double epsilon, double trueVal){
	// Initialize the call count variable
	int calls = 0;
	// Perform integration
	double Q = integrateTridivision(f, a, b, delta, epsilon);
	// Calculate the estimated and actual error
	double estimatedError = delta + fabs(Q)*epsilon;
	double actualError = fabs(Q - trueVal);
	printf("====================\n");
	printf("The integration of %s from %g to %g using recursive adaptive integration with tri-division:\n", functag, a, b);
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
	// Absolute and relative accuracy goals
	double delta = 0.001, epsilon = 0.001;
	
	// fun1
	printIntegrate(fun1, "√(x)", 0, 1, delta, epsilon, 2./3);
	
	// fun2
	printIntegrate(fun2, "4√(1-x²)", 0, 1, delta, epsilon, M_PI);
	
	// fun3
	//printIntegrate(fun3, "exp(-x)", 0, INFINITY, delta, epsilon, 1);
	
	return 0;
}
