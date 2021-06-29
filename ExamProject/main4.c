#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_integration.h>

#include "integrator4.h"

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
/*
double fun2GSL(double x, void* params){
	calls++;
	double alpha = *(double*) params;
	return alpha*4*sqrt(1 - pow(x, 2));
}
*/

double fun3(double x){ // fun3(x) = exp(-x)
	calls++;
	return exp(-x);
}
double fun3GSL(double x, void* params){
	calls++;
	double alpha = *(double*) params;
	return alpha*exp(-x);
}

double fun4(double x){ // fun4(x) = exp(x)
	calls++;
	return exp(x);
}
double fun4GSL(double x, void* params){
	calls++;
	double alpha = *(double*) params;
	return alpha*exp(x);
}

double fun5(double x){ // fun5(x) = 1/(1+x^2)
	calls++;
	return 1./(1 + pow(x, 2));
}
double fun5GSL(double x, void* params){
	calls++;
	double alpha = *(double*) params;
	return alpha*1./(1 + pow(x, 2));
}

void testAndPrintIntegrate(double f(double), char* functag, double a, double b, double delta, double epsilon, double trueVal, char* integrationType, int manySignificantDigits){
	// Initialize the call count variable
	calls = 0;
	// Perform integration
	double Q;
	Q = generalisedIntegrator(f, a, b, delta, epsilon);
	// Calculate the estimated and actual error
	double estimatedError = delta + fabs(Q)*epsilon;
	double actualError = fabs(Q - trueVal);
	// Print the found values
	printf("----------\n");
	printf("The integration of %s from %g to %g using %s:\n", functag, a, b, integrationType);
	printf("\t delta = %g, and epsilon = %g.\n", delta, epsilon);
	if (manySignificantDigits) {// 1 = true, 0 = false
		printf("\t Found value     = %.25lg.\n", Q);
	} else {
		printf("\t Found value     = %g.\n", Q);
	}
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
	double delta = 1e-4, epsilon = 1e-4;
	
	// EXERCISE A
	printExercise("A");
	printSubtext("Cauclating som integrals with finite limits using own implementation");
	// fun1
	testAndPrintIntegrate(fun1, "√(x)", 0, 1, delta, epsilon, 2./3, "adaptive and recursive integration with bi-division", 0);
	// fun2
	testAndPrintIntegrate(fun2, "4√(1-x²)", 0, 1, delta, epsilon, M_PI, "adaptive and recursive integration with bi-division", 0);
	
	//EXERCISE C
	printf("\n");
	printExercise("C");
	printSubtext("Calculating (converging) integrals with infinite limits using own implementation");
	// fun3
	printSubtext("Calculating integral of exp(-x) from 0 to infinity using own and GSL rutines");
	testAndPrintIntegrate(fun3, "exp(-x)", 0, INFINITY, delta, epsilon, 1, "adaptive and recursive integration with bi-division", 0);
	// fun4
	printSubtext("Calculating integral of exp(x) from -infinity to 0 using own and GSL rutines");
	testAndPrintIntegrate(fun4, "exp(x)", -INFINITY, 0, delta, epsilon, 1, "adaptive and recursive integration with bi-division", 0);
	// fun5
	printSubtext("Calculating integral of 1/(1+x²) from -infinity to infinity using own and GSL rutines");
	testAndPrintIntegrate(fun5, "1/(1+x²)", -INFINITY, INFINITY, delta, epsilon, M_PI, "adaptive and recursive integratioon with bi-division", 0);
	
	return 0;
}
