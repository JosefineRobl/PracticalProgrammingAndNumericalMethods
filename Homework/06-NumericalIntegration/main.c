#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_integration.h>

#include "adaptiveIntegration.h"
#include "adaptiveIntegration2.h"

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
double fun2GSL(double x, void* params){
	calls++;
	double alpha = *(double*) params;
	return alpha*4*sqrt(1 - pow(x, 2));
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
double fun6(double x){ // fun6(x) = 1/sqrt(x)
	calls++;
	return 1./sqrt(x);
}
double fun7(double x){ // fun7(x) = ln(x)/sqrt(x)
	calls++;
	return log(x)/sqrt(x);
}

void testAndPrintIntegrate(double f(double), char* functag, double a, double b, double delta, double epsilon, double trueVal, char* integrationType, int integrator, int manySignificantDigits){
	assert((integrator == 1) || (integrator == 2));
	// Initialize the call count variable
	calls = 0;
	// Perform integration
	double Q;
	switch (integrator) {
		case 1:
			Q = generalisedIntegrator2(f, a, b, delta, epsilon);
			break;
		case 2:
			Q = ClenshawCurtisIntegrate(f, a, b, delta, epsilon);
	}
	// Calculate the estimated and actual error
	double estimatedError = delta + fabs(Q)*epsilon;
	double actualError = fabs(Q - trueVal);
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
	testAndPrintIntegrate(fun1, "√(x)", 0, 1, delta, epsilon, 2./3, "adaptive and recursive integration with bi-division", 1, 0);
	// fun2
	testAndPrintIntegrate(fun2, "4√(1-x²)", 0, 1, delta, epsilon, M_PI, "adaptive and recursive integration with bi-division", 1, 0);
	
	// EXERCISE B
	printf("\n");
	printExercise("B");
	printSubtext("Calculating som integrals with finite limits using Clenshaw-Curtis variable transformation");
	// fun6
	testAndPrintIntegrate(fun6, "1/√(x)", 0, 1, delta, epsilon, 2, "Clenshaw-Curtis variable transformation", 2, 0);
	// fun7
	testAndPrintIntegrate(fun7, "ln(x)/√(x)", 0, 1, delta, epsilon, -4, "Clenshaw-Curtis variable transformation", 2, 0);
	printSubtext("Calculating the integral ∫_0^1 dx 4√(1-x²) = π with normal, Clenshaw-Curtis and GSL integrator");
	// fun2: Normal
	testAndPrintIntegrate(fun2, "4√(1-x²)", 0, 1, delta, epsilon, M_PI, "adaptive and recursive integration with bi-division", 1, 1);
	// fun2: Clenshaw-Curtis
	testAndPrintIntegrate(fun2, "4√(1-x²)", 0, 1, delta, epsilon, M_PI, "Clenshaw-Curtis variable transformation", 2, 1);
	// fun2: GSL
	double limit = 10000,
	       alpha = 1.0,
	       error, result;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
	gsl_function G;
	G.function = &fun2GSL;
	G.params = &alpha;
	gsl_integration_qags(&G, 0, 1, delta, epsilon, limit, w, &result, &error);
	printf("----------\n");
	printf("The integration of 4√(1-x²) from 0 to 1 using GSL QAGS integration rutine:\n");
	printf("\t Found value     = %.25lg.\n", result);
	printf("\t Exact value     = %.25lg.\n", M_PI);
	gsl_integration_workspace_free(w);
	
	//EXERCISE C
	printf("\n");
	printExercise("C");
	// fun3
	printSubtext("Calculating (converging) integrals with infinite limits using own implementation");
	testAndPrintIntegrate(fun3, "exp(-x)", 0, INFINITY, delta, epsilon, 1, "adaptive and recursive integration with bi-division", 1, 0);
	//testAndPrintIntegrate(fun4, "exp(x)", -INFINITY, 0, delta, epsilon, 1, "adaptive and recursive integration with bi-division", 1, 0);
	//testAndPrintIntegrate(fun5, "1/(1+x²)", -INFINITY, INFINITY, delta, epsilon, M_PI, "adaptive and recursive integratioon with bi-division", 1, 0);

	printSubtext("Calculating (converging) integrals with infinite limits using GSL's inegration rutines");
	
	return 0;
}
