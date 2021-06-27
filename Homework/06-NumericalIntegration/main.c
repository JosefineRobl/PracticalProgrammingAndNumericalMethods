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

void testAndPrintIntegrateGSL(double f(double, void*), char* functag, double a, double b, double delta, double epsilon, double trueVal, int integrator, int manySignificantDigits){
	assert((integrator >= 0) && (integrator <= 3));
	// Initialize call variable
	calls = 0;
	// Initialize variables for GSL integration rutine
	double limit = 10000,
	       alpha = 1.0,
	       result, error;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
	gsl_function G;
	G.function = f;
	G.params = &alpha;
	char* integrationType;
	// Perfom integration
	switch (integrator) {
		case 0:
			gsl_integration_qags(&G, a, b, delta, epsilon, limit, w, &result, &error);
			integrationType = "GSL QAGS rutine";
			break;
		case 1:
			gsl_integration_qagil(&G, b, delta, epsilon, limit, w, &result, &error);
			integrationType = "GSL QAGI rutine (lower limit inf: QAGIL)";
			break;
		case 2:
			gsl_integration_qagiu(&G, a, delta, epsilon, limit, w, &result, &error);
			integrationType = "GSL QAGI rutine (upper limit inf: QAGIU)";
			break;
		case 3:
			gsl_integration_qagi(&G, delta, epsilon, limit, w, &result, &error);
			integrationType = "GSL QAGI rutine (both limits inf: QAGI)";
	}
	// Calculate the estimated error
	double estimatedError = delta + fabs(result)*epsilon,
	       actualError = fabs(result - trueVal);
	// Print the found values
	printf("----------\n");
	printf("The integration of %s from %g to %g using %s:\n", functag, a, b, integrationType);
	if (manySignificantDigits) {// 1 = true, 0 = false
		printf("\t Found value     = %.25lg.\n", result);
	} else {
		printf("\t Found value     = %g.\n", result);
	}
	printf("\t Exact value     = %g.\n", trueVal);
	printf("\t Estimated error = %g.\n", estimatedError);
	printf("\t Actual error    = %g.\n", actualError);
	printf("\t Number of calls = %d.\n", calls);
	// Free the GSL workspace
	gsl_integration_workspace_free(w);
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
	testAndPrintIntegrateGSL(&fun2GSL, "4√(1-x²)", 0, 1, delta, epsilon, M_PI, 0, 1);
	
	//EXERCISE C
	printf("\n");
	printExercise("C");
	printSubtext("Calculating (converging) integrals with infinite limits using own implementation");
	// fun3
	printSubtext("Calculating integral of exp(-x) from 0 to infinity using own and GSL rutines");
	testAndPrintIntegrate(fun3, "exp(-x)", 0, INFINITY, delta, epsilon, 1, "adaptive and recursive integration with bi-division", 1, 0);
	testAndPrintIntegrateGSL(&fun3GSL, "exp(-x)", 0, INFINITY, delta, epsilon, 1, 2, 0);
	// fun4
	printSubtext("Calculating integral of exp(x) from -infinity to 0 using own and GSL rutines");
	testAndPrintIntegrate(fun4, "exp(x)", -INFINITY, 0, delta, epsilon, 1, "adaptive and recursive integration with bi-division", 1, 0);
	testAndPrintIntegrateGSL(&fun4GSL, "exp(-x)", -INFINITY, 0, delta, epsilon, 1, 1, 0);
	// fun5
	printSubtext("Calculating integral of 1/(1+x²) from -infinity to infinity using own and GSL rutines");
	testAndPrintIntegrate(fun5, "1/(1+x²)", -INFINITY, INFINITY, delta, epsilon, M_PI, "adaptive and recursive integratioon with bi-division", 1, 0);
	testAndPrintIntegrateGSL(&fun5GSL, "1/(1+x²)", -INFINITY, INFINITY, delta, epsilon, M_PI, 3, 0);
	
	return 0;
}
