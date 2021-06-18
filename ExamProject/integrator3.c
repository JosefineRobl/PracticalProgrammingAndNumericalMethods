#include <stdio.h>
#include <math.h>
#include <assert.h>

double recursiveAdaptiveIntegrate(double f(double), double a, double b, double delta, double epsilon, double func2, double func3, int recursionLimit){
	//
	double func1 = f(a + (b - a)/6);
	double func4 = f(a + 5*(b - a)/6);
	//
	double Q = (2*func1 + func2 + func3 + 2*func4)*(b - a)/6;
	double q = (func1 + func2 + func3 + func4)*(b - a)/4;
	//
	double error = fabs(Q - q);
	double tolerance = delta + epsilon*fabs(Q);
	
	//assert(recursionLimit > 0);
	if (recursionLimit == 0) {
		fprintf(stderr, "Recursion limit reached\n");
		return Q;
	}
	if (error < tolerance) {
		return Q;
	} else {
		return recursiveAdaptiveIntegrate(f, a, (a + b)/2, delta/sqrt(2), epsilon, func1, func2, recursionLimit - 1)
			+ recursiveAdaptiveIntegrate(f, (a + b)/2, b, delta/sqrt(2), epsilon, func3, func4, recursionLimit - 1);
	}
}

double integrate(double f(double), double a, double b, double delta, double epsilon){
	//
	double func2 = f(a + 2*(b - a)/6);
	double func3 = f(a + 4*(b - a)/6);
	//
	int recursionLimit = 99;
	//
	return recursiveAdaptiveIntegrate(f, a, b, delta, epsilon, func2, func3, recursionLimit);
}
