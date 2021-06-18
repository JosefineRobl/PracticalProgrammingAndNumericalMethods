#include <stdio.h>
#include <math.h>
#include <assert.h>

double recursiveAdaptiveIntegrate(double f(double), double a, double b, double delta, double epsilon, double func2, int recursionLimit){
	//
	double func1 = f(a + (b - a)/6);
	double func3 = f(a + 5*(b - a)/6);
	//
	double Q = (3*func1 + 2*func2 + 3*func3)*(b - a)/8;
	double q = (func1 + func2 + func3)*(b - a)/3;
	//
	double error = fabs(Q - q);
	double tolerance = delta + epsilon*fabs(Q);
	
	assert(recursionLimit > 0);
	/*
	if (recursionLimit == 0) {
		fprintf(stderr, "Recursion limit reached\n");
		return Q;
	}
	*/
	if (error < tolerance) {
		return Q;
	} else {
		return recursiveAdaptiveIntegrate(f, a, (a + b)/3, delta/sqrt(3), epsilon, func1, recursionLimit - 1)
			+ recursiveAdaptiveIntegrate(f , (a + b)/3, 2*(a + b)/3, delta/sqrt(3), epsilon, func2, recursionLimit - 1)
			+ recursiveAdaptiveIntegrate(f, 2*(a + b)/3, b, delta/sqrt(3), epsilon, func3, recursionLimit - 1);
	}
}

double integrate(double f(double), double a, double b, double delta, double epsilon){
	//
	double func2 = f(a + 3*(b - a)/6);
	//
	int recursionLimit = 99;
	//
	return recursiveAdaptiveIntegrate(f, a, b, delta, epsilon, func2, recursionLimit);
}
