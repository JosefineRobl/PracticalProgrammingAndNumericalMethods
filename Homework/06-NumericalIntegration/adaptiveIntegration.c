#include <stdio.h>
#include <math.h>
#include <assert.h>

double recursiveIntegrate(double f(double), double a, double b, double delta, double epsilon, double func2, double func3, int recursionLimit){
	// Initialization of functions surrounding func2 and func3 using 'Numerical Integration' PDF eq. 51 and 48
	double func1 = f(a + (b - a)/6);
	double func4 = f(a + 5*(b - a)/6);
	// Initialization of higher and lower order rules (Q and q respectively)
	double Q = (2*func1 + func2 + func3 + 2*func4)*(b - a)/6;
	double q = (func1 + func2 + func3 + func4)*(b - a)/4;
	// Initialization of the error and the tolerance
	double error = fabs(Q - q);
	double tolerance = delta + epsilon*fabs(Q);
	if (recursionLimit == 0) {
		fprintf(stderr, "Function 'recursiveIntegrate' have reached the recursion limit.\n");
		return Q;
	}
	if (error < tolerance) {
		return Q;
	} else {
		return recursiveIntegrate(f, a, (a + b)/2, delta/sqrt(2), epsilon, func1, func2, recursionLimit - 1) + recursiveIntegrate(f, (a + b)/2, b, delta/sqrt(2), epsilon, func3, func4, recursionLimit - 1);
	}
}

double integrate(double f(double), double a, double b, double delta, double epsilon){
	// Initialized func2 and func3 from eq. 51 and 48 in the 'Numerical Integration' PDF
	double func2 = f(a + 2*(b - a)/6);
	double func3 = f(a + 4*(b - a)/6);
	// Initialize the limit of recursions to 99 times - the 100th time should resolve in an error
	int recursionLimit = 99;
	// Begin recursion
	return recursiveIntegrate(f, a, b, delta, epsilon, func2, func3, recursionLimit);
}

/*
void testIntegrate(double delta, double epsilon){
	int calls = 0;
	double f(double x){
		calls++;
		return sqrt(x);
	}
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
	double h(double x){
		calls++;
		return 4*sqrt(1 - pow(x, 2));
	}
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
	// assert(sqrtXIntegral == 2./3);
	// assert(4sqrt1MinusXSquaredIntegral == M_PI);
	printf("Test of function 'integrate' was successful.\n");
}
*/

/*
double generalisedIntegrator(double f(double), double a, double b, double delta, double epsilon){
	double g(double);
	// Check if a and/or b is INFINITY (isinf returns 0 for not infinity and nonzero for argument being infinity)
	if (isinf(a) != 0) {
		if (isinf(b) != 0) {
			// If both a and b are INFINITY: Using converted integral from table in notes (eq. 58)
			g(double t){return ( f((1 - t)/t) + f(- (1 - t)/t ))/pow(t, 2);}
			integrate(g, 0, 1, delta, epsilon);
		} else {
			// If only a is INFINITY: Using converted integral from table in notes (eq. 62)
			g(double t){return ( f(b - (1 - t)/t) )/pow(t, 2);}
			integrate(g, 0, 1, delta, epsilon);
		}
	} else if (isinf(b) != 0) {
		// If only b is INFINITY: Using converted integral from table in notes (eq. 60)
		g(double t){return ( f(a + (1 - t)/t) )/pow(t, 2);}
		integrate(g, 0, 1, delta, epsilon);
	} else {
		// If neither a nor b are INFINITY: Regular integrate function from question A
		integrate(f, a, b, delta, epsilon);
	}

}
*/

int calls;
double f(double x){
	calls++;
	return sqrt(x);
}

double h(double x){
	calls++;
	return 4*sqrt(1 - pow(x, 2));
}

int main(){
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
	//
	//testIntegrate(0.001, 0.001);
	return 0;
}
