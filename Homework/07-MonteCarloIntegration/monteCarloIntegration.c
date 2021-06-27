// Inclusions
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <assert.h>
// Definitions
#define RANDOM (double) rand()/RAND_MAX

/*
 * Prints the exercise title surrounded by equal signs (===).
 *
 * exerciseTitle: String which shall be printed as the exercise title.
 */
void printExerciseTitle(char* exerciseTitle){
	printf("========================= %s =========================\n", exerciseTitle);
}

/*
 * Prints the exercise subtitle surrounded by dashes (---).
 *
 * exerciseSubtitle: String which shall be printed as the exercise subtitle.
 */
void printExerciseSubtitle(char* exerciseSubtitle){
	printf("---------- %s ----------\n", exerciseSubtitle);
}

/*
 * ...
 *
 * dim:
 * f:
 * x:
 * a:
 * b:
 * N:
 *
 * returns ...
 */
complex plainMonteCarloIntegration(int dim, double f(int dim, double* x), double* a, double* b, int N){
	// Initialization of ...
	double V = 1;
	double sum = 0;
	double sumSquared = 0;
	double x[dim];
	//
	for (int i = 0; i < dim; i++) {
		V *= b[i] - a[i];
	}
	//
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < dim; j++) {
			x[j] = a[j] + RANDOM*(b[j] - a[j]);
		}
		// Initializing fx by the function value
		double fx = f(dim, x);
		// Updating the sum by adding the function value
		sum += fx;
		// Updating the squared sum by adding the square of the function value
		sumSquared += fx*fx;
	}
	// Initializing the mean and the variance
	double mean = sum/N;
	double sigma = sqrt(sumSquared/N - mean*mean);
	// Returns the value of the integral
	return mean*V + I*sigma*V/sqrt(N);
}

/*
 * Corput sequence
 */
double corput(int n, int base){
	double q = 0;
	double bk = 1./base;
	while (n > 0) {
		q += (n % base)*bk; // % is the modulo operator
		n /= base;
		bk /= base;
	}
	return q;
}

/*
 * The two Halton sequences with different bases. We use the difference between these two as error estimation.
 */
void halton(int n, int dim, double* x, double* a, double* b, int haltonVersion){
	// The base numbers must be coprimes, thus we choose prime numbers
	int base1[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 47, 53, 59, 61, 67, 71, 73, 79},
	    base2[] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 47, 53, 59, 61, 67, 71, 73, 79, 83},
	    base[21];
	assert(haltonVersion == 1 || haltonVersion == 2);
	switch (haltonVersion) {
		case 1:
			for (int i = 0; i < sizeof(base1)/sizeof(int); i++) {
				base[i] = base1[i];
			}
			break;
		case 2:
			for (int i = 0; i < sizeof(base2)/sizeof(int); i++) {
				base[i] = base2[i];
			}
	}
	int maxDim = sizeof(base)/sizeof(int);
	// The dimension shall be at max the maximal dimension
	assert(dim <= maxDim);
	for (int i = 0; i < dim; i++) {
		x[i] = a[i] + corput(n + 1, base[i]) *(b[i] - a[i]);
	}
}

/*
 * Quasi-random Monte-Carlo integration.
 */
complex quasiMonteCarloIntegration(int dim, double f(int dim, double* x), double* a, double* b, int N){
	// Initialize list with dimension dim
	double x[dim], x2[dim];
	// Initializations
	double V = 1;
	double sum1 = 0, sum2 = 0; // Sums for halton1 and halton2 function respectively
	// Setting V
	for (int i = 0; i < dim; i++) {
		V *= b[i] - a[i];
	}
	// Setting sums
	for (int i = 0; i < N/2; i++) {
		halton(i, dim, x, a, b, 1);
		halton(i, dim, x2, a, b, 2);
		sum1 += f(dim, x);
		sum2 += f(dim, x2);
	}
	double integral = V*(sum1 + sum2)/N;
	double error = V*fabs(sum1 - sum2)/N;

	return integral + I*error;
}


/*
 * Prints the integral alongside the calculated and the theoretical solution.
 *
 * integral: String containing the integral to be solved.
 * calculatedSolution: Complex double containint the theoretical result. Can also take a regular double value instead.
 * assumedSolution: Complex double containing the theoretical result. Can also take a regular double value instead.
 */
void printSolutionToIntegral(char* integral, complex calculatedSolution, complex assumedSolution){
	// Prints the integral to be solved
	printf("The integral to solve is: %s.\n", integral);
	// Initializes the text for the print statements
	char* textCalculatedSolution = "The calculated solution to the integral is";
	char* textAssumedSolution = "The theoretical solution to the integral is";
	// Prints the calculated solution to the integral. Prints without imaginary part if this is zero.
	if (cimagf(calculatedSolution) == 0) {
		printf("%s: %g.\n", textCalculatedSolution, creal(calculatedSolution));
	} else {
		printf("%s: %g%+gi.\n", textCalculatedSolution, crealf(calculatedSolution), cimagf(calculatedSolution));
	}
	// Prints the theoretical solution to the integral. Prints without imaginary part if this is zero.
	if (cimagf(assumedSolution) == 0){
		printf("%s: %g.\n", textAssumedSolution, creal(assumedSolution));
	} else {
		printf("%s: %g%+gi.\n", textAssumedSolution, crealf(assumedSolution), cimagf(assumedSolution));
	}
}

/*
 * Calculating some interesting integrals with the Monte-Carlo routine.
 */
void testIntegralsForPartAB(int monteCarloMethod){
	assert(monteCarloMethod == 1 || monteCarloMethod == 2);
	
	// Initializes two strings with the integral and its solution respectively
	char* integralOverSqrtX = "∫0π  dx/√(x)";
	double integralOverSqrtXApproxSolution = 2*sqrt(M_PI);
	// Initializes the function for integration, the integration limits, the number of points and the dimensions
	double functionOverSqrtX(int d, double* t) {
		return 1./sqrt(t[0]);
	}
	double a[] = {0};
	double b[] = {M_PI};
	int N = 1e7;
	int dim = 1;
	// Prints the integral as the exercise subtitle
	printExerciseSubtitle(integralOverSqrtX);
	// Calculating the integral
	complex integralOverSqrtXCalculatedSolution;
	switch (monteCarloMethod) {
		case 1:
			integralOverSqrtXCalculatedSolution = plainMonteCarloIntegration(dim, functionOverSqrtX, a, b, N);
			break;
		case 2:
			integralOverSqrtXCalculatedSolution = quasiMonteCarloIntegration(dim, functionOverSqrtX, a, b, N);
	}
	// Prints the integral, the calculated solution to the integral and the theoretical solution to the integral
	printSolutionToIntegral(integralOverSqrtX, integralOverSqrtXCalculatedSolution, integralOverSqrtXApproxSolution);

	// Initializes two strings with the integral and its solution respectively
	char* integralComplexExponential = "∫0π exp(ix) dx";
	complex integralComplexExponentialApproxSolution = 2*I;
	// Initializes the function for integration, the integration limits, the number of points and the dimensions
	double functionComplexExponential(int d, double* t) {
		return cexp(I*t[0]);
	}
	double aComplexExponential[] = {0};
	double bComplexExponential[] = {M_PI};
	N = 1e7;
	dim = 1;
	// Prints the integral as the exercise subtitle
	printf("\n");
	printExerciseSubtitle(integralComplexExponential);
	// Calculating the integral
	complex integralComplexExponentialCalculatedSolution;
	switch (monteCarloMethod) {
		case 1:
			integralComplexExponentialCalculatedSolution = plainMonteCarloIntegration(dim, functionComplexExponential, aComplexExponential, bComplexExponential, N);
			break;
		case 2:
			integralComplexExponentialCalculatedSolution = quasiMonteCarloIntegration(dim, functionComplexExponential, aComplexExponential, bComplexExponential, N);
	}
	// Prints the integral, the calculated solution to the integral and the theoretical solution to the integral
	printSolutionToIntegral(integralComplexExponential, integralComplexExponentialCalculatedSolution, integralComplexExponentialApproxSolution);
}

void testPartAB(int monteCarloMethod){
	assert(monteCarloMethod == 1 || monteCarloMethod == 2);

	// Calculating som interesting integrals with the Monte-Carlo routine
	testIntegralsForPartAB(monteCarloMethod);

	printf("\n");

	// Initializes two strings with the singular integral and its solution respectively
	char* singularIntegral = "∫0π dx/π ∫0π dy/π ∫0π dz/π 1/[1 - cos(x)cos(y)cos(z)]";
	double singularIntegralApproxSolution = pow(tgamma(1./4), 4)/(4*pow(M_PI, 3));
	// Initializes the function for integration, the integration limits, the number of points and the dimensions
	double functionSingularIntegral(int d, double* t) {
		return 1./((1 - cos(t[0])*cos(t[1])*cos(t[2]))*pow(M_PI, d));
	}
	double a[] = {0, 0, 0};
	double b[] = {M_PI, M_PI, M_PI};
	int N = 1e7;
	int dim = 3;
	// Prints the singular integral as the exercise subtitle
	printExerciseSubtitle(singularIntegral);
	// Calculating the integral
	complex singularIntegralCalculatedSolution;
	switch (monteCarloMethod) {
		case 1:
			singularIntegralCalculatedSolution = plainMonteCarloIntegration(dim, functionSingularIntegral, a, b, N);
			break;
		case 2:
			singularIntegralCalculatedSolution = quasiMonteCarloIntegration(dim, functionSingularIntegral, a, b, N);
	}
	// Prints the integral, the calculated solution to the integral and the theoretical solution to the integral
	printSolutionToIntegral(singularIntegral, singularIntegralCalculatedSolution, singularIntegralApproxSolution);
}

double funForErrorComparison(int dim, double* r){
	assert(dim == 2);
	double x = r[0];
	double y = r[1];
	return pow(cos(x), 2)*pow(sin(y), 2);
}

/*
 * Testing error scaling.
 */
void errorComparisonAB(void){	
	FILE* error = fopen("error.txt", "w");
	int dim = 2,
	    N = 1e6;
	double a[dim],
	       b[dim];
	for (int i = 0; i < dim; i++) {
		a[i] = 0;
		b[i] = 2*M_PI;
	}
	for(int n = 1000; n < 100000; n += 1000){
		complex resultPlainMC = plainMonteCarloIntegration(dim, funForErrorComparison, a, b, N);
		complex resultQuasiMC = quasiMonteCarloIntegration(dim, funForErrorComparison, a, b, N);
		fprintf(error, "%20d %20g %20g\n", n, cimag(resultPlainMC), cimag(resultQuasiMC));
	}
	printf("The scaling of the errors for the plain monte-carlo vs the quasi-random monte-carlo as a function of N are shown in error.png\n");
	fclose(error);
}

int main(int argc, char** argv){
	// If argument passed update error-file. If argument not passed create out.txt as normal
	if (argc > 1) {
		errorComparisonAB();
	} else {
		printExerciseTitle("Exercise A");
		printf("Integrals calculated using the plain Monte-Carlo integrator.\n");
		testPartAB(1);
		
		printf("\n");
		printExerciseTitle("Exercise B");
		printf("The same integrals as those from exercise A is tested using the quasi-random Monte-Carlo integrator.\n");
		testPartAB(2);
		printExerciseSubtitle("Error comparison");
		printf("Comparison of the scaling of the error for the quasi-random Monte-Carlo integrator (exercise B) with the pseudo-random Monte-Carlo integrator (exercise A) can be seen at error.png (with the data in error.txt). The error scaling is done with the integral of cos^2(x)*sin^2(y) from 0 to 2*PI for both limits.\n");
		
		printf("\n");
		printExerciseTitle("Exercise C");
	}
	return 0;
}
