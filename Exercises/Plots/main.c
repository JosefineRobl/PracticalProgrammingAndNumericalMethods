#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include <assert.h>

#include "erfApproximation.h"
#include "gammaApproximation.h"

void printGamma(double x){
	printf("%10g\t%10g\t%10g\t%10g\n", x, tgamma(x), gsl_sf_gamma(x), gammaStirlingApproximation(x));
}

int main(int argc, char** argv){
	// The main function shall take 1 argument being 1 for error function and 2 for gamma function
	assert(argc == 2);
	long int arg = strtol(argv[1], NULL, 10); // Converts the string number from the argument to an integer (long)
	if (arg == 1) {
		for (double x = 0; x <= 3; x += 1./8) {
			printf("%10g\t%10g\t%10g\t%10g\n", x, erf(x), gsl_sf_erf(x), erfApproximation(x));
		}
	} else if (arg == 2) {
		double x;
		for (x = -4.9; x <= -4.1; x += 1./8) {
			printGamma(x);
		}
		for (x = -3.9; x <= -3.1; x += 1./8) {
			printGamma(x);
		}
		for (x = -2.9; x <= -2.1; x += 1./8) {
			printGamma(x);
		}
		for (x = -1.9; x <= -1.1; x += 1./8) {
			printGamma(x);
		}
		for (x = -0.82; x <= -0.18; x += 1./8) {
			printGamma(x);
		}
		for (x = 0.2; x <= 4.1; x += 1./8) {
			printGamma(x);
		}
	}
	return 0;
}
