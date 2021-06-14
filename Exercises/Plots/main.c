#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include <assert.h>

#include "erfApproximation.h"
#include "gammaApproximation.h"

int main(int argc, char** argv){
	// The main function shall take 1 argument being 1 for error function and 2 for gamma function
	assert(argc == 2);
	long int arg = strtol(argv[1], NULL, 10); // Converts the string number from the argument to an integer (long)
	if (arg == 1) {
		for (double x = 0; x <= 3; x += 1./8) {
			printf("%10g\t%10g\t%10g\t%10g\n", x, erf(x), gsl_sf_erf(x), erfApproximation(x));
		}
	} else if (arg == 2) {
		void printGamma(double x){
			printf("%10g\t%10g\t%10g\t%10g\n", x, tgamma(x), gsl_sf_gamma(x), gammaStirlingApproximation(x));
		}
		double x;
		for (x = -4.998; x <= -4.000001; x += 1./8) {
			printGamma(x);
		}
		for (x = -3.992; x <= -3.001; x += 1./8) {
			printGamma(x);
		}
		for (x = -2.97; x <= -2.04; x += 1./8) {
			printGamma(x);
		}
		for (x = -1.91; x <= -1.15; x += 1./8) {
			printGamma(x);
		}
		for (x = -0.81; x <= -0.14; x += 1./8) {
			printGamma(x);
		}
		for (x = 0.15; x <= 4; x += 1./8) {
			printGamma(x);
		}
	}
	return 0;
}
