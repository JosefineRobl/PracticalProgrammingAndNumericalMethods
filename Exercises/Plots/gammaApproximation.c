#include <stdio.h>
#include <math.h>

double gammaStirlingApproximation(double x){
	// single precision gamma function (Gergo Nemes, from Wikipedia)
	if (x < 0) {
		return M_PI/sin(M_PI*x)/gammaStirlingApproximation(1 - x);
	} else if (x < 9) {
		return gammaStirlingApproximation(x + 1)/x;
	}
	double lnGammaStirlingApproximation = x*log(x + 1/(12*x - 1/x/10)) - x + log(2*M_PI/x)/2;
	return exp(lnGammaStirlingApproximation);
}
