#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "adaptiveIntegration2.h"

static double A, B; // the left and right integration limit

/*
static double gOnlyLowerLimitInf(double f(double), double t){
	return f( B - (1 - t)/t ) / pow(t, 2);
}
*/

static double gOnlyUpperLimitInf(double f(double),double t){// variable transformation formula
	return f( A + (1 - t)/t ) / pow(t, 2);
}

static double gOnlyLowerLimitInf(double f(double), double t){
	return -gOnlyUpperLimitInf(f, t);
}

static double gBothLimitsInf(double f(double), double t){
	return ( f((1-t)/t) + f(-(1-t)/t) ) / pow(t, 2);
}

double wrap24( double f(double),double a, double b,
	double delta, double epsilon, double f2, double f3, int nrec, int variableTransformationFormula){//wrapper calls F(f,t)
	assert(nrec > 0);
	double xValLowerLimit = a + (b - a)/6,
	       xValUpperLimit = a + 5*(b - a)/6;
	double f1, f4;
	switch (variableTransformationFormula) {
		case 1:
			f1 = gOnlyLowerLimitInf(f, xValLowerLimit);
			f4 = gOnlyLowerLimitInf(f, xValUpperLimit);
			break;
		case 2:
			f1 = gOnlyUpperLimitInf(f, xValLowerLimit);
			f4 = gOnlyUpperLimitInf(f, xValUpperLimit);
			break;
		case 3:
			f1 = gBothLimitsInf(f, xValLowerLimit);
			f4 = gBothLimitsInf(f, xValUpperLimit);
			break;
		default:
			f1 = f(xValLowerLimit);
			f4 = f(xValUpperLimit);
	}
	double Q = (2*f1 + f2 + f3 + 2*f4) / 6 * (b - a);
	double q = (f1 + f2 + f3 + f4) / 4 * (b - a);
	double tolerance = delta + epsilon*fabs(Q);
	double error = fabs(Q - q);
	if(error < tolerance) return Q;
	else {
		double Q1 = wrap24(f, a, (a+b)/2, delta/sqrt(2), epsilon, f1, f2, nrec - 1, variableTransformationFormula);
		double Q2 = wrap24(f, (a+b)/2, b, delta/sqrt(2), epsilon, f3, f4, nrec - 1, variableTransformationFormula);
		return Q1 + Q2;
	}
}

double integrate2(double f(double),double a, double b, double delta,double epsilon, int variableTransformationFormula){
	//A=a; B=b; // save left limit in file-scope variable
	double xValLowerLimit = a + 2*(b - a)/6,
	       xValUpperLimit = a + 4*(b - a)/6;
	double f2,f3;
	switch (variableTransformationFormula) {
		case 1:
			A = -b;
			f2 = gOnlyLowerLimitInf(f, xValLowerLimit);
			f3 = gOnlyLowerLimitInf(f, xValUpperLimit);
			break;
		case 2:
			A = a;
			f2 = gOnlyUpperLimitInf(f, xValLowerLimit);
			f3 = gOnlyUpperLimitInf(f, xValUpperLimit);
			break;
		case 3:
			A = a; B = b;
			f2 = gBothLimitsInf(f, xValLowerLimit);
			f3 = gBothLimitsInf(f, xValUpperLimit);
			break;
		default:
			f2 = f(xValLowerLimit);
			f3 = f(xValUpperLimit);
	}
	int nrec = 99;
	return wrap24(f, a, b, 2*delta, 2*epsilon, f2, f3, nrec, variableTransformationFormula);
}

double generalisedIntegrator2(double f(double), double a, double b, double delta, double epsilon){
	int variableTransformationFormula;
	if (isinf(a) != 0) {
		if (isinf(b) != 0) {
			variableTransformationFormula = 3;
			a = 0; b = 1;
		} else {
			variableTransformationFormula = 1;
			a = 0; b = 1;
		}
	} else if (isinf(b) != 0) {
		variableTransformationFormula = 2;
		a = 0; b = 1;
	} else {
		variableTransformationFormula = 0;
	}
	return integrate2(f, a, b, delta, epsilon, variableTransformationFormula);
}
