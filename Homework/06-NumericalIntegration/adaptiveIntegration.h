#ifndef HAVE_ADAPTIVEINTEGRATION_H
#define HAVE_ADAPTIVEINTEGRATION_H

#include <stdio.h>
#include <math.h>
#include <assert.h>

/*
 * Integrates a function with limits possibly being infinity, using two sub-divisions of the integration interval, and an adaptive and recursive technique.
 * 
 * f: Function to be integrated. Taks a double as input and returns a double (the function value at that point).
 * a: Double containing the lower integration limit.
 * b: Double containing the upper integration limit.
 * delta: Double containing the absolute accuracy goal.
 * epsilon: Double containing the relative accuracy goal.
 * variableTransformationFormula: Integer determining the variable transformation formula to use. 0 if no integration limit is infinity, 1 if lower integration limit is infinity, 2 if upper integration limit is infinity, and 3 if both lower and upper integration limit is infinity.
 * 
 * returns a double containing the recursively and adaptively integrated value of the function in the integration limit.
 */
double integrate(double f(double), double a, double b, double delta, double epsilon, int variableTransformationFormula);

/*
 * Integrates a function with one or more limits possibly being infinity, using two sub-divisions of the integration interval, and an adaptive and recursive technique.
 * 
 * f: Function to be integrated. Takes a double as input and returns a double (the function value at the point).
 * a: Double containing the lower integration limit.
 * b: Double containing the upper integration limit.
 * delta: Double containing the absolute accuracy goal.
 * epsilon: Double containing the relative accuracy goal.
 * 
 * returns a double containing the recursively and adaptively integrated value of the function in the integration limit.
 */
double generalisedIntegrator(double f(double), double a, double b, double delta, double epsilon);

#endif
