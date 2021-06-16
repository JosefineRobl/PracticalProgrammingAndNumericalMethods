#include <stdio.h>
#include <math.h>

/*
 * Integrates a function with one or more limits possibly being infinity, using three sub-divisions of the integration interval instead of two, and an adaptive and recursive technique.
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
