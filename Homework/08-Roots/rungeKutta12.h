#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*
 * Stepper.
 * 
 * f: The function from dy/dt = f(t, y).
 * n: Integer containing the size of the vector y.
 * x: The current value of the variable.
 * yx: The current value y(t) of the sought function.
 * h: Double containing the step to be taken.
 * yh: The y(t+h)-value.
 * dy: The error esimate.
 */
void rungeKuttaStep23(void f(int n, double x, double* y, double* dydx), int n, double x, double* yx, double h, double* yh, double* dy);

/*
 * Driver.
 * 
 */
void rungeKuttaDrive23(void f(int n, double x, double* y, double* dydx), int n, double a, double b, double* ya, double* yb, double h, double delta, double epsilon, char* path)
