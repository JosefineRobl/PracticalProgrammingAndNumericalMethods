#include <stdio.h>
#include <math.h>

/*
 * Stepper.
 * 
 * n: Integer containing the size of the vector y.
 * f: The function from dy/dt = f(t, y).
 * t: The current value of the variable.
 * y: The current value y(t) of the sought function.
 * h: Double containing the step to be taken.
 * step: The y(t+h)-value.
 * dy: The error esimate.
 */
void rungeKuttaStep12(int n, void f(double t, double y[], double dydt[]), double t, double y[], double h, double step[], double dy[]);

/*
 * Driver.
 * 
 * n: Integer containing the size of the vector y.
 * f: The function from dy/dt = f(t, y).
 * a: The start of y.
 * b: The end of y.
 * y: The current value y(t) of the sought function.
 * h: Double containing the step to be taken.
 * delta: Double containing the absolute accuracy goal.
 * epsilon: Double containing the relative accuracy goal.
 * out: Pointer to FILE object to store the points calculated along the way.
 */
void rungeKuttaDrive12(int n, void f(double t, double y[], double dydt[]), double a, double b, double y[], double h, double delta, double epsilon);
