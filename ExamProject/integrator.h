#include <stdio.h>
#include <math.h>
#include <assert.h>

double adaptiveRecursiveIntegrate(double f(double), double a, double b, double delta, double epsilon, double func2, int recursionLimit);

double integrateTridivision(double f(double), double a, double b, double delta, double epsilon);

double generalisedIntegrator(double f(double), double a, double b, double delta, double epsilon);
