#include <stdio.h>
#include <math.h>
#include <assert.h>

double ClenshawCurtisIntegrate(double f(double), double a, double b, double delta, double epsilon);

double generalisedIntegrator(double f(double), double a, double b, double delta, double epsilon);
