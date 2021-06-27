#include <stdio.h>
#include <math.h>
#include <assert.h>

double integrate2(double f(double), double a, double b, double delta, double epsilon, int variableTransformationFormula);

double ClenshawCurtisIntegrate(double f(double), double a, double b, double delta, double epsilon);

double generalisedIntegrator2(double f(double), double a, double b, double delta, double epsilon);
