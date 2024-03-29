#include<stdio.h>
#include<assert.h>
#include<gsl/gsl_vector.h>
#include<math.h>
#include<gsl/gsl_interp.h>
#include<gsl/gsl_spline.h>
#include"binarySearch.h"

typedef struct{gsl_vector* x;
               gsl_vector* y;
               gsl_vector* b;
               gsl_vector* c;
               gsl_vector* d;
} cubicSpline;

/*
 * Allocates memory for the cubic spline.
 */
cubicSpline* cubicSplineAlloc(gsl_vector* x, gsl_vector* y);

/*
 * Evaluating cubic spline at given z.
 */
double cubicSplineEvaluate(cubicSpline* s, double z);

/*
 * Evaluate the derivative of the cubic spline at given z.
 */
double cubicSplineEvaluateDerivative(cubicSpline* s, double z);

/*
 * Evaluates the integral from x[0] to z.
  */
double cubicSplineIntegrate(cubicSpline* s, double z);

//Function for freeing allocated memory
void cubicSplineFree(cubicSpline* s);
