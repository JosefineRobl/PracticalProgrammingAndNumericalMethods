#include <assert.h>
#include <gsl/gsl_vector.h>
#include "binarySearch.h"

typedef struct{
	gsl_vector* x;
	gsl_vector* y;
	gsl_vector* b;
	gsl_vector* c;
} quadraticSpline;

quadraticSpline* quadraticSplineAlloc(gsl_vector* x, gsl_vector* y);

/*
 * Evaulates interpolation in point z (see table in notes).
 */
double quadraticSplineEval(quadraticSpline* s, double z);

/*
 * Function that evaluates the derivative in point z.
 */
double quadraticSplineEvalDerivative(quadraticSpline* s, double z);

/*
 * Function that evaluates integral from x[0] to z.
 */
double quadraticSplineIntegrate(quadraticSpline* s, double z);

//Function for freeing allocated memory
void quadraticSplineFree(quadraticSpline* s);
