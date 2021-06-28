#include <assert.h>
#include <gsl/gsl_vetor.h>
#include "binarySearch.h"

typedef struct{
	gsl_vector* x;
	gsl_vector* y;
	gsl_vector* b;
	gsl_vector* c;
} quadraticSpline;

quadraticSpline* quadraticSplineAlloc(gsl_vector* x, gsl_vector* y){
	// Initialize struct
	quadraticSpline* spline = (quadraticSpline*)malloc(sizeof(quadraticSpline)); // See theory in the book, table 3
	// Gets the size of the vectors
	int N = x->size;
	// Allocate space for vectors
	spline->x = gsl_vector_alloc(N); // N points
	spline->y = gsl_vector_alloc(N); // N points
	spline->b = gsl_vector_alloc(N - 1); // N points => N-1 intervals with polynomials
	spline->c = gsl_vector_alloc(N - 1); // N points => N-1 intervals with polynomials
	gsl_vector* h = gsl_vector_alloc(N - 1);
	gsl_vector* p = gsl_vector_alloc(N - 1);
	// Filling vectors spline->x and spline->y
	for (int i = 0; i < N; i++) {
		// Getting the i'th element of vectors x and y
		double xi = gsl_vector_get(x, i);
		double yi = gsl_vector_get(y, i);
		// Sets i'th number in vector spline->x to xi and in vector spline->y to yi
		gsl_vector_set(spline->x, i, xi);
		gsl_vector_set(spline->y, i, yi);
	}
	// Filling vectors h and p
	for (int i = 0; i < N - 1; i++) {
		double deltaX = gsl_vector_get(x, i + 1) - gsl_vector_get(x, i);
		double deltaY = gsl_vector_get(y, i + 1) - gsl_vector_get(y, i);
		gsl_vector_set(h, i, deltaX);
		gsl_vector_set(p, i, deltaY/deltaX);
	}
	// Finding c_i by recursion
	// Recursion up
	gsl_vector_set(spline->c, 0, 0); // This can be chosen freely (see notes)
	for(int i = 0; i < N-2; i++){
		double diffp = gsl_vector_get(p, i+1) - gsl_vector_get(p, i);
		double ci = gsl_vector_get(spline->c, i);
		double hi = gsl_vector_get(h, i);
		double hi1 = gsl_vector_get(h, i+1);
		double ci1 = (diffp - ci*hi)/hi1;
		gsl_vector_set(spline->c, i+1, ci1);
	}
	// Recursion down
	gsl_vector_set(spline->c, N-2, gsl_vector_get(s->c, N-2)/2);
	for(int i = N - 3; i >= 0; i--){
		double diffp = gsl_vector_get(p, i+1) - gsl_vector_get(p, i);
		double ci1 = gsl_vector_get(spline->c, i+1);
		double hi1 = gsl_vector_get(h, i);
		double hi = gsl_vector_get(h, i+1);
		double ci = (diffp - ci1*hi1)/hi;
		gsl_vector_set(spline->c, i, ci);
	}
	// Finding b_i
	for(int i = 0; i < N-1; i++){
		double pi = gsl_vector_get(p, i);
		double ci = gsl_vector_get(spline->c, i);
		double hi = gsl_vector_get(h, i);
		double bi = pi - ci*hi;
		gsl_vector_set(spline->b, i, bi);
	}
	return s;
}

/*
 * Evaulates interpolation in point z (see table in notes).
 */
double qubicSplineEval(quadraticSpline* s, double z){
	int i = binsearch(s->x, z);
	double h = z-gsl_vector_get(s->x, i);
	double yi = gsl_vector_get(s->y, i);
	double bi = gsl_vector_get(s->b, i);
	double ci = gsl_vector_get(s->c, i);
	double eval = yi + h*(bi + h*ci);
	return eval;
}

/*
 * Function that evaluates the derivative in point z.
 */
double qinterp_der(quadraticSpline* s, double z){
	int i = binsearch(s->x, z);
	double h = z - gsl_vector_get(s->x, i);
	double bi = gsl_vector_get(s->b, i);
	double ci = gsl_vector_get(s->c, i);
	double derivative = bi + 2*ci*h;
	return derivative;
}

/*
 * Function that evaluates integral from x[0] to z.
 */
double qinterp_integ(quadraticSpline* s, double z){
	int j = binsearch(s->x, z);
	double area = 0;
	// Summing over all intervals except the last one, which is done outside the loop
	for(int i = 0; i < j; i++){
		double h = gsl_vector_get(s ->x, i+1)-gsl_vector_get(s->x, i);
		double yi = gsl_vector_get(s ->y, i);
		double bi = gsl_vector_get(s ->b, i);
		double ci = gsl_vector_get(s ->c, i);
		area += yi*h+1./2*bi*pow(h,2)+1./3*ci*pow(h,3);
	}
	// Summing over the last interval
	double hj = z - gsl_vector_get(s->x, j);
	double yj = gsl_vector_get(s ->y, j);
	double bj = gsl_vector_get(s ->b, j);
	double cj = gsl_vector_get(s ->c, j);
	area += yj*hj + 1./2 * bj*pow(hj, 2) + 1./3 * cj*pow(hj, 3);
	return area;
}

//Function for freeing allocated memory
void qinterp_free(quadraticSpline* s){
	gsl_vector_free(s->x);
	gsl_vector_free(s->y);
	gsl_vector_free(s->b);
	gsl_vector_free(s->c);
	free(s);
}

