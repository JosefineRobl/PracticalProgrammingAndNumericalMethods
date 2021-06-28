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
	// ....
	}
}

