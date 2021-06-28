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
} cubicSline;

/*
 * Allocates memory for the cubic spline.
 */
cubicSline* cubicSplineAlloc(gsl_vector* x, gsl_vector* y){
	// Initializing a cubic spline element
	cubicSline* spline = (cinterp*)malloc(sizeof(cinterp));
	int N = x->size;
	spline -> x = gsl_vector_alloc(N);
	spline -> y = gsl_vector_alloc(N);
	spline -> b = gsl_vector_alloc(N);
	spline -> c = gsl_vector_alloc(N-1);
	spline -> d = gsl_vector_alloc(N-1);
	// Initialize more vectors
	gsl_vector *h = gsl_vector_alloc(N-1);
	gsl_vector *p = gsl_vector_alloc(N-1);
	// Set the values of the spline
	for(int i = 0; i < N; i++){
		gsl_vector_set(spline -> x, i, gsl_vector_get(x, i));
		gsl_vector_set(spline -> y, i , gsl_vector_get(y, i));
	}
	//
	for(int i = 0; i < N - 1; i++){
		double diffx = gsl_vector_get(x,i+1) - gsl_vector_get(x,i);
		double diffy = gsl_vector_get(y, i+1) - gsl_vector_get(y,i);
		gsl_vector_set(h, i, diffx);
		gsl_vector_set(p, i, diffy/diffx);
	}
	// Initialize tridiagonal elements as in the notes (the sizes follows by the table)
	gsl_vector *D = gsl_vector_alloc(N);
	gsl_vector *Q = gsl_vector_alloc(N - 1);
	gsl_vector *B = gsl_vector_alloc(N);
	// Updating above initialized elements
	gsl_vector_set(D, 0, 2);
	gsl_vector_set(Q, 0, 1);
	gsl_vector_set(B, 0, 3*gsl_vector_get(p, 0));
	for(int i = 0; i < N - 2; i++){
		double hi = gsl_vector_get(h, i);
		double hi1 = gsl_vector_get(h, i + 1);
		gsl_vector_set(D, i + 1, 2*hi/hi1 + 2);
		gsl_vector_set(Q, i + 1, hi/hi1);
		gsl_vector_set(B, i + 1, 3*(gsl_vector_get(p, i) + gsl_vector_get(p, i + 1)*hi/hi1));
	}
	gsl_vector_set(D, N-1, 2);
	gsl_vector_set(B, N - 1, 3*gsl_vector_get(p, N - 2));

	for(int i = 1; i < N; i++){
		double Di_1 = gsl_vector_get(D, i-1);
		gsl_vector_set(D, i, gsl_vector_get(D, i) - gsl_vector_get(Q, i-1)/Di_1);
		gsl_vector_set(B, i, gsl_vector_get(B,i) - gsl_vector_get(B, i-1)/Di_1);
	}
	// Setting now the spline values with the rest of the vectors
	gsl_vector_set(spline->b, N-1, gsl_vector_get(B, N-1)/gsl_vector_get(D, N-1));
	for(int i = N - 2; i >= 0; i--){
		gsl_vector_set(spline->b, i, (gsl_vector_get(B, i) - gsl_vector_get(Q, i)*gsl_vector_get(spline->b, i+1))/gsl_vector_get(D, i));
	}
	for(int i=0; i<N-1; i++){
		double pi = gsl_vector_get(p,i);
		double hi = gsl_vector_get(h,i);
		double bi = gsl_vector_get(spline->b, i);
		double bi1 = gsl_vector_get(spline->b, i+1);
		gsl_vector_set(spline->c, i, (-2*bi - bi1 + 3*pi)/hi);
		gsl_vector_set(spline->d, i, (bi + bi1 - 2*pi)/pow(hi, 2));
	}
	return s;
}

/*
 * Evaluating cubic spline at given z.
 */
double cubicSplineEvaluate(cubicSline* s, double z){
	int i = binarySearch(s->x, z);
	double h = z-gsl_vector_get(s->x, i);
	double yi = gsl_vector_get(s->y, i);
	double bi = gsl_vector_get(s->b, i);
	double ci = gsl_vector_get(s->c, i);
	double di = gsl_vector_get(s->d, i);
	double eval = yi + h*(bi + h*(ci + h*di));
	return eval;
}

/*
 * Evaluate the derivative of the cubic spline at given z.
 */
double cubicSplineEvaluateDerivative(cubicSline* s, double z){
	int i = binarySearch(s->x, z);
	double h = z-gsl_vector_get(s->x, i);
	double bi = gsl_vector_get(s->b, i);
	double ci = gsl_vector_get(s->c, i);
	double di = gsl_vector_get(s->d, i);
	double derivative = bi + 2*ci*h + 3*di*pow(h, 2);
	return derivative;
}

/*
 * Evaluates the integral from x[0] to z.
  */
double cubicSplineIntegrate(cubicSline* s, double z){
	int j = binarySearch(s->x, z);
	double area = 0;
	for(int i = 0; i < j; i++){
		double h = gsl_vector_get(s->x, i+1) - gsl_vector_get(s->x, i);
		double yi = gsl_vector_get(s->y, i);
		double bi = gsl_vector_get(s->b, i);
		double ci = gsl_vector_get(s->c, i);
		double di = gsl_vector_get(s->d, i);
		area += yi*h+1./2*bi*pow(h,2)+1./3*ci*pow(h,3)+1./4*di*pow(h,4);
	}
	double hj = z - gsl_vector_get(s->x, j);
	double yj = gsl_vector_get(s->y, j);
	double bj = gsl_vector_get(s->b, j);
	double cj = gsl_vector_get(s->c, j);
	double dj = gsl_vector_get(s->d, j);
	area += yj*hj + 1./2 * bj*pow(hj, 2) + 1./3 * cj*pow(hj, 3) + 1./4 * dj*pow(hj,4);
	return area;
}


//Function for freeing allocated memory
void cubicSplineFree(cubicSline* s){
	gsl_vector_free(s->x);
	gsl_vector_free(s->y);
	gsl_vector_free(s->b);
	gsl_vector_free(s->c);
	gsl_vector_free(s->d);
	free(s);
}
