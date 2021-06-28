#include <stdio.h>
#include <math.h>
#include <assert.h>
// GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
// Own files
#include "linearSpline.h"
#include "cubicSpline.h"
#include "quadraticSpline.h"
#include "binarySearch.h"

/*
 * EXERCISE A - Linear spline.
 */
void exerciseA(void){
	// First the points to interpolate
	int N = 9; //number of tabulated points
	gsl_vector* x = gsl_vector_alloc(N); 
	gsl_vector* y = gsl_vector_alloc(N);
	
	// Open files for the x and y points
	FILE* x_file = fopen("xPoints.txt","r");
	FILE* y_file = fopen("yPoints.txt","r");
	gsl_vector_fscanf(x_file, x);
	gsl_vector_fscanf(y_file, y);
	FILE* xy_file = fopen("xyPoints.txt","w");

	// For GSL comparison we must use ordinary arrays:
	double xa[x -> size];
	double ya[y -> size];
	for(int i=0; i< x -> size; i++){
		 xa[i] = gsl_vector_get(x,i);
		 ya[i] = gsl_vector_get(y,i);
	}
	for(int i=0; i<=N-1; i++){
		fprintf(xy_file, "%10g %10g\n", gsl_vector_get(x,i), gsl_vector_get(y,i));
	}
	
	// For GSL functions
	gsl_interp* linear = gsl_interp_alloc(gsl_interp_linear, N);
	gsl_interp_init(linear, xa, ya, N);
	// Open file to write result
	FILE* linterp_file = fopen("exerciseA.txt", "w");
	// Then the interpolation function
	int z=0;
	double fine = 0.1;
	while(z*fine <= gsl_vector_get(x, N-1)){
		double interp_l_gsl = gsl_interp_eval(linear, xa, ya, z*fine, NULL); // GSL lin. interp.
		double interp_integ_gsl = gsl_interp_eval_integ(linear, xa, ya, xa[0], z*fine, NULL); //GSL integral
		fprintf(linterp_file, "%10g %10g %10g %10g %10g\n",z*fine, linearInterpolation(x, y, z*fine), linearInterpolationIntegration(x, y, z*fine), interp_l_gsl, interp_integ_gsl);
		z++;

	}
	
	// Free GSL vectors and interpreters
	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_interp_free(linear);
	
	// Close files
	fclose(x_file);
	fclose(y_file);
	fclose(xy_file);
	fclose(linterp_file);
}

/*
 * EXERCISE B - Cubic spline.
 */
void exerciseB(void){
	int N = 9;
	gsl_vector* x = gsl_vector_alloc(N);
	gsl_vector* y = gsl_vector_alloc(N);
	FILE* x_file = fopen("xPoints.txt", "r");
	FILE* y_file = fopen("yPoints.txt", "r");
	gsl_vector_fscanf(x_file,x);
	gsl_vector_fscanf(y_file,y);
	
	// For GSL comparison we must have ordinary arrays:
	double xa[x->size];
	double ya[y->size];
	for(int i = 0; i < x->size; i++){
		xa[i] = gsl_vector_get(x, i);
		ya[i] = gsl_vector_get(y, i);
	}
	cubicSpline* s = cubicSplineAlloc(x,y);
	
	// For GSL functions
	gsl_interp* cspline = gsl_interp_alloc(gsl_interp_cspline,N);
	gsl_interp_init(cspline, xa, ya, N);
	// File to write exeercise B
	FILE* cinterp_file = fopen("exerciseB.txt", "w");
	// Comparison
	int z=0;
	double fine = 0.1;
	while(z*fine <= gsl_vector_get(x,N-1)){
		double interp_eval_gsl = gsl_interp_eval(cspline, xa, ya, z*fine, NULL);
		double interp_der_gsl = gsl_interp_eval_deriv(cspline, xa, ya, z*fine, NULL);
		double interp_integ_gsl = gsl_interp_eval_integ(cspline, xa, ya,xa[0], z*fine, NULL);
		fprintf(cinterp_file, "%10g %10g %10g %10g %10g %10g %10g\n", z*fine, cubicSplineEvaluate(s, z*fine), cubicSplineEvaluateDerivative(s, z*fine), cubicSplineIntegrate(s, z*fine), interp_eval_gsl, interp_der_gsl, interp_integ_gsl);
		z++;
	}

	// Closing files
	fclose(x_file);
	fclose(y_file);
	fclose(cinterp_file);
	// Freeing GSL vectors
	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_interp_free(cspline);
	// Freeing memory of cubic spline
	cubicSplineFree(s);
}

/*
 * EXERCISE C - Quadratic spline.
 */
void exerciseC(void){
	int N = 9;
	// Make GSL vectors
	gsl_vector* x = gsl_vector_alloc(N);
	gsl_vector* y = gsl_vector_alloc(N);
	// Open files for writing
	FILE* x_file = fopen("xPoints.txt","r");
	FILE* y_file = fopen("yPoints.txt", "r");
	gsl_vector_fscanf(x_file,x);
	gsl_vector_fscanf(y_file,y);
	
	// Create quadratic spline
	quadraticSpline* s = quadraticSplineAlloc(x,y);
	// Open file to write to
	FILE* qinterp_file = fopen("exerciseC.txt", "w");
	// GSL
	int z = 0;
	double fine = 0.1;
	while(z* fine <= gsl_vector_get(x,N-1)){
		fprintf(qinterp_file, "%10g %10g %10g %10g\n", z*fine, quadraticSplineEval(s, z*fine), quadraticSplineEvalDerivative(s, z*fine), quadraticSplineIntegrate(s,z*fine));
		z++;
	}
	
	// Closing files
	fclose(x_file);
	fclose(y_file);
	fclose(qinterp_file);
	// Freeing memory
	quadraticSplineFree(s);
	gsl_vector_free(x);
	gsl_vector_free(y);
}

int main(void){
	exerciseA();
	exerciseB();
	exerciseC();
	return 0;
}
