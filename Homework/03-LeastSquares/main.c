// General
#include <stdio.h>
#include <math.h>
#include <assert.h>
// GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
// Own .h files
#include "leastSquareFunctions.h"

int main(){
	FILE* outFile = fopen("out.txt", "w");
	// Print info
	fprintf(outFile, "GENERAL INFO\n\nFor the answer to exercise A, the fit and the half-life can be found later in this file, alongside the reference value half-life from Wikipedia, and the plot of this is seen in leastSquarePlot.png. Also later in this file, the covariance matrix can be found (exercise B), and on the png-file, the fit is plottet both with and without its uncertainties (exercise C).\n\n\n");
	
	int m = 2;
	double t[] = {1, 2, 3, 4, 6, 9, 10, 13, 15};
	double y[] = {117,100,88,72,53,29.5,25.2,15.2,11.1};
	int n = sizeof(y)/sizeof(y[0]);
	double dy[n];
	double ylog[n];
	double dlogy[n];
	for(int i=0;i<n;i++){
		dy[i] = 0.05 * y[i];
		dlogy[i] = dy[i]/y[i];
		ylog[i] = log(y[i]);
	}
	double func(int i, double x){
		switch(i){
			case 0: return 1; break;
			case 1: return x; break;
			case 2: return x*x; break;
			case 3: return pow(x,3); break;
			case 4: return pow(x,4); break;
			case 5: return pow(x,5); break;
			default: return NAN;
		}
	}
	gsl_vector* lt = gsl_vector_alloc(n);
	gsl_vector* ly = gsl_vector_alloc(n);
	gsl_vector* ldy = gsl_vector_alloc(n);
	for(int i=0;i<n;i++){
		gsl_vector_set(lt,i,t[i]);
		gsl_vector_set(ly,i,ylog[i]);
		gsl_vector_set(ldy,i,dlogy[i]);
	}
	// Allocates GSL matrix and vector for covariance matrix and fit coefficient
	gsl_matrix* cov = gsl_matrix_alloc(m,m);
	gsl_vector* c = gsl_vector_alloc(m);
	// Making the least square fit
	leastSquareFit(m, func, lt, ly, ldy, c, cov);
	gsl_vector* dfit  = gsl_vector_alloc(m);
	// Set the vectors
	for(int k = 0; k < m; k++){
		gsl_vector_set(dfit, k, sqrt(gsl_matrix_get(cov, k, k)));
	}
	// Perform the last part of fit
	double fit(double x){
		double sum =0;
		for(int i=0;i<m;i++){
			sum += gsl_vector_get(c, i)*func(i, x);
		}
		return sum;
	}
	// Writing to the out file about the fit
	fprintf(outFile, "RESULTS\n\nThe polynomial for the fit is given by:\n\t");
	for(int i = 0; i < m; i++) {
		fprintf(outFile, "%g*x^%i)", gsl_vector_get(c, i), i);
	}
	fprintf(outFile, "and thus the value of c1 is %g.\n",gsl_vector_get(c, 1));
	double t0 = -log(2)/gsl_vector_get(c, 1);
	// Writing to the out file about the half-life and the covariance matrix
	fprintf(outFile, "The estimation of the value for the half-life of the THX, which is based upon a least square fit, is %g days.\n Resolving to Wikipedia for the correct value, one gets %g days. This is slightly in agreement with one another.\n\n", t0, 3.6319);
	fprintf(outFile, "Furthermore one could note that the covariance matrix for the calculation is\n");
	fprintf(outFile, "[");
	for (int i = 0; i < m ; i++){
		for(int j = 0; j < m; j++){
		fprintf(outFile, "  %3.g", gsl_matrix_get(cov, i, j));
		}
		fprintf(outFile, "\n");
	}
	fprintf(outFile, " ]");
	// Closing file
	fclose(outFile);
	
	// Write to result file
	FILE* resultFile = fopen("results.txt", "w");
	// Data points
	fprintf(resultFile,"# data points\n");
	for(int i = 0; i < n; i++) {
		fprintf(resultFile,"%g %g %g\n", t[i], ylog[i], dlogy[i]);
	}
	// The fit itself
	fprintf(resultFile,"\n\n# fit\n");
	for(int i = 0; i < n; i++) {
		fprintf(resultFile,"%g %g \n", t[i], fit(t[i]));
	}
	// The fit with the above-fit uncertainties
	fprintf(resultFile,"\n\n# fit with uncertainties= +delta(c0)+delta(c1)\n");
	for(int i = 0; i < n; i++) {
		fprintf(resultFile,"%g %g \n", t[i], fit(t[i]) + gsl_vector_get(dfit, 0)*func(0, t[i]) + gsl_vector_get(dfit, 1)*func(1, t[i]));
	}
	// The fit with the below-fit uncertainties
	fprintf(resultFile,"\n\n# fit with uncertainties = -delta(c0) - delta(c1)\n");
	for(int i = 0; i < n; i++) {
		fprintf(resultFile,"%g %g \n", t[i], fit(t[i]) - gsl_vector_get(dfit, 0)*func(0, t[i]) - gsl_vector_get(dfit, 1)*func(1, t[i]));
	}
	// Close result file
	fclose(resultFile);
	// Free matrices and vectors
	gsl_matrix_free(cov);
	gsl_vector_free(c);
	gsl_vector_free(dfit);
	
	return 0;
}
