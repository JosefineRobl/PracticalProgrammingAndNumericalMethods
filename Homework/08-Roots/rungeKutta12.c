#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*
 * Stepper.
 * 
 * f: The function from dy/dt = f(t, y).
 * n: Integer containing the size of the vector y.
 * x: The current value of the variable.
 * yx: The current value y(t) of the sought function.
 * h: Double containing the step to be taken.
 * yh: The y(t+h)-value.
 * dy: The error esimate.
 */
void rungeKuttaStep23(void f(int n, double x, double* y, double* dydx), int n, double x, double* yx, double h, double* yh, double* dy) {
	double k1[n], k2[n], k3[n], k4[n];
	double yt[n];
	// Evaluate to find k1
	f(n, x, yx, k1);
	for (int i = 0; i < n; i++) {
		yt[i] = yx[i] + 1./2 * k1[i]*h;
	}
	// Evaluate to find k2
	f(n, x+1./2 * h, yt, k2);
	for (int i = 0; i < n; i++) {
		yt[i] = yx[i] + 3./4 * k2[i]*h;
	}
	// Evaluate to find k3
	f(n, x + 3./4 * h, yt, k3);
	for (int i=0; i < n; i++) {
		yh[i] = yx[i] + (2./9 * k1[i] + 1./3 * k2[i] + 4./9 * k3[i])*h;
	}
	// Evaluate to find k4
	f(n, x + h, yh, k4);
	for (int i=0; i<n; i++){
		yt[i] = yx[i] + (7./24 * k1[i] + 1./4 * k2[i] + 1./3 * k3[i] + 1./8 * k4[i])*h;
		dy[i] = yh[i] - yt[i];
	}
}

/*
 * Driver.
 * 
 */
void rungeKuttaDrive23(void f(int n, double x, double* y, double* dydx), int n, double a, double b, double* ya, double* yb, double h, double delta, double epsilon, char* path){
	// Create file to write to
	FILE* file = fopen(path, "w");
	// Initialize variables (mostly) for the step
	double x;
	double y[n];
	double sumNorm, sumError,
	       error, norm, tolerance;
	double yh[n];
	double dy[n];
	// Initialize the x as the starting point
	x = a;
	// Print values of x and y before run
	fprintf(file, "%20g", x);
	for(int i = 0; i < n; i++){
		y[i] = ya[i];
		fprintf(list, "%20g", y[i]);
	}
	fprintf(file, "\n");
	// Calculate as long as we have not reached the last point yet
	while (x < b){
		if (x + h > b) {
			h = b - x;
		}
		// Perform step with stepper from above
		rkstep23(f, n, x, y, h, yh, dy);
		// Calculate the norm of y, the error and the tolerance
		sumNorm = 0;
		sumError = 0;
		for (int i = 0; i < n; i++){
			sumError += pow(dy[i], 2);
			sumNorm += pow(yh[i], 2);
		}
		error = sqrt(sumError);
		norm = sqrt(sumError);
		tolerance = (norm*epsilon + delta)*sqrt(h/(b - a));
		// Compare the error with the tolerance, and just the error by itself
		if (error < tolerance){
			x += h;
			fprintf(list, "%20g", x);
			for (int i = 0; i < n; i++) {
				// Update the y-list
				y[i] = yh[i];
				// Print values of t and y after run
				fprintf(file, "%20g", y[i]);
			}
			fprintf(file, "\n");
		}
		if (err > 0) {
			h *= 0.95*pow(tol/err, 0.25);
		} else {
			h *= 2;
		}
	}
	for (int i = 0; i < n; i++) {
		yb[i] = yh[i];
	}
	// Close file
	fclose(file);
}
