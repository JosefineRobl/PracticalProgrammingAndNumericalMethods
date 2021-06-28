#include <stdio.h>
#include <math.h>

// Definition
#define FORALL(i) for (int i = 0; i < n; i++)
#define PRINT(t, y) fprintf(out, "%15g\t", t); FORALL(i) fprintf(out, "%15g\t", y[i]); fprintf (out, "\n")

/*
 * Stepper.
 * 
 * n: Integer containing the size of the vector y.
 * f: The function from dy/dt = f(t, y).
 * t: The current value of the variable.
 * y: The current value y(t) of the sought function.
 * h: Double containing the step to be taken.
 * step: The y(t+h)-value.
 * dy: The error esimate.
 */
void rungeKuttaStep12(int n, void f (double t, double y[], double dydt[]), double t, double y[], double h, double step[], double dy[]){
	double k0[n], k1[n], y1[n];
	// Evaluate the function to determine k0
	f(t, y, k0);
	FORALL(i) y1[i] = y[i] + (0.5*h)*k0[i];
	// Evaluate the function to determine k1
	f(t + 0.5*h, y1, k1);
	FORALL(i) {
		step[i] = y[i] + h*k1[i];
		dy[i] = (k1[i] - k0[i])*h;
	}
}

/*
 * Driver.
 * 
 * n: Integer containing the size of the vector y.
 * f: The function from dy/dt = f(t, y).
 * a: The start of y.
 * b: The end of y.
 * y: The current value y(t) of the sought function.
 * h: Double containing the step to be taken.
 * delta: Double containing the absolute accuracy goal.
 * epsilon: Double containing the relative accuracy goal.
 * out: Pointer to FILE object to store the points calculated along the way.
 */
void rungeKuttaDrive12(int n, void f (double t, double y[], double dydt[]), double a, double b, double y[], double h, double delta, double epsilon, FILE* out){
	// Initialize the time t as the starting point
	double t = a;
	// Print values of t and y before run (exercise part B)
	PRINT(t, y);
	// Calculate as long as we have not reached the last point yet
	while (t < b) {
		if (t + h > b) {
			h = b - t;
		}
		// Initialize vectors and summation variables
		double steps[n], dy[n],
		       sumForNorm = 0, sumForError = 0;
		// Perform step with stepper from above
	      	rungeKuttaStep12 (n, f, t, y, h, steps, dy);
		// Calculate the norm of y, the error and the rolerance
		FORALL(i) {
			sumForNorm += pow(y[i], 2);
			sumForError += pow(dy[i], 2);
		}
		double error = sqrt(sumForError); 
		double tolerance = (delta + epsilon*sqrt(sumForNorm)) * sqrt(h / (b - a));
		// Compare the error with the tolerance, and just the error by itself
		if (error < tolerance) {
			t = t + h;
			// Update the y-list
			FORALL(i) y[i] = steps[i];
			// Print values of t and y after run (exercise part B)
			PRINT(t, y);
		}
		if (error > 0) {
			h *= 0.95*pow(tolerance/error, 0.25); 
		} else 	{
			h *= 2;
		}
	}
}
