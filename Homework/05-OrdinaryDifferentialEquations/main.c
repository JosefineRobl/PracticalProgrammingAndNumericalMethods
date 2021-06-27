#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rungeKutta12.h"

/*
 * Print exercise title/description.
 */
void printExercise(char* exercise, FILE* out){
	fprintf(out, "=============== %s ===============\n", exercise);
}

/*
 * Harmonic oscillator.
 */
void harmonicOscillator(double x, double* y, double* dydx){
	dydx[0] = y[1];
	dydx[1] = -y[0];
}

/*
 * The SIR model of epidemic development, disregarding natural birth and death.
 */
void SIR(double time_between_contacts, double t, double y[], double dydt[]){
	double population_size = 6*1e6;
	double recovery_time = 14;

	dydt[0] = -y[0]*y[1]/population_size/time_between_contacts;
	dydt[1] = y[0]*y[1]/population_size/time_between_contacts - y[1]/recovery_time;
	dydt[2] = y[1]/recovery_time;
}

/*
 * File-scope variable for the function SIRDependingOnTimeBetweenContacts.
 */
static double timeBetweenContactsFileVar;

/*
 * The SIR model depending on different values of the time between contacts.
 */
void SIRDependingOnTimeBetweenContacts(double t, double y[], double dydt[]){
	SIR(timeBetweenContactsFileVar, t, y, dydt);
}

/*
 * File-scope variables for the function threeBodyProblem.
 */
static double G, M1, M2, M3;

void threeBodyProblem(double x, double* y, double* dydx){
	// Vector y has 12 entries: position (x,y) and velocity (vx,vy) for the three different masses
	double x1 = y[0], y1 = y[1],
		x2 = y[2], y2 = y[3],
		x3 = y[4], y3 = y[5],
		vx1 = y[6], vy1 = y[7],
		vx2 = y[8], vy2 = y[9],
		vx3 = y[10], vy3 = y[11];
	
	double dx1dt = vx1, dy1dt = vy1,
	       dx2dt = vx2, dy2dt = vy2,
	       dx3dt = vx3, dy3dt = vy3;
	
	double x12 = x2 - x1, y12 = y2 - y1,
	       x13 = x3 - x1, y13 = y3 - y1,
	       x23 = x3 - x2, y23 = y3 - y2;
	
	double r12 = sqrt(pow(x12, 2) + pow(y12, 2)),
	       r13 = sqrt(pow(x13, 2) + pow(y13, 2)),
	       r23 = sqrt(pow(x23, 2) + pow(y23, 2));
	
	double f12 = G*M1*M2/pow(r12, 2),
	       f13 = G*M1*M3/pow(r13, 2),
	       f23 = G*M2*M3/pow(r23, 2);
	
	// Newton's Second Law
	double dvx1dt = 1./M1 * (f12*x12/r12 + f13*x13/r13),
	       dvy1dt = 1./M1 * (f12*y12/r12 + f13*y13/r13),
	       dvx2dt = 1./M2 * (-f12*x12/r12 + f23*x23/r23),
	       dvy2dt = 1./M2 * (-f12*y12/r12 + f23*y23/r23),
	       dvx3dt = 1./M3 * (-f13*x13/r13 - f23*x23/r23),
	       dvy3dt = 1./M3 * (-f13*y13/r13 - f23*y23/r23);
	
	dydx[0] = dx1dt; dydx[1] = dy1dt;
	dydx[2] = dx2dt; dydx[3] = dy2dt;
	dydx[4] = dx3dt; dydx[5] = dy3dt;
	dydx[6] = dvx1dt; dydx[7] = dvy1dt;
	dydx[8] = dvx2dt; dydx[9] = dvy2dt;
	dydx[10] = dvx3dt; dydx[11] = dvy3dt;
}

/*
 * Main function.
 */
int main(void){
	// Open file for writing
	FILE* out = fopen("out.txt", "w");
	
	// Testing implementation using u''=-u
	printExercise("Testing implementation using u''=-u", out);
	// Initial conditions
	int n = 2, // Dimension of array
	    a = 0, // Starting point for ODE
	    b = 6; // Ending point for ODE
	double yHarmonic[] = {1, 0}, // Value of u(a) and u'(a) respectively
	       h = 0.2,
	       delta = 1e-3,
	       epsilon = 1e-3;
	// Create file for writing data to
	FILE* filePointerHarmonic = fopen("harmonic.txt", "w");
	// Perform ODE
	rungeKuttaDrive12(n, harmonicOscillator, (double) a, (double) b, yHarmonic, h, delta, epsilon, filePointerHarmonic);
	// Close file for writing
	fclose(filePointerHarmonic);
	// Printing to the out.txt file
	fprintf(out, "Solution to the ODE u''(x)=-u(x) can be seen in harmonicPlot.png alongside its derivative. The initial values are u(0) = 1, u'(0) = 0. The data can be found in the file harmonic.txt.\n");
	
	// SIR model of epidemic development
	printExercise("SIR model of epidemic development", out);
	// Different times between contacts
	double time_between_contacts[] = {1, 3, 5, 10}; // days
	n = sizeof(time_between_contacts) / sizeof(time_between_contacts[0]); // Length of the array
	// Reasonable parameters for Denmark
	a = 0; // We start from 0 days
	b = 100; // And goes to 100 days (first lockdown time)
	double y[n];
	h = 0.1;
	delta = 1e-3;
	epsilon = 1e-3;
	// Open files to write
	FILE* filePointerSusceptible = fopen("susceptible.txt", "w");
	FILE* filePointerInfectious = fopen("infectious.txt", "w");
	FILE* filePointerRemoved = fopen("removed.txt", "w");
	FILE* exerciseBSIR = fopen("exerciseBSIR.txt", "w");
	// Calculate ODE during time
	for (int t = 1; t < b; t++) {
		// Set the appropriate initial values for Denmark
		y[0] = 6*1e6; // Population of Denmark
		y[1] = 10; // 10 people are infected from the beginning (whose who came how from a bar in Tyrol)
		y[2] = 0; // 0 people are dead or recovered from COVID-19 at beginning of epidemi
		// Printing the day number in each file
		fprintf(filePointerSusceptible, "%i\t", t);
		fprintf(filePointerInfectious, "%i\t", t);
		fprintf(filePointerRemoved, "%i\t", t);
		// For different values of contact time
		for (int i = 0; i < n; i++) {
			// Set file-scope variable to the i'th T_C value
			timeBetweenContactsFileVar = time_between_contacts[i];
			// Perform ODE
			rungeKuttaDrive12(n, SIRDependingOnTimeBetweenContacts, (double) a, (double) t, y, h, delta, epsilon, exerciseBSIR);
			// Print the found values from the ODE
			fprintf(filePointerSusceptible, "%15g\t", y[0]);
			fprintf(filePointerInfectious, "%15g\t", y[1]);
			fprintf(filePointerRemoved, "%15g\t", y[2]);
		}	
		// New line in all files
		fprintf(filePointerSusceptible, "\n");
		fprintf(filePointerInfectious, "\n");
		fprintf(filePointerRemoved, "\n");
	}
	// Close files
	fclose(filePointerSusceptible);
	fclose(filePointerInfectious);
	fclose(filePointerRemoved);
	fclose(exerciseBSIR);
	// Print exercise considerations
	fprintf(out, "The graph of the SIR model for the COVID-19 epidemi in Denmark can be seen on sirPlot.png, and the data in the susceptible.txt, infected.txt, and removed.txt files.\n");
	fprintf(out, "The parameters used is:\n \t Population of Denmark: %g\n \t Number of infected from beginning: %d\n \t Number of dead and recoverede from beginning: %d\n", 6*1e6, 10, 0);
	fprintf(out, "From the plot it can be seen, that an decreasing contact time (T_C) give rise to faster and steaper peaks for the infected and number of dead/recovered graphs. This makes perfect sense due to the fact, that there for shorter times between contacts will be more contact and thus a larger virus spreading.\n");
	
	// Newtonian gravitational three-body problem
	printExercise("Newtonian gravitational three-body problem", out);
	// Define variables
	n = 12;
	a = 0;
	b = 6;
	double yThreeBodyProblem[n];
	h = 0.1;
	delta = 1e-3;
	epsilon = 1e-3;
	// Set file-scope variables
	G = 1; M1 = 1; M2 = 1; M3 = 1;
	// Set initial conditions
	yThreeBodyProblem[0] = -0.97000436; // x1
	yThreeBodyProblem[1] = 0.24308753; // y1
	yThreeBodyProblem[2] = 0; // x2
	yThreeBodyProblem[3] = 0; // y2
	yThreeBodyProblem[4] = 0.9700436; // x3
	yThreeBodyProblem[5] = -0.24308753; // y3
	yThreeBodyProblem[6] = 0.4662036850; // vx1
	yThreeBodyProblem[7] = 0.4323657300; // vy1
	yThreeBodyProblem[8] = -0.93240737; // vx2
	yThreeBodyProblem[9] = -0.86473146; // vy2
	yThreeBodyProblem[10] = 0.4662036850; // vx3
	yThreeBodyProblem[11] = 0.4323657300; // vy3
	// Open file for writing
	FILE* filePointerThreeBodyProblem = fopen("threeBodyProblem.txt", "w");
	// Perform ODE solving
	rungeKuttaDrive12(n, threeBodyProblem, (double) a, (double) b, yThreeBodyProblem, h, delta, epsilon, filePointerThreeBodyProblem);
	// Close file for writing
	fclose(filePointerThreeBodyProblem);
	// Printing to the out.txt file
	fprintf(out, "The graph can be found in the file threeBodyProblemPlot.png, while the data can be found in threeBodyProblem.txt.\n");
	
	// Closing file for writing
	fclose(out);

	return 0; 
}
