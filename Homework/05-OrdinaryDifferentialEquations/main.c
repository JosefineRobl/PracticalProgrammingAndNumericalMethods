#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "rungeKutta12.h"

void printExercise(char* exercise){
	printf("=============== %s ===============\n", exercise);
}

/*
 * Harmonic oscillator.
 */
void f(int n, double x, double* y, double* dydx){
	dydx[0] = y[1];
	dydx[1] = -y[0];
}

/*
 * The SIR model of epidemic development, disregarding natural birth and death.
 */
void SIR(double timeBetweenContacts, double t, double y[], double dydt[]) {
	double populationSize = 6*1e6; // Number of people in Denmark
	double recoveryTime = 14; // Number of days for recovery from COVID-19
	// Susceptible (those who can be infected)
	dydt[0] = -y[0]*y[1] / (populationSize * timeBetweenContacts);
	// Infected
	dydt[1] = y[0]*y[1] / (populationSize * timeBetweenContacts) - y[1] / recoveryTime;
	// Removed (either recovered with aquired immunity, or dead)
	dydt[2] = y[1] / recoveryTime;
}

static double timeBetweenContactsFileVar;

/*
 * The SIR model depending on different values of the time between contacts.
 */
void SIRDependingOnTimeBetweenContacts(double t, double y[], double dydt[]){
	SIR(timeBetweenContactsFileVar, t, y, dydt);
}

/*
 * Main function.
 */
int main(void){
	// Testing implementation using u''=-u
	printExercise("Testing implementation using u''=-u");
	// ...
	
	// SIR model of epidemic development
	printExercise("SIR model of epidemic development");
	// Different times between contacts
	double timeBetweenContacts[] = {1, 3, 5, 10}; // days
	int n = sizeof(timeBetweenContacts)/sizeof(timeBetweenContacts[0]); // Length of the array
	// Reasonable parameters for Denmark
	int a = 0, // We start from 0 days
	    b = 100; // And goes to 100 days (first lockdown)
	double y[n],
	       h = 0.1,
	       delta = 1e-2,
	       epsilon = 1e-2;
	// Open files to write
	FILE* filePointerSusceptible = fopen("susceptible.txt", "w");
	FILE* filePointerInfectious = fopen("infectious.txt", "w");
	FILE* filePointerRemoved = fopen("removed.txt", "w");
	FILE* exerciseBSIR = fopen("exerciseBSIR.txt", "w");
	for (int t = 1; t < b; t++) {
		// Set the appropriate values for Denmark
		y[0] = 6*1e6; // Population of Denmark
		y[1] = 10; // 10 people are infected from the beginning (those who came from bar in Tyrol)
		y[2] = 0; // 0 people are dead or recovered from COVID-19 at beginning of epidemi
		// For different values of contact times
		for (int i = 0; i < n; i++) {
			timeBetweenContactsFileVar = timeBetweenContacts[i]; // Setting the file-scope variable
			rungeKuttaDrive12(n, SIRDependingOnTimeBetweenContacts, (double) a, (double) t, y, h, delta, epsilon);
			fprintf(filePointerSusceptible, "%g\t", y[0]);
			fprintf(filePointerInfectious, "%g\t", y[1]);
			fprintf(filePointerRemoved, "%g\t", y[2]);
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
	// Print out information to general txt-file
	printf("Hej");
	
	
	// Newtonian gravitational three-body problem
	printExercise("Newtonian gravitational three-body problem");
	// ...
	
	return 0;
}
