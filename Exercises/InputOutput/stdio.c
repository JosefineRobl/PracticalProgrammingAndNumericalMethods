#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "printExercise.h"
#include "formatPrint.h"

/*
 * Reads inputs from the standard input and prints the them alongside their sines and cosines in a table. If no inputs are provided an error message is produces.
 */
void sinCosFromStandardInput(void){
	// Initialize variables
	double x;
	int input = scanf("%lg", &x);
	// Considering if there is or isn't inputs
	if (input != EOF) {
		// If there is inputs, print the headline and while there is still inputs left print the inputs alongside their sines and cosines and update the input variable
		formatHeaderXSinCos();
		//printf("\t x \t sin(x) \t cos(x)\n");
		while (input != EOF) {
			formatTableXSinCos(x, sin(x), cos(x));
			//printf("\t %g \t %g \t %g\n", x, sin(x), cos(x));
			input = scanf("%lg", &x);
		}
	} else {
		// If no inputs an error message is produces
		fprintf(stderr, "There were no inputs given the function, thus it cannot print the inputs alongside their sines and cosines.\n");
	}
}

/*
 * Main fucntion. Calls the sinCosFromStandardInput function.
 */
int main(void){
	printExercise("B");
	sinCosFromStandardInput();
	return 0;
}
