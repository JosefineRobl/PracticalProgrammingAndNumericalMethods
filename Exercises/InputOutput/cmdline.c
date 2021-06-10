#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "printExercise.h"
#include "formatPrint.h"

/*
 * Prints a table of the arguments along with their sines and cosines. If no arguments are given, the function resolves in an error message.
 *
 * argc: Integer containing the argument count.
 * argv: List of the arguments given the function.
 */
void sinCosOfArguments(int argc, char** argv){
	// Prints the name of the program
	printf("The name of the program is %s\n", argv[0]);
	// If there is passed no arguments (argc == 1, since function name is always passed), the program writes out an error message
	if (argc < 2) {
		fprintf(stderr, "No arguments were given, thus the program cannot print the sines and cosines.\n");
	} else {
		// There is arguments, thus these are printed alongside their sines and cosines
		formatHeaderXSinCos();
		//printf("\t x \t sin(x) \t cos(x)\n");
		for(int i = 1; i < argc; i++) {
			double x = atof(argv[i]);
			formatTableXSinCos(x, sin(x), cos(x));
			//printf("\t %g \t %g \t %g\n", x, sin(x), cos(x));
		}
	}
}

/*
 * Main function. Calls the sinCosOfArguments function with the arguments given the main function itself.
 *
 * argc: Integer containing the number of arguments.
 * argv: List containing the arguments themselves.
 */
int main(int argc, char** argv) {
	printExercise("A");
	
	printSubtext("Explanation");
	printf("If no arguments are given, the program will write out an error message saying, that it is not possible for the function to calculate the sines and cosines since no arguments are given.\n");
       printf("If arguments are given, these are printed in a table alongside their sines and cosines.\n");
	
	printSubtext("Result of run");
	sinCosOfArguments(argc, argv);

	return 0;
}
