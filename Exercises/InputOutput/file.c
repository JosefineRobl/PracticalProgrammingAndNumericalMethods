#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "printExercise.h"
#include "formatPrint.h"

/*
 * Reads inputs from an input file and then prints the inputs alongside their sines and cosines in a table. If no inputs are present in the file, the function prints out a message saying, that the file is empty.
 *
 * argv: List where send item (argv[1]) contains the path and name of the input file.
 */
void sinCosFromInputFile(char** argv){
	// Loads the file
	FILE *filePointer = fopen(argv[1], "r");
	// Initializes variables
	double x;
	int input = fscanf(filePointer, "%lg", &x);
	// Checks if the file is empty
	if (input != EOF) {
		// If the file is not empty print the header and print the inputs alongside their sines and cosines, updating the input after printing
		formatHeaderXSinCos();
		//printf("\t x \t sin(x) \t cos(x)\n");
		while (input != EOF) {
			formatTableXSinCos(x, sin(x), cos(x));
			//printf("\t %g \t %g \t %g\n", x, sin(x), cos(x));
			input = fscanf(filePointer, "%lg", &x);
		}
	} else {
		// If file is empty, print message.
		printf("Input file is empty, thus a table cannot be created.\n");
	}
	fclose(filePointer);
}

/*
 * Main function taking one argument. Passes its arguments to the sinCosFromInputFile function, which it calls.
 *
 * argc: Integer containing the number of arguments.
 */
int main(int argc, char** argv){
	printExercise("C");
	if (argc < 2) {
		fprintf(stderr, "Too few arguments. 0 arguments have been specified, thus no file name is specified, and the reader cannot read from an unspecified file. Plase specify a single filename as an argument for the function.\n");
	} else if (argc > 2) {
		fprintf(stderr, "Too many arguments. %d arguments have been specified, thus there is more arguments than only the filename (1 argument) and the function cannot be sure to choose the correct file. Please specify only 1 arguemnt for the function, that being the file name.\n", argc-1);
	} else {
		sinCosFromInputFile(argv);
	}
	return 0;
}
