#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "printExercise.h"
#include "formatPrint.h"

/*
 * Reads inputs from an input file and then prints the inputs alongside their sines and cosines in a table. If no inputs are present in the file, the function prints out a message saying, that the file is empty.
 *
 * argv: List containing where second item (argv[1]) contains the path and name of the input file and third item (argv[2]) contains the path and name of the output file.
 * outputFilePointer: Pointer to FILE object that identifies with the stream for the output file.
 */
void sinCosFromInputFile(char** argv, FILE *outputFilePointer){
	// Loads the file
	FILE *inputFilePointer = fopen(argv[1], "r");
	// Initializes variables
	double x;
	int input = fscanf(inputFilePointer, "%lg", &x);
	// Checks if the file is empty
	if (input != EOF) {
		// If the file is not empty print the header and print the inputs alongside their sines and cosines, updating the input after printing
		formatHeaderXSinCosFprint(outputFilePointer);
		//printf("\t x \t sin(x) \t cos(x)\n");
		while (input != EOF) {
			formatTableXSinCosFprint(x, sin(x), cos(x), outputFilePointer);
			//printf("\t %g \t %g \t %g\n", x, sin(x), cos(x));
			input = fscanf(inputFilePointer, "%lg", &x);
		}
	} else {
		// If file is empty, print message.
		printf("Input file is empty, thus a table cannot be created.\n");
	}
	fclose(inputFilePointer);
}

/*
 * Main function taking one argument. Passes its arguments to the sinCosFromInputFile function, which it calls.
 *
 * argc: Integer containing the number of arguments.
 */
int main(int argc, char** argv){
	// Loads the output file for writing (and creates it if it doen't already exist)
	FILE *outputFilePointer = fopen(argv[2], "w");

	// Prints the exercise title
	printExerciseFprint("C", outputFilePointer);
	
	// Checks if the correct amount of arguments are present or else writes out an error message.
	if (argc < 2) {
		fprintf(outputFilePointer, "Too few arguments. 0 arguments have been specified, thus no file name is specified, and the reader cannot read from an unspecified file. Plase specify a single filename as an argument for the function.\n");
	} else if (argc > 3) {
		fprintf(outputFilePointer, "Too many arguments. %d arguments have been specified, thus there is more arguments than just the filenames of the input and output files (2 arguments) and the function cannot be sure to choose the correct files. Please specify only 2 arguemnts for the function, that being the filenames of firstly the input file and then of the output file.\n", argc-1);
	} else {
		sinCosFromInputFile(argv, outputFilePointer);
	}
	
	// Closes the file
	fclose(outputFilePointer);

	return 0;
}
