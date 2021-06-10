#include <stdio.h>
#include <math.h>

#include "formatPrint.h"

/*
 * Prints the formatted header for the table.
 */
void formatHeaderXSinCos(void){
	printf("\t x \t sin(x) \t cos(x)\n");
	printf("      ------------------------------------\n");
}

/*
 * Prints the formatted header for the table to a specific file.
 *
 * fileStream: Pointer to the FILE object, that identifies the stream to which the printing should occur.
 */
void formatHeaderXSinCosFprint(FILE *fileStream){
	fprintf(fileStream, "\t x \t sin(x) \t cos(x)\n");
	fprintf(fileStream, "      ------------------------------------\n");
}

/*
 * Prints the formatted table entries.
 *
 * x: Double containing the x value.
 * sin: Double containing the sin(x) value.
 * cos: Double containing the cos(x) value.
 */
void formatTableXSinCos(double x, double sin, double cos){
	if ((sin == fabs(1)) | (sin == fabs(0))) {
		printf("\t %g \t %g \t\t %g\n", x, sin, cos);
	} else {
		printf("\t %g \t %g \t %g\n", x, sin, cos);
	}
}

/*
 * Prints the formatted table entries to a specific file.
 * 
 * x: Double containing the x value.
 * sin: Double containing the sin(x) value.
 * cos: Double containing the cos(x) value.
 * fileStream: Pointer to the FILE object, that identifies the stream to which the printing should occur.
 */
void formatTableXSinCosFprint(double x, double sin, double cos, FILE *fileStream){
	if ((sin == fabs(1)) | (sin == fabs(0))) {
		fprintf(fileStream, "\t %g \t %g \t\t %g\n", x, sin, cos);
	} else {
		fprintf(fileStream, "\t %g \t %g \t %g\n", x, sin, cos);
	}
}
