#include <stdio.h>
#include <math.h>

#include "formatPrint.h"

/*
 * Prints the formattet header for the table.
 */
void formatHeaderXSinCos(void){
	printf("\t x \t sin(x) \t cos(x)\n");
	printf("      ------------------------------------\n");
}

/*
 * Prints the formattet table entries.
 *
 * x: Double containing the x value.
 * sin: Double containing the sin(x).
 * cos: Double containing the cos(x).
 */
void formatTableXSinCos(double x, double sin, double cos){
	if ((sin == fabs(1)) | (sin == fabs(0))) {
		printf("\t %g \t %g \t\t %g\n", x, sin, cos);
	} else {
		printf("\t %g \t %g \t %g\n", x, sin, cos);
	}
}
