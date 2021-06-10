#include <stdio.h>
#include <math.h>
#include "exponential.h"

int main(void){
	FILE *filePointer = fopen("data.txt", "w");

	double min = -5.01; 
	double max = 5.01; 
	double epsilon = 1./8; 
	
	for (double x = min; x <= max; x += epsilon) {
		fprintf(filePointer, "%10g %10g %10g\n", x, exp(x), exponential(x));
	}

	fclose(filePointer);
	
	return 0; 
}
