#include <stdio.h>
#include "komplex.h"

void printSubtext(char* subtext){
	printf("---------- %s ----------\n", subtext);
}

int main(void){
	// Initialize complex numbers
	komplex z1 = komplexNew(1, 2);
	komplex z2 = komplexNew(4, 6);
	
	// Print the chosen values
	printSubtext("By initialization we should have z1 = 1+2i and z2 = 4+6i:");
	komplexPrint("z1 =", z1);
	komplexPrint("z2 =", z2);
	
	// Set a complex number to a new value
	komplexSet(&z2, 3, 7);
	printSubtext("Setting z2 = 3+7i:");
	komplexPrint("z2 =", z2);

	// TESTING THE ALGEBRAIC KOMPLEX FUNCTIONS
	
	// Norm
	printSubtext("The norm of z1 should be |z1| = sqrt(5) = 2.236067:");
	printf("|z1| = %g\n", komplexNorm(z1));
	
	// Complex conjugate
	printSubtext("The complex conjugate of z2 should be conj(z2) = 3-7i:");
	komplexPrint("conj(z2) =", komplexConjugate(z2));
	
	// Adding
	printSubtext("Adding z1 and z2 should yield z1 + z2 = 4+9i:");
	komplexPrint("z1 + z2 =", komplexAdd(z1, z2));
	
	// Subtracting
	printSubtext("Subtracting z1 from z2 should yield z2 - z1 = 2+5i:");
	komplexPrint("z2 - z1 =", komplexSubtract(z2, z1));
	
	// Multiplication
	printSubtext("Multiplying z1 by z2 should yield z1 * z2 = -11+13i:");
	komplexPrint("z1 * z2 =", komplexMultiply(z1, z2));
	
	// Division
	printSubtext("Dividing z2 by z1 should yield z2 / z1 = (17+i)/5 = 3.4+0.2i:");
	komplexPrint("z2 / z1 =", komplexDivide(z2, z1));
	
	return 0;
}
