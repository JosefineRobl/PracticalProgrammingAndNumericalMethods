#include <stdio.h>
#include <math.h>
#include "komplex.h"

/*
 * Generates a new komplex number x+yI.
 * 
 * x: Double containing the real part of the komplex number.
 * y: Double containing the imaginary part of the komplex number.
 */
komplex komplexNew(double x, double y){
	/*
	// My implementation such that it actually initializes a komplex number.
	struct komplex z;
	z.re = x;
	z.im = y;
	return z;
	//return (komplex){.re = x, .im = y};
	*/
	// Dmitri's implementation, which won't let me use it for initialization. Here one should stille write the below line for initialization in the main script, so why have this function?
	komplex z = {x, y};
	return z;
}

/*
 * Sets a komplex number to x+yI.
 * 
 * z: Pointer to a komplex number.
 * x: Double containing the new value for the real part of the komplex number.
 * y: Double containing the new value for the imaginary part of the komplex number.
 */
void komplexSet(komplex* z, double x, double y){
	(*z).re = x;
	(*z).im = y;
}

/*
 * Prints a komplex number as "<s> (a.re,a.im)".
 * 
 * s: String containing what shall be printet in front of the complex number.
 * a: Komplex number.
 */
void komplexPrint(char* s, komplex a){
	printf("%s (%g,%g).\n", s, a.re, a.im);
}

/*
 * Adds two komplex numbers together.
 * 
 * a: Komplex number to be added together with another komplex number.
 * b: Komplex number to be added to another komplex number.
 */
komplex komplexAdd(komplex a, komplex b){
	return komplexNew(a.re + b.re, a.im + b.im);
}

/*
 * Subtracts two komplex numbers from each other.
 * 
 * a: Komplex number to have another komplex number subtract from it.
 * b: Komplex number to be subtracted from another komplex number.
 */
komplex komplexSubtract(komplex a, komplex b){
	komplexAdd(a, -b);
}

/*
 * Takes the norm of a komplex number.
 * 
 * z: Komplex number to be taken the norm of.
 */
double komplexNorm(komplex z){
	return sqrt(pow(z.re, 2) + pow(z.im, 2));
}

/*
 * Conjugates a komplex number.
 *
 * z: Komplex number to be conjugated.
 */
komplex komplexConjugate(komplex z){
	return komplexNew(z.re, -z.im);
}

/*
 * Multiplies two komplex numbers together.
 * 
 * a: Komplex number to be multiplied together with another complex number.
 * b: Komplex number to be multiplied to another complex number.
 */
komplex komplexMultiply(komplex a, komplex b){
	return komplexNew(a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re);
}

/*
 * Divides a two komplex numbers.
 * 
 * a: Komplex number to be divided by another komplex number.
 * b: Komplex number to divide another komplex number with.
 */
komplex komplexDivide(komplex a, komplex b){
	double normSquareOfDenominator = komplexNorm(b);
	komplex newNumerator = komplexMultiply(a, komplexConjugate(b));
	return komplexNew(newNumerator.re / normSquareOfDenominator, newNumerator.im / normSquareOfDenominator);
}
