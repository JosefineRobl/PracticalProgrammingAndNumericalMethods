#ifndef HAVE_KOMPLEX_H // This is necessary for multiple includes
#define HAVE_KOMPLEX_H

/*
 * Komplex number.
 */
typedef struct{
	double re; // Real part of komplex number
	double im; // Imaginary part of komplex number
} komplex;

/*
 * Generates a new komplex number x+yI.
 * 
 * x: Double containing the real part of the komplex number.
 * y: Double containing the imaginary part of the komplex number.
 */
void komplexNew(double x, double y);

/*
 * Sets a komplex number to x+yI.
 * 
 * z: Pointer to a komplex number.
 * x: Double containing the new value for the real part of the komplex number.
 * y: Double containing the new value for the imaginary part of the komplex number.
 */
void komplexSet(komplex* z, double x, double y);

/*
 * Prints a komplex number as "<s> (a.re,a.im)".
 * 
 * s: String containing what shall be printet in front of the complex number.
 * a: Komplex number.
 */
void komplexPrint(char* s, komplex a);

/*
 * Adds two komplex numbers together.
 * 
 * a: Komplex number to be added together with another komplex number.
 * b: Komplex number to be added to another komplex number.
 */
komplex komplexAdd(komplex a, komplex b);

/*
 * Subtracts two komplex numbers from each other.
 * 
 * a: Komplex number to have another komplex number subtract from it.
 * b: Komplex number to be subtracted from another komplex number.
 */
komplex komplexSubtract(komplex a, komplex b);

/*
 * Takes the norm of a komplex number.
 * 
 * z: Komplex number to be taken the norm of.
 */
double komplexNorm(komplex z);

/*
 * Conjugates a komplex number.
 *
 * z: Komplex number to be conjugated.
 */
komplex komplexConjugate(komplex z);

/*
 * Multiplies two komplex numbers together.
 * 
 * a: Komplex number to be multiplied together with another complex number.
 * b: Komplex number to be multiplied to another complex number.
 */
komplex komplexMultiply(komplex a, komplex b);

/*
 * Divides a two komplex numbers.
 * 
 * a: Komplex number to be divided by another komplex number.
 * b: Komplex number to divide another komplex number with.
 */
komplex komplexDivide(komplex a, komplex b);

#endif
