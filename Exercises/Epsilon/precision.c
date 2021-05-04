/*
#ifndef "precision.h" // Once-Only Header (ensures that file is not read multiple times resulting in erros). See https://www.tutorialspoint.com/cprogramming/c_header_files.htm
#define "precision.h" // Once-Only Header
#endif // Once-Only Header
*/

#include <stdio.h>
#include <math.h>
#include "precision.h"

/*
 * Function dettermining if numbers 'a' and 'b' are equal with absolute or relative precision.
 *
 * a: Double
 * b: Double
 * tau: Double
 * epsilon: Double
 */
int equal(double a, double b, double tau, double epsilon){
	return ((fabs(a - b) < tau) || ((fabs(a - b) / ( fabs(a) + fabs(b))) < (epsilon / 2))) ? 1 : 0;
}

/*
 * Prints the cases for the tests along with the result.
 *
 * a: Double
 * b: Double
 * tau: Double
 * epsilon: Double
 * testText: String containing explanatory text about the predicted result of the test case.
 */
void printEqualTest(double a, double b, double tau, double epsilon, char testText[]){
	printf("%s: a = %g, b = %g, tau = %g, epsilon = %g, equal(%g, %g, %g, %g) = %i.\n", testText, a, b, tau, epsilon, a, b, tau, epsilon, equal(a, b, tau, epsilon));
}

/*
 * Tests cases for the fucntion 'equal'.
 */
void testEqual(void){
	printf("We test equal(a, b, tau, epsilon):\n");
	printEqualTest(2, 1, 1, 1, "Only relative precision should be true"); // Only relative precision (epsilon).
	printEqualTest(2, 1, 2, 1, "Both absolute and relative precision should be true"); // Both absolute (tau) and relative precision (epsilon).
	printEqualTest(2, 1, 2, 1./4, "Only absolute precision should be true"); // Only absolute precision (tau).
	printEqualTest(2, 1, 1./2, 1./4, "Neither precision should be true, i.e. (both should be false)"); // Both should fail.
}
