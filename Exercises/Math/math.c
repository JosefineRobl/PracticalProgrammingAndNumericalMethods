#include <stdio.h>
#include <math.h>
#include <tgmath.h>
#include <complex.h>
#include <string.h>

// Definitions were already defined, thus now out-commented
// #define M_E 2.7182818284590452354 // e
// #define M_PI 3.1415926535897932384 // pi

/*
 * Printing dashes for separation.
 */
void printDash(void){
	printf("----------\n");
}

/* 
 * Printing the exercise number with dash separation.
 */
void printExercise(int n){
	printf("---------- Exercise %d ", n);
}

/*
 * Printing the real numbers as shall be done in Exercise 1.
 *
 * funName: String containing the function name
 * valType: String containing the type of value, i.e. "found" or (computed) or "real" (found on WolframAlpha)
 * val: Double containing the function to evaluate or the correct value found online.
 */
void printReal(char funName[], char valType[], double val){
	printf("%s %s value: %g.\n", funName, valType, val);
}

/*
 * Printing the complex numbers as shall be done in Exercise 1.
 *
 * funName: String containing the function name.
 * valType: String containing the type of value, i.e. "found" (computed) or "real" (found on WolframAlpha).
 * val: Complex double containing th function to evaluate or the correct value found online.
 */
void printComplex(char funName[], char valType[], double complex val){
	printf("%s %s value: %g%+gi.\n", funName, valType, crealf(val), cimagf(val));
}

/*
 * Printing the numbers as shall be done in Exercise 1. Chooses which of the functions printReal and printComplex to use.
 * 
 * funName: String containing the function name.
 * fun: Complex double containing the real or complex function to be computed.
 * realVal: Complex double containing the correct value of the computation of the function fun found online.
 */
void print(char funName[], double complex fun, double complex realVal){
	printDash(); // Printing dashed above all different functions for separations
	// Caculating which print function to use for the function
	if (cimagf(fun) == 0*I)
		printReal(funName, "found", fun);
	else
		printComplex(funName, "found", fun);
	// Calculating which print fucntion to use for the real value of the computation found online
	if (cimagf(realVal) == 0*I)
		printReal(funName, "real", realVal);
	else
		printComplex(funName, "real", realVal);
}


/*
 * Calculated the precision of a computation of 1./9.
 *
 * val: Imutable string containing the value of the calculation of 1./9.
 *
 * Return the integer number of digits after decimal point which are still precise.
 */
int precisionNumber(char* str){
	// Initialization of counter
	int counter = 0;
	// Updating the counter with the number of 1's in the number
	while (str[counter] == '1') counter++;
	// Returns the counter
	return counter;
}

/*
 * Converts a float, a double and a long double to an array of strings containing the three numbers digits after the decimal point, which of is taken to be to the 25th decimal.
 *
 * f: Float containing the number as a float.
 * d: Double containing the number as a double.
 * ld: Long double containing the number as a long double.
 *
 * Returns an array of the three strings formed from the decimals of the three given numbers.
 */
/* NO LONGER IN USE!!!
const char toString(float f, double d, long double ld){
	// Initialize list of strings and strings
	char floatString[28];
	char doubleString[28];
	char longDoubleString[28];
	// Translate the numbers to strings
        snprintf(*floatString, 28, "%.25g", f);
	snprintf(*doubleString, 28, "%.25lg", d);
	snprintf(*longDoubleString, 28, "%.25Lg", ld);
	// Getting strings containing only the numbers after the decimal point, which strings will also be "searchable" by str[i]
	char* floatStringAfterDecimal = &floatString[2];
	char* doubleStringAfterDecimal = &doubleString[2];
	char* longDoubleStringAfterDecimal = &longDoubleString[2];
	// return the array of strings
	const char* strings[] = {floatStringAfterDecimal, doubleStringAfterDecimal, longDoubleStringAfterDecimal};
	return strings;
}
*/

/*
 * Hash function.
 * 
 * str: String which shall be hashed.
 * 
 * returns the hash of str.
 */
unsigned long hash(char* str){
	unsigned long hash = 5381;
	int c;
	while ((c = *str++))
		hash = ((hash << 5) + hash) + c; // hash * 33 + c
	return hash;
}

/*
 * Prints the result of a calculation of 1./9 as well as the precision of the calculation.
 *
 * dataType: A string containing the data type in question, i.e. "float", "double" and "long double".
 * result: String containing the result of the calculation of 1./9.
 * precisionIndex: Integer of the number of digits after the decimal seperator still precise.
*/
void printPrecision(char* dataType, char* result, int precisionIndex){
	char* num;
	// String comparison without hash function may also be viable, but to be sure of the equality, the hash function have been used. This is also the remains from the time where the if-statement were switch-cases instead, thus the hash function was needed to distinguish the cases
	if (hash(dataType) == hash("float")) {
		num = "1.f/9";
	} else if (hash(dataType) == hash("double")) {
		num = "1./9";
	} else if (hash(dataType) == hash("long double")) {
		num = "1.L/9";
	}
	printDash();
	printf("1./9 for %s: %s = 0.%s.\n", dataType, num, result);
	printf("For %s the precision is down to digit %d after decimal separator.\n", dataType, precisionIndex);
}

/*
 *  * Calculates the value og 1./9 as float, double and long double, and prints the result along with the number of digits after the decimal separator, which are still precise.
 *   */
void calculatePrecision(void){
	// Calculating the values of 1./9 as float, double and long double
	float xFloat = 1.f/9;
	double xDouble = 1./9;
	long double xLongDouble = 1.L/9;
	// Initialize list of strings
	char floatString[28];
	char doubleString[28];
	char longDoubleString[28];
	// Translate the numbers to strings
	snprintf(floatString, 28, "%.25g", xFloat);
	snprintf(doubleString, 28, "%.25lg", xDouble);
	snprintf(longDoubleString, 28, "%.25Lg", xLongDouble);
	// Getting strings containing only the numbers after the decimal point, which will also be searchable by str[i]
	char* floatStringAfterDecimal = &floatString[2];
	char* doubleStringAfterDecimal = &doubleString[2];
	char* longDoubleStringAfterDecimal = &longDoubleString[2];
	// Index of precision in the computation (index of the first number !=1, thus the actual number of 1's in the calculation)
	int indexFloat = precisionNumber(floatStringAfterDecimal);
	int indexDouble = precisionNumber(doubleStringAfterDecimal);
	int indexLongDouble = precisionNumber(longDoubleStringAfterDecimal);
	// Calls the print function for the precisions
	printPrecision("float", floatStringAfterDecimal, indexFloat);
	printPrecision("double", doubleStringAfterDecimal, indexDouble);
	printPrecision("long double", longDoubleStringAfterDecimal, indexLongDouble);
}

/*
 * Main fucntion containing all the function calls for Exercise 1 and Exercise 2
 *
 * Returns 0 for correct execution
 */
int main(void){
	// Exercise 1
	printExercise(1);
	print("Gamma(5)", expf(gamma(5)), 24);
	print("Bessel j_1(0.5)", j1(0.5), 0.242268);
	print("sqrt(-2)", csqrt(-2), 1.41421356*I);
	print("exp(i*pi)", cexp(I*M_PI), -1);
	print("exp(i)", cexp(I), 0.54030231+0.84147098*I);
	print("i^e", pow(I, M_E), -0.42821977-0.90367462*I);
	print("i^i", pow(I, I), 0.20787958);
	
	// Exercise 2
	printf("\n");
	printExercise(2);
	calculatePrecision();
	
	return 0;
}
