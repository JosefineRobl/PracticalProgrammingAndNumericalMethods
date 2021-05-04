#include <stdio.h>
#include <math.h>
#include <tgmath.h>
#include <complex.h>
#include <string.h>

// Definitions were already defined, thus now out commented
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
	printf("%s %s value: %g\n", funName, valType, val);
}

/*
 * Printing the complex numbers as shall be done in Exercise 1.
 *
 * funName: String containing the function name.
 * valType: String containing the type of value, i.e. "found" (computed) or "real" (found on WolframAlpha).
 * val: Complex double containing th function to evaluate or the correct value found online.
 */
void printComplex(char funName[], char valType[], double complex val){
	printf("%s %s value: %g%+gi\n", funName, valType, crealf(val), cimagf(val));
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
int precisionNumber(const char val){
	// Length of the initial segment containing "0.1", i.e. until it does not have 0.1111111... anymore.
	return strspn(val, "0.1");
}

/*
 * Converts a float, a double and a long double to an array of strings the three numbers' digits after the decimal point, which of is taken to be to the 25th decimal.
 *
 * f: Float containing the number as a float.
 * d: Double containing the number as a double.
 * ld: Long double containing the number as a long double.
 *
 * Returns an array of the three strings formed from the decimals of the three given numbers.
 */
const char toString(float f, double d, long double ld){
	// Initialize list of strings and strings
	const char strings[3];
	char floatString[25];
	char doubleString[25];
	char longDoubleString[25];
	// Translate the numbers after the decimal poit
        snprintf(floatString, 25, "%.25g", f*pow(10, 25));
	snprintf(doubleString, 25, "%.25lg", d*pow(10, 25));
	snprintf(longDoubleString, 25, "%.25Lg", ld*pow(10, 25));
	// return the array of strings
	return strings;
}

/*
 * Calculates the value of 1./9 as float, double and long double and prints the result along with the number of digits after decimal separator, which are still precise.
 */
void calculatePrecision(void){
	// Calculating the values of 1./9 as float, double and long double
	float x_float = 1./9;
	double x_double = 1./9;
	long double x_long_double = 1./9;
	// Converts the numbers to strings
	const char strings[3];
	strcpy(strings, toString(x_float, x_double, x_long_double));
	// Index of precision in the computation (index of first number !=1 and -2 due to "0." is of length 2)
	int index_float = precisionNumber(strings[0]);
	int index_double = precisionNumber(strings[1]);
	int index_long_double = precisionNumber(strings[2]);
	// Calls the print function for the precisions
	printPrecision("float", strings[0], index_float);
	printPrecision("double", strings[1], index_double);
	printPrecision("long double", strings[2], index_long_double);
}

unsigned long hash(char *str){
	unsigned long hash = 5381;
	int c;
	while ((c = *str++))
		hash = ((hash << 5) + hash) + c;
	return hash;
}

/*
 * Prints the result of a calculation of 1./9 as well as the precision of the calculation.
 *
 * dataType: A string containing the data type in question, i.e. "float", "double" and "long double".
 * formatSpecifier: A string containing the format specifier for the given dataType, i.e. "g", "lg" and "Lg" for float, double and long double respectively..
 * result: String containing the result of the calculation of 1./9.
 * precisionIndex: Integer of the number of digits after the decimal seperator still precise.
 */
void printPrecision(char dataType[], char result[], int precisionIndex){
	char formatSpecifier[];
	switch ((int)hash(dataType)){
		case (int)hash("float"):
			formatSpecifier = "%.25g";
			float num = strtof(result);
			break;
		case (int)hash("double"):
			formatSpecifier = "%.25lg";
			double num = strtod(result);
			break;
		case (int)hash("long double"):
			formatSpecifier = "%.25Lg";
			long double num = strtold(result);
	}
	printDash();
	printf("1./9 for %s: 1./9 = %s", dataType, formatSpecifier, num*pow(10, -25));
	printf("For %s the precision is down to digit %d after decimal separator", dataType, precisionIndex);
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
	printExercise(2);
	calculatePrecision();
	return 0;
}
