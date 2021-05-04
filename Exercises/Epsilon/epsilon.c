#include <stdio.h> // Standard In/Out contains printf()
#include <limits.h> // Contains INT_MAX and INT_MIN
#include <float.h> // Contains FLT_EPSILON, DBL_EPSILON and LDBL_EPSILON
#include <string.h> // Contains e.g. strcpy()
#include <stdlib.h> // Contains e.g. exit(), EXIT_FAILURE and EXIT_SUCCESS
#include "precision.h" // File containing Exercise 3

/*
 * Prints a string of dashes. No new line at the end.
 */
void printDashes(void){
	printf("----------");
}

/*
 * Prints the exercise number surrounded by dashes.
 *
 * exercise: String containing the exercise name and number.
 */
void printExercise(char exercise[]){
	printDashes();
	printf(" %s ", exercise);
	printDashes();
	printf("\n");
}

/*
 * Printing "My <max/min> int for the <loop> is <integer value>".
 *
 * maxOrMin: String denoting either "max" or "min" to determine the calculation done.
 * loop: String denoting the kind of loop.
 * integer: The integer found as max/min int. 
 */
void printMaxMinInt(char maxOrMin[], char loop[], int integer){
	printf("My %s int for the %s is %i.\n", maxOrMin, loop, integer);
}

/*
 * Function pointer to an integer operation.
 *
 * x: Integer value to take part in the integer operation.
 * y: Integer value to take part in the integer operation.
 */
int (*intOperation)(int x, int y);

/*
 * Calculates if one integer is less than another.
 *
 * x: Integer value to compare with y.
 * y: Integer value which x is compared to.
 *
 * returns 1 if true and 0 if false
 */
int lessThan(int x, int y){
	return x < y;
}

/*
 * Calculates if one integer is larger than another.
 *
 * x: Integer value to compare with y.
 * y: Integer value which x is compared to.
 *
 * returns 1 if true and 0 if false.
 */
int greatherThan(int x, int y){
	return x > y;
}

/*
 * Calculates the maximum or minimum integer for the three loops - for loop, while loop, and do while loop - and prints the result along with the actual maximum or minimum integer (the one being a system constant and not found in loops).
 *
 * posNeg: Integer value being either +1 or -1 depending on if there should be calculated a maxminum or a minimum integer respectively.
 */
void calculateMaxMinInt(int posNeg){
	// Sets variables corresponding to max or min depending on the variable posNeg
	char maxMin[10];
	char intMaxMinName[10];
	int intMaxMin;
	if (posNeg == 1) {
		strcpy(maxMin, "max");
		strcpy(intMaxMinName, "INT_MAX");
		intMaxMin = INT_MAX;
		intOperation = greatherThan;
	} else if (posNeg == -1) {
		strcpy(maxMin, "min");
		strcpy(intMaxMinName, "INT_MIN");
		intMaxMin = INT_MIN;
		intOperation = lessThan;
	} else {
		printf("ERROR: posNeg shall be either +1 or -1 but were %i.\n", posNeg);
		exit(EXIT_FAILURE);
	}
	// Actual MAX_INT or MIN_INT
	printMaxMinInt(maxMin, intMaxMinName, intMaxMin);
	// Initialise integer variable
	int i;
	// In the following all the i's start from posNeg, since it cannot be known if it can count as high or low as i.e. INT_MAX-10 or INT_MIN+10 respectively.
	// for loop
	for (i = posNeg; intOperation(i+(posNeg*1), i); i = i + posNeg) {}
	printMaxMinInt(maxMin, "for loop", i);
	// while loop
	i = posNeg;
	while (intOperation(i+(posNeg*1), i)) {
		i = i + posNeg;
	}
	printMaxMinInt(maxMin, "while loop", i);
	// do while loop
	i = posNeg;
	do {
		i = i + posNeg;
	} while (intOperation(i+(posNeg*1), i));
	printMaxMinInt(maxMin, "do while loop", i);
}
/*
 * Prints "Machine epsilon for <type> for <loop> is <value>."
 *
 * type: String containing the type of number. Supported types are "float", "double", "long double" and the same with capital first letter only.
 * loop: String containing the loop type, e.g. "for loop".
 * floatValue: Float containing the value of the calculated machine epsilon as a float. Just write an arbitrary value if type is not "float" or "Float".
 * doubleValue: Double containing the value of the calculated machine epsilon as a double. Just write an arbitrary value if type is not "double" or "Double".
 * longDoubleValue: Long double containing the value of the calculated machine epsilon as a long double. Just write an arbitrary value if type is not "long double" or "Long double".
 */
void printMachineEpsilon(char type[], char loop[], float floatValue, double doubleValue, long double longDoubleValue){
	// Using string comparison strcmp() to compare the strings. Returns 0 strings are equal.
	if (strcmp(type, "float") == 0 || strcmp(type, "Float") == 0) {
		printf("Machine epsilon for floats for %s is: %g.\n", loop, floatValue);
	} else if (strcmp(type, "double") == 0 || strcmp(type, "Double") == 0) {
		printf("Machine epsilon for doubles for %s is: %lg.\n", loop, doubleValue);
	} else if (strcmp(type, "long double") == 0 || strcmp(type, "Long double") == 0) {
		printf("Machine epsilon for long doubles for %s is: %Lg.\n", loop, longDoubleValue);
	} else {
		printf("ERROR: The chosen type '%s' is not one of the three supported: float, double, long double, or the same with only first letter capital.\n", type);
		exit(EXIT_FAILURE);
	}
}

/*
 * Calculates the machine epsilon for the three loops - for loop, while loop, and do while loop - for floats, and prints the results and the actual system value.
 */
void machineEpsilonFloat(void){
	// System value
	printMachineEpsilon("float", "FLT_EPSILON", FLT_EPSILON, 0, 0);
	// Initialisation
	float machEps;
	// for loop
	for (machEps = 1.0f; 1 + (machEps/2) > 1; machEps /= 2) {}
	printMachineEpsilon("float", "for loop", machEps, 0, 0);
	// while loop
	machEps = 1.0f;
	while (1 + (machEps/2) > 1) {
		machEps /= 2;
	}
	printMachineEpsilon("float", "while loop", machEps, 0, 0);
	// do while loop
	machEps = 1.0f;
	do {
		machEps /= 2;
	} while (1 + (machEps/2) > 1);
	printMachineEpsilon("float", "do while loop", machEps, 0, 0);
}

/*
 * Calculates the machine epsilon for the three loops - for loop, while loop, and do while loop - for doubles, and prints the results and the actual system value.
 */
void machineEpsilonDouble(void){
	// System value
	printMachineEpsilon("double", "DBL_EPSILON", 0, DBL_EPSILON, 0);
	// Initialisation
	double machEps;
	// for loop
	for (machEps = 1; 1 + (machEps/2) > 1; machEps /= 2) {}
	printMachineEpsilon("double", "for loop", 0, machEps, 0);
	// while loop
	machEps = 1;
	while (1 + (machEps/2) > 1) {
		machEps  /= 2;
	}
	printMachineEpsilon("double", "while loop", 0, machEps, 0);
	// do while loop
	machEps = 1;
	do {
		machEps /=2;
	} while (1 + (machEps/2) > 1);
	printMachineEpsilon("double", "do while loop", 0, machEps, 0);
}

/*
 * Calculates the machine epsilon for the three loops - for loop, while loop, and do while loop - for long doubles, and prints the results and the actual system value.
 */
void machineEpsilonLongDouble(void){
	// System value
	printMachineEpsilon("long double", "LDBL_EPSILON", 0, 0, LDBL_EPSILON);
	// Initialisation
	long double machEps;
	// for loop
	for (machEps = 1.0L; 1 + (machEps/2) > 1; machEps /= 2) {}
	printMachineEpsilon("long double", "for loop", 0, 0, machEps);
	// while loop
	machEps = 1.0L;
	while (1 + (machEps/2) > 1) {
		machEps /= 2;
	}
	printMachineEpsilon("long double", "while loop", 0, 0, machEps);
	// do while loop
	machEps = 1.0L;
	do {
		machEps /= 2;
	} while (1 + (machEps/2) > 1);
	printMachineEpsilon("long double", "do while loop", 0, 0, machEps);
}

/*
 * Calculates the machine epsilon for the three loops - for loop, while loop, and do while loop - for both floats, doubles , and long doubles, and prints the results. These are compared to the actual system values.
 */
void machineEpsilon(void){
	// Floats
	printDashes();
	printf(" Float\n");
	machineEpsilonFloat();
	// Doubles
	printDashes();
	printf(" Double\n");
	machineEpsilonDouble();
	// Long doubles
	printDashes();
	printf(" Long double\n");
	machineEpsilonLongDouble();
}

/*
 * Calculates and prints Exercise 1 and its subexercises.
 */
void solutionExercise1(void){
	printExercise("Exercise 1");
	// Part i
	printExercise("i");
	calculateMaxMinInt(1);
	// Part ii
	printExercise("ii");
	calculateMaxMinInt(-1);
	// Parti iii
	printExercise("iii");
	machineEpsilon();
}

/*
 * Calculates and prints the float sum 1.0f + 1.0f/2 + 1.0f/3 + ... + 1.0f/max.
 *
 * max: Integer representing the number of iterations.
 */
void sumUpFloat(int max){
	float sum = 0;
	for (int i = 1; i <= max; i++) {
		sum += 1.0f/i;
	}
	printf("sumUpFloat(max=%i) resulted in %g.\n", max, sum);
}

/*
 * Calculates and prints the float sum 1.0f/max + 1.0f/(max-1) + 1.0f/(max-2) + ...  +1.0f.
 *
 * max: Integer representing the number of iterations.
 */
void sumDownFloat(int max){
	float sum = 0;
	for (int i = 0; i < max; i++) {
		sum += 1.0f/(max - i);
	}
	printf("sumDownFloat(max=%i) resulted in %g.\n", max, sum);
}

/*
 * Calculates and prints the double sum 1.0 + 1.0/2 + 1.0/3 + ... + 1.0/max.
 *
 * max: Integer representing the number of iterations.
 */
void sumUpDouble(int max){
	double sum = 0;
	for (int i = 1; i <= max; i++) {
		sum += 1.0/i;
	}
	printf("sumUpDouble(max=%i) resulted in %lg.\n", max, sum);
}

/*
 * Calculates and prints the double sum 1.0/max + 1.0/(max-1) + 1.0/(max-2) + ...  +1.0.
 *
 * max: Integer representing the number of iterations.
 */
void sumDownDouble(int max){
	double sum = 0;
	for (int i = 0; i < max; i++) {
		sum += 1.0/(max - i);
	}
	printf("sumDownDouble(max=%i) resulted in %lg.\n", max, sum);
}

/*
 * Calculates and prints Exercise 2 and its subexercises.
 */
void solutionExercise2(void){
	printExercise("Exercise 2");
	// Initialise max
	int divideBy = 2;
	int max = INT_MAX / divideBy;
	printf("max is set to INT_MAX/%i.\n", divideBy);
	// Part i
	printExercise("i");
	sumUpFloat(max);
	sumDownFloat(max);
	// Part ii
	printExercise("ii");
	printf("The difference between the two sums is due to the fact, that the significant digits .\n");
	// Part iii
	printExercise("iii");
	printf("Yes/No, the sum converge as a function of 'max'.\n");
	// Part iv
	printExercise("iv");
	sumUpDouble(max);
	sumDownDouble(max);
	printf("Explanation ... .\n");
}

/*
 * Calculates and prints Exercise 3.
 */
void solutionExercise3(void){
	printExercise("Exercise 3");
	testEqual(); // Command in precision.c
}

/*
 * Writes the integer digit as a string with capital first letter, if the digit is a single digit, i.e. [0, ..., 9], while it prints "not a digit" for all other integers.
 *
 * i: Integer to write out name for. Name of the digit only prints if i is in range [0, 9].
 */
void nameDigit(int i){
	switch (i) {
		case 0:
			printf("Zero\n");
			break;
		case 1:
			printf("One\n");
			break;
		case 2:
			printf("Two\n");
			break;
		case 3:
			printf("Three\n");
			break;
		case 4:
			printf("Four\n");
			break;
		case 5:
			printf("Five\n");
			break;
		case 6:
			printf("Six\n");
			break;
		case 7:
			printf("Seven\n");
			break;
		case 8:
			printf("Eight\n");
			break;
		case 9:
			printf("Nine\n");
			break;
		default:
			printf("%i is not a single digit between 0 and 9.\n", i);
	}
}

void solutionExercise4(void){
	printExercise("Exercise 4");
	printf("Now we test the function nameDigit with the integers from -1 to 10:\n");
	for (int i = -1; i < 11; i++) {
		printf("nameDigit for input digit %i is: ", i);
		nameDigit(i);
	}
}

/*
 * The main function.
 */
int main(void){
	solutionExercise1();
	solutionExercise2();
	solutionExercise3();
	solutionExercise4();
	return 0;
}
