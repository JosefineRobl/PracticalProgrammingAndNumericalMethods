#include <stdio.h>

#include "printExercise.h"

/*
 * Prints the exercise as "=== Exercise <str> ===".
 *
 * str: String containing the exercise number or letter.
 */
void printExercise(char* str){
	printf("=============== Exercise %s ===============\n", str);
}

/*
 * Prints the subtext for the exercise as "--- <subtext> ---".
 *
 * subtext: String containing the subtext to the exercise.
 */
void printSubtext(char* subtext){
	printf("---------- %s ----------\n", subtext);
}
