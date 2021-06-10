#include <stdio.h>

/*
 * Prints the exercise as "=== Exercise <str> ===".
 *
 * str: String containing the exercise number or letter.
 */
void printExercise(char* str);

/*
 * Prints the exercise as "=== Exercise <str> ===" to the specified file.
 *
 * str: String containing the exercise number or letter.
 * fileStream: Pointer to the FILE object that identifies the stream to which the printing should occur.
 */
void printExerciseFprint(char* str, FILE *fileStream);

/*
 * Prints the subtext for the exercise as "--- <subtext> ---".
 *
 * subtext: String containing the subtext to the exercise.
 */
void printSubtext(char* subtext);
