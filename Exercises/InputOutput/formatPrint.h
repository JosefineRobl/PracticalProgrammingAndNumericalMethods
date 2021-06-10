#include <stdio.h>
#include <math.h>

/*
 * Prints the formatted header for the table.
 */
void formatHeaderXSinCos(void);

/*
 * Prints the formatted header for the table to a specific file.
 *
 * fileStream: Pointer to the FILE object, that identifies the stream to which the printing should occur.
 */
void formatHeaderXSinCosFprint(FILE *fileStream);

/*
 * Prints the formatted table entries.
 *
 * x: Double containing the x value.
 * sin: Double containing the sin(x) value.
 * cos: Double containing the cos(x) value.
 */
void formatTableXSinCos(double x, double sin, double cos);

/*
 * Prints the formatted table entries to a specific file.
 * 
 * x: Double containing the x value.
 * sin: Double containing the sin(x) value.
 * cos: Double containing the cos(x) value.
 * fileStream: Pointer to the FILE object, that identifies the stream to which the printing should occur.
 */
void formatTableXSinCosFprint(double x, double sin, double cos, FILE *fileStream);
