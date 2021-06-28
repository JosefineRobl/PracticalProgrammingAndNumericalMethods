#include<stdio.h>
#include<math.h>
#include<float.h>

/*
 * Calculates the numeric gradient.
 */
void gradient(int d, double f(int d, double* x), double* x, double* grad);

/*
 * Different quasiNewton that from Homework 9, since I could not get GSL vectors and non-nested functions (due to WSL) to collaborate.
 *
 * F: Objective function.
 * x: Pointer to vector; on input is the starting point, while on exit is the approximation to the root.
 * epsilon: Double containing the accuracy goal; on exit the absolute value of the gradient should be less than epsilon.
 *
 * returns the number of steps.
 */
int quasiNewton(int d, double f(int, double*), double* x, double acc);
