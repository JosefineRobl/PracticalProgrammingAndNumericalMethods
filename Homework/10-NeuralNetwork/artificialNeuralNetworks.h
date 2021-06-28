#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

/*
 * Include guard
 */
#ifndef haveANNh
#define haveANNh

/*
 * Structure containing the network.
 */
typedef struct{
	int n; // Number of neurons
	double (*f)(double t); // Activation function
	double (*dfdt)(double t); // Derivative of activation function 
	double (*F)(double t); // Antiderivative of activation function
	double* params; // Parameters to be tuned during training
} artificialNeuralNetwork;

/*
 * Allocates space for the artificial neural network.
 * 
 * n: Number of neurons in the artificial neural network.
 * f: Function.
 * dfdt: Derivative of function.
 * F: Integral of funtion.
 */
artificialNeuralNetwork* annAlloc(int n, double (*f)(double t), double (*dfdt)(double t), double (*F)(double t));

/*
 * Free memory allocated by artificial neural network.
 */
void annFree(artificialNeuralNetwork* network);

/*
 * Respons function.
 */
double annResponse(artificialNeuralNetwork* network, double x);

/*
 * Derivative of respons function. Not needed for tuning due to use of same parameters for both function, derivative and antiderivative, but still need it to approximate derivative of tabulated function.
 */
double annDerivative(artificialNeuralNetwork* network, double x);

/*
 *
 */
double annIntegral(artificialNeuralNetwork* network, double x);

/*
 * Train the artificial neural network on xs, ys.
 */
void annTrain(artificialNeuralNetwork* network, int nx, double* xs, double* ys);

/*
 * Ending definition inside include guard
 */
#endif
