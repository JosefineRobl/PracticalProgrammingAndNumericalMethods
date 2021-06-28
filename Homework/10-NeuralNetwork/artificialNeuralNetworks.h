#include <gsl/gsl_vector.h> // Implements vectors

/*
 * Include guard
 */
#ifndef haveANNh
#define haveANNh

/*
 *
 */
typedef struct{
	int n;
	double (*f)(double t);
	double (*dfdt)(double t);
	double (*F)(double t);
	gsl_vector* params;
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
 *
 */
double annResponse(artificialNeuralNetwork* network, double x);

/*
 *
 */
double annDerivative(artificialNeuralNetwork* network, double x);

/*
 *
 */
double annIntegral(artificialNeuralNetwork* network, double x);

/*
 * Train the artificial neural network on xs, ys.
 */
void annTrain(artificialNeuralNetwork* network, gsl_vector* xs, gsl_vector* ys);

/*
 * Ending definition inside include guard
 */
#endif
