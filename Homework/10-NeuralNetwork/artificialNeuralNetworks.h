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
	double (*f)(double);
	gsl_vector* params;
} artificialNeuralNetwork;

/*
 *
 */
artificialNeuralNetwork* annAlloc(int n, double (*f)(double));

/*
 *
 */
void annFree(artificialNeuralNetwork* network);

/*
 *
 */
double annResponse(artificialNeuralNetwork* network, double x);

/*
 *
 */
void annTrain(artificialNeuralNetwork* network, gsl_vector* xs, gsl_vector* ys);

/*
 * Ending definition inside include guard
 */
#endif
