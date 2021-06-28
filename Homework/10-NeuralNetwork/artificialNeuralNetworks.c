#include <gsl/gsl_vector.h> // Implements vectors
#include "artificialNeuralNetworks.h"
#include "quasiNewton.h"

/*
 * Allocates space for the artificial neural network.
 * 
 * n: Number of neurons in the artificial neural network.
 * f: Function.
 * dfdt: Derivative of function.
 * F: Integral of funtion.
 */
artificialNeuralNetwork* annAlloc(int n, double (*f)(double t), double (*dfdt)(double t), double (*F)(double t)){
	artificialNeuralNetwork* network = malloc(sizeof(artificialNeuralNetwork));
	network->n = n;
	network->f = f;
	network->dfdt = dfdt;
	network->F = F;
	network->params = gsl_vector_alloc(3*n);
	return network;
}

/*
 * Free memory allocated by artificial neural network.
 */
void annFree(artificialNeuralNetwork* network){
	gsl_vector_free(network->params);
	free(network);
}

/*
 *
 */
double annResponse(artificialNeuralNetwork* network, double x){
	double sum = 0;
	// Initialization of weights
	double a, b, w;
	// Set weights and the suum
	for (int i = 0; i < network->n; i++) {
		a = gsl_vector_get(network->params, 3*i + 0);
		b = gsl_vector_get(network->params, 3*i + 1);
		w = gsl_vector_get(network->params, 3*i + 2);
		sum += network->f((x - a)/b)*w;
	}
	return sum;
}

/*
 *
 */
double annDerivative(artificialNeuralNetwork* network, double x){
	double sum = 0; 
	// Initialization of weights
	double a, b, w;
	// Set weights and the sum
	for (int i = 0; i < network->n; i++) {
		a = gsl_vector_get (network->parameters, 3*i + 0); 
		b = gsl_vector_get (network->parameters, 3*i + 1); 
		w = gsl_vector_get (network->parameters, 3*i + 2); 
		sum += network->derivative((x - a)/b)*w;
	}
	return sum; 
}

/*
 *
 */
double annIntegral(artificialNeuralNetwork* network, double x){
	double sum = 0; 
	// Initialization of weights
	double a, b, w;
	// Set weights and the sum
	for (int i = 0; i < network->n; i++) {
		a = gsl_vector_get (network->parameters, 3*i + 0); 
		b = gsl_vector_get (network->parameters, 3*i + 1); 
		w = gsl_vector_get (network->parameters, 3*i + 2); 
		sum += network->integral((x - a)/b)*w;
	}
	return sum; 
}

/*
 * Train the artificial neural network on xs, ys.
 */
void annTrain(artificialNeuralNetwork* network, gsl_vector* xs, gsl_vector* ys){
	double costFunction(gsl_vector* params){
		gsl_vector_memcpy(network->params, params);
		double sum = 0;
		double xi, yi, fi;
		for (int i = 0; i < xs->size; i++) {
			xi = gsl_vector_get(xs, i);
			yi = gsl_vector_get(ys, i);
			fi = annResponse(network, xi);
			sum += pow(fi - yi, 2);
		}
		// Average
		return sum/(x->size);
	}
	gsl_vector* params = gsl_vector_alloc(network->params->size);
	gsl_vector_memcpy(params, network->params);
	quasiNewton(costFunction, network->params, 1e-3);
	gsl_vector_memcpy(params, network->params);
	gsl_vector_free(params);
}
