#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
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
	network->params = malloc(3*n*sizeof(double)); // Allocation of space for 3 parameters (the weights, a,b,w) for each neuron
	return network;
}

/*
 * Free memory allocated by artificial neural network.
 */
void annFree(artificialNeuralNetwork* network){
	free(network->params);
	free(network);
}

/*
 * Respons function.
 */
double annResponse(artificialNeuralNetwork* network, double x){
	double sum = 0;
	// Initialization of weights
	double a, b, w;
	// Set weights and the suum
	for (int i = 0; i < network->n; i++) {
		a = network->params[3*i];
		b = network->params[3*i + 1];
		w = network->params[3*i + 2];
		sum += network->f((x - a)/b)*w;
	}
	return sum;
}

/*
 * Derivative of respons function. Not needed for tuning due to use of same parameters for both function, derivative and antiderivative, but still need it to approximate derivative of tabulated function.
 */
double annDerivative(artificialNeuralNetwork* network, double x){
	double sum = 0; 
	// Initialization of weights
	double a, b, w;
	// Set weights and the sum
	for (int i = 0; i < network->n; i++) {
		a = network->params[3*i];
		b = network->params[3*i + 1];
		w = network->params[3*i + 2];
		sum += network->dfdt((x - a)/b)*w;
	}
	return sum; 
}

/*
 *
 */
double annIntegral(artificialNeuralNetwork* network, double x, double x0){
	double sum = 0; 
	// Initialization of weights
	double a, b, w;
	// Set weights and the sum
	for (int i = 0; i < network->n; i++) {
		a = network->params[3*i];
		b = network->params[3*i + 1];
		w = network->params[3*i + 2];
		sum += network->F((x - a)/b)*w*b - network->F((x0 - a)/b)*w*b;
	}
	return sum; 
}

/*
 * Cost function as non-nested function due to WSL.
 */
static int N;
static double* X;
static double* Y;
static artificialNeuralNetwork* NETWORK;
double costFunction(int d, double* p){
	assert(d == 3*NETWORK->n);
	for (int i = 0; i < d; i++) {
		NETWORK->params[i]=p[i];
	}
	double sum = 0;
	for(int i = 0; i < N; i++){
		double fi = ann_response(NETWORK, X[i]);
		sum += pow(fi-Y[i], 2);
	}
	return sum/N;
}

/*
 * Train the artificial neural network on xs, ys.
 */
void annTrain(artificialNeuralNetwork* network, int nx, double* xs, double* ys){
	// Set file-scope variables
	N = nx; X = xs; Y = ys;
	NETWORK = network;
	// Set accuracy
	double acc = 1e-3;
	// Define the size
	int d = 3*network->n;
	double p[d];
	for (int i = 0; i < d; i++) {
		p[i] = network->params[i];
	}
	quasiNewton(d, cost_func, p, acc);
	for (int i = 0; i < d; i++) {
		network->params[i] = p[i];
	}
}
