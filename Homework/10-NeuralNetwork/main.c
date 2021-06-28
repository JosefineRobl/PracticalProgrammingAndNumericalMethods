#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include "artificialNeuralNetworks.h"

/*
 * The activation function used for the neural network.
 */
double activationFunction(double x){
	return x*exp(-pow(x, 2));
}

/*
 * Derivative of activation function for the artificial neural network.
 */
double derivativeOfActivationFunction(double x){
	return exp(-pow(x, 2))*(1 - 2*pow(x, 2)*exp(-pow(x, 2)));
}

/*
 * Integrated activation function for the artificial neural network.
 */
double integralOfActivationFunction(double x){
	return - 1./2 * exp(-pow(x, 2));
}

/*
 * Function for the artificial neural network to hopefully learn.
 */
double function(double x){
	return sin(x)*exp(-x);
}

/*
 * Derivative of function for the artificial neural network to hopefully learn.
 */
double function(double x){
	return (cos(x) - sin(x))*exp(-x);
}

/*
 * Antiderivative of function for the artificial neural network to hopefully learn.
 */
double antiderivativeOfFunction(double x, double x0){
	return -1./2 * (sin(x) + cos(x))*exp(-x) + 1./2 * (sin(x0) + cos(x0)) * exp(-x0);
}

/*
 * Main function.
 */
int main(void){
	// Generate network
	int n=6; //number of neurons
	artificialNeuralNetwork* network = annAlloc(n, activationFunction, derivativeOfActivationFunction, integralOfActivationFunction);
	// Initialize the interval on the x-axis
	double xMin = -1, xMax = 1; 
	// Initialize number of points
	int m = 50;
	
	// Set parameters
	for(int i = 0; i < network->n; i++){
		network->params[3*i] = xMin + (xMax - xMin)*i/(network->n - 1);
		network->params[3*i + 1]=1;
		network->params[3*i + 2]=1;
	}
	// Generate data points
	double xs[m];
	double ys[m]; // Tabulated points
	double yms[m]; // Derivative
	double Ys[m]; // Antiderivative
	for(int i = 0; i < m; i++){
		xs[i] = xMin + (xMax - xMin)*i/(m - 1);
		ys[i] = function(xs[i]);
		yms[i] = derivativeOfFunction(xs[i]);
		Ys[i] = OfFunction(xs[i], xMin);
	}
	// Train the artificial neural network
	annTrain(network, x, y);
	// Print the found optimized patameters
	FILE* foundOptimizedParameters = fopen("foundOptimizedParameters.txt", "w");
	for (int i = 0; i < network->n; i++){
		double ai = network->params[3*i];
		double bi = network->params[3*i + 1];
		double wi = network->params[3*i + 2];
		fprintf(foundOptimizedParameters, "i = %d \t ai = %g \t bi = %g \t wi = %g\n", i, ai, bi, wi);
	}
	fclose(foundOptimizedParameters);
	// Generate file with the generated points
	FILE* pointsFile = fopen("generatedPoints.txt", "w");
	for (int i = 0; i < nx; i++) {
		fprintf(pointsFile, "%g\t%g\t%g\t%g\n", xs[i], ys[i], yms[i], Ys[i]);
	}
	fclose(pointsFile);
	// Generate file with the data from the functions
	FILE* data = fopen("dataFunctions.txt", "w");
	for (double z = xMin; z < xMax; z += 0.2) {
		fprintf(data, "%g\t%g\t%g\t%g\n", z, annResponse(network, z), annDerivative(network,z), annIntegral(network, z, xmin)); 
	}
	fclose(data);
	// Free network
	annFree(network);
	return 0;
}
