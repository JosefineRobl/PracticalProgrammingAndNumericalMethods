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
	return cos(5*x - 1)*exp(-pow(x, 2)); 
}

/*
 * Main function.
 */
int main(void){	
	// Number of neurons in the artificial neural network
	int n = 6;
	// Initialize the artificial neural network using the activation fucntion
	artificialNeuralNetwork* network = annAlloc(n, activationFunction, derivativeOfActivationFunction, integralOfActivationFunction); 
	// Initialize the interval on the x-axis
	double xMin = -1, xMax = 1; 
	// Initialize number of points
	int m = 50;
	// Allocate memory for data points
	gsl_vector* x = gsl_vector_alloc(m);
	gsl_vector* y = gsl_vector_alloc(m);
	gsl_vector* derivative = gsl_vector_alloc(m);
	gsl_vector* antiderivative = gsl_vector_alloc(m);
	// Generate data points
	for (int i = 0; i < m; i++) {
		gsl_vector_set(x, i, xMin + (xMax - xMin)*i/(m - 1));
		gsl_vector_set(y, i, function(gsl_vector_get(x, i))); // y = function(xi)
		gsl_vector_set(derivative, i, derivativeOfActivationFunction(gsl_vector_get(x, i))); // derivative(xi) = f'(xi)
		gsl_vector_set(antiderivative, i, integralOfActivationFunction(gsl_vector_get(x, i))); // antiderivative(xi) = F(xi)
	}
	// Setting the parameters for the network
	for (int i = 0; i < network->n; i++) {
		gsl_vector_set(network->params, 3*i + 0, xMin + (xMax - xMin)*i/(network->n - 1));
		gsl_vector_set(network->params, 3*i + 1, 1.); 
		gsl_vector_set(network->params, 3*i + 2, 1.); 
	}
	// Train the artificial neural network
	annTrain(network, x, y);
	// Print the found optimized patameters
	FILE* foundOptimizedParameters = fopen("foundOptimizedParameters.txt", "w");
	for (int i = 0; i < network->n; i++){
		double ai = gsl_vector_get(network->params, 3*i);
		double bi = gsl_vector_get(network->params, 3*i + 1);
		double wi = gsl_vector_get(network->params, 3*i + 2);
		fprintf(foundOptimizedParameters, "i = %i \t ai = %g \t bi = %g \t wi = %g\n", i, ai, bi, wi);
	}
	// Generate file with the generated points
	FILE* pointsFile = fopen("generatedPoints.txt", "w"); 
	for (int i = 0; i < m; i++) {
		fprintf(pointsFile, "%g\t%g\t%g\t%g\n", gsl_vector_get(x, i), gsl_vector_get(y, i), gsl_vector_get(derivative, i), gsl_vector_get(antiderivative, i));
	}
	fclose(pointsFile);
	// Generate file with the data from the functions
	FILE* data = fopen("dataFunctions.txt", "w");
	for (double d = xMin; d < xMax; d += 0.2) {
		fprintf(data, "%g\t%g\t%g\t%g\n", d, annResponse(network, d), annDerivative(network, d), annIntegral(network, d)); 
	}
	fclose(data);
	
	return 0;
}
