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
 * Function for the artificial neural network to hopefully learn.
 */
double function(double x){
	return cos(5*x - 1)*exp(-pow(x, 2)); 
}

/*
 * Exercise A.
 */
void exerciseA(void){
	// Number of neurons in the artificial neural network
	int n = 6;
	// Initialize the artificial neural network using the activation fucntion
	artificialNeuralNetwork* network = annAlloc(n, activation_function); 
	// Initialize the interval on the x-axis
	double xMin = -1, xMax = 1; 
	// Initialize number of points
	int m = 50;
	// Allocate memory for data points
	gsl_vector* x = gsl_vector_alloc(m)
	gsl_vector* y = gsl_vector_alloc(m); 
	// Generate data points
	for (int i = 0; i < m; i++) {
		gsl_vector_set (x, i, xMin + (xMax - xMin)*i/(m - 1)); 
		gsl_vector_set (y, i, function (gsl_vector_get (x, i))); 
	}
	// Setting the parameters for the network
	for (int i = 0; i < network->n; i++) {
		gsl_vector_set (network->parameters, 3*i + 0, a + (b - a)*i/(network->n - 1));
		gsl_vector_set (network->parameters, 3*i + 1, 1.); 
		gsl_vector_set (network->parameters, 3*i + 2, 1.); 
	}
	// Train the artificial neural network
	NN_train (network, x, y); 
	// Generate file with the generated points
	FILE* pointsFile = fopen("generatedPoints.txt", "w"); 
	for (int i = 0; i < m; i++) {
		fprintf (points, "%g\t%g\n", gsl_vector_get(x, i), gsl_vector_get(y, i));
	}
	fclose(pointsFile);
	// Generate file with the data from the functions
	FILE* data = fopen ("dataFunctions.txt", "w");
	for (double d = xMin; d < xMax; d += 1./64) {
		fprintf (data, "%g\t%g \n", d, annResponse(network, d)); 
	}
	fclose(data);
}

/*
 * Exercise B.
 */
void exerciseB(void){
}

/*
 * Main function.
 */
int main(void){
	exerciseA();
	exerciseB();
	return 0;
}
