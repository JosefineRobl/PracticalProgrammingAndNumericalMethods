#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

// Print exercise title
void printExercise(char* exercise){
	printf("=============== Exercise %s ===============\n", exercise);
}

// Integration function for exercise A
double f(double x, void* params){
	double z = *(double*)params;
	double f = log(x*z)/sqrt(x);
	return f;
}

// Integration function for exercise B
double errFun(double x, void* params){
	double z = *(double*)params;
	double errFun = 2/sqrt(M_PI)*exp(-z*x*x);
	return errFun;
}

/*
// Struct definition for exercise C
typedef struct{
	double x;
	int n;
} bessel;

// Integration function for exercise C
double besselFun(double t, void* params){
	bessel b = *(bessel*) params;
	double x = (*b).x;
	int n = (*b).n;
	// Definition of function
	double besselFun = 1/M_PI * cos(n*t - x*sin(t));
	return besselFun;
}
*/

// Integration of error function
double errFunIntegration(double x){
	gsl_function F;
	F.function = &errFun;
	
	double z = 1.0; // Since the function is just 2/sqrt(pi)*exp(-x^2)
	F.params = (void*)&z;

	// Allocation of space for the integration
	int allocationLimit = 999;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(allocationLimit);

	// Initialize variables
	double a = 0, // Integration limit
	       acc = 1e-6, epsilon = 1e-6, // Acceptance/tolerance
	       result, error; // For collecting the data about the run

	// Perform rutine
	gsl_integration_qags(&F, a, x, acc, epsilon, allocationLimit, w, &result, &error);

	// Free allocated workspace
	gsl_integration_workspace_free(w);
	
	return result;
}

/*
// Integration of the bessel function
double besselFunIntegration(double x, int n){
	gsl_function F;
	F.function = &besselFun;
	
	F.params = &b;
	
	// Allocation of space for the integration
	int allocationLimit = 999;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(allocationLimit);
	
	// Initialize variables
	double a = 0, b = M_PI, // Integration limits
	       acc = 1e-6, epsilon = 1e-6, // Acceptance/tolerance
	       result, error; // For collecting the data about the run
	
	// Perform rutine
	gsl_integration_qags(&F, a, b, acc, epsilon, allocationLimit, w, &result, &error);
	
	return result;
}
*/

// Exercise A
void exerciseA(void){
	printExercise("A");
	
	gsl_function F;
	F.function = &f;
	
	double z = 1.0; // Since function is ln(x)/sqrt(x)
	F.params = (void*)&z;
	
	// Allocation of space for the integration
	int allocationLimit = 999;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(allocationLimit);
	
	// Initialize variables
	double a = 0, b = 1, // Integration limits
	       acc = 1e-6, epsilon = 1e-6, // Acceptance/tolerance
	       result, error; // For collecting the data about the run
	
	// Perfom rutine
	gsl_integration_qags(&F, a, b, acc, epsilon, allocationLimit, w, &result, &error);
	printf("∫_0^1 dx (ln(x)/√x) = %g.\n", result);
	
	// Free workspace
	gsl_integration_workspace_free(w);
}

// Exercise B
void exerciseB(void){
	printExercise("B");
	
	// Open file to write data to
	FILE* file = fopen("dataB.txt", "w");
	// Writing the data from the integration to the file
	for (double x = -2; x <= 2; x += 1./8) {
		fprintf(file, "%10g \t %10g\n", x, errFunIntegration(x));
	}
	// Close the file
	fclose(file);
	
	printf("The data can be seen in dataB.txt, where the rows are the different runs and the collumns are x and the function value (errFun(x)) respectively.\n");
	printf("The plot can be seen in figure errfun.png.\n");
}

/*
// Exercise C
void exerciseC(void){
	printExercise("C");

	// Open file to write data to
	FILE* file = fopen("dataC.txt", "w");
	// Writing the data from the integration to the file
	for (double x = -2; x <= 2; x += 1./8) {
		fprintf(file, "%10g \t %10g \t %10g \t %10g\n", x, besselFunIntegration(x, 0), besselFunIntegration(x, 1), besselFunIntegration(x, 2));
	}
	// Close the file
	fclose(file);
	
	printf("The data can be seen in dataC.txt, where the rows are the different runs and the collumns are x, bessel(x, 0), bessel(x, 1) and bessel(x,2) respectively.\n");
	printf("The plot can be seen in figure bessel.png.\n");
}
*/

// Main function
int main(void){
	exerciseA();
	printf("\n");
	exerciseB();
	/*
	printf("\n");
	exerciseC();
	*/
	printf("==========================================\n");
	return 0;
}
