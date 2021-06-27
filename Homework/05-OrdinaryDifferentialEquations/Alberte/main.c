#include <stdio.h>
#include "ode_rk12.h"

void sir_model (
		double time_between_contacts, 
		double t, double y[], double dydt[]
		)
{
	double population_size = 6*1e6;
	double recovery_time = 14;

	dydt[0] = -y[0]*y[1]/population_size/time_between_contacts;
	dydt[1] = y[0]*y[1]/population_size/time_between_contacts - y[1]/recovery_time;
	dydt[2] = y[1]/recovery_time;
}

static double time_between_contacts_file_var;

void sir_model_given_time_between_contacts (double t, double y[], double dydt[])
{
	sir_model (time_between_contacts_file_var, t, y, dydt);
}

int main () 
{
	double time_between_contacts[] = {1, 3, 5, 10}; 
	int n = sizeof (time_between_contacts)/sizeof (time_between_contacts[0]);

	int a = 0; 
	int b = 100; 
	double y[n]; 
	double h = 0.1; 
	double acc = 1e-2;
	double eps = 1e-2;

	FILE * susceptible = fopen ("susceptible.txt", "w"); 
	FILE * infectious = fopen ("infectious.txt", "w");
	FILE * removed = fopen ("removed.txt", "w"); 

	for (int t = 1; t < b; t++)
	{
		y[0] = 6*1e6;
		y[1] = 10;
		y[2] = 0; 

		fprintf (susceptible, "%i ", t);
		fprintf (infectious, "%i ", t);
		fprintf (removed, "%i ", t);

		for (int i = 0; i < n; i++)
		{
			time_between_contacts_file_var = time_between_contacts[i];

			driver_rk12 (n, sir_model_given_time_between_contacts, (double) a, (double) t, y, h, acc, eps);
			
			fprintf (susceptible, "%g ", y[0]);
			fprintf (infectious, "%g ", y[1]); 
			fprintf (removed, "%g ", y[2]); 
		}	

		fprintf (susceptible, "\n");
		fprintf (infectious, "\n");
		fprintf (removed, "\n");
	}

	fclose (susceptible); 
	fclose (infectious); 
	fclose (removed); 

	return 0; 
}
