#include <math.h>

#define FORALL(i) for (int i = 0; i < n; i++)

void stepper_rk12 (
		int n, 
		void f (double t, double y[], double dydt[]),
		double t, double y[], double h, 
		double step[], double dy[]
		)
{
	double k0[n]; 
	f(t, y, k0); 

	double y1[n]; 
	FORALL(i) y1[i] = y[i] + (0.5*h)*k0[i]; 

	double k1[n]; 
	f(t + 0.5*h, y1, k1); 

	FORALL(i) step[i] = y[i] + h*k1[i]; 
	FORALL(i) dy[i] = (k1[i] - k0[i])*h;
}

void driver_rk12 (
		int n, 
		void f (double t, double y[], double dydt[]),
		double a, double b, 
		double y[], double h, 
		double acc, double eps
		)
{
		double t = a; 

		while (t < b)
		{
			if (t + h > b)
			{
				h = b - t;
			}

			double steps[n], dy[n];
		       	stepper_rk12 (n, f, t, y, h, steps, dy);

			double sum = 0; 
			FORALL(i) sum += y[i]*y[i];
			double norm_y = sqrt(sum);

			sum = 0;
			FORALL(i) sum += dy[i]*dy[i];
			double error = sqrt(sum); 

			double tolerance = (acc + eps*norm_y)*sqrt(h/(b - a));

			if (error < tolerance)
			{
				t = t + h; 
				FORALL(i) y[i] = steps[i];
			}

			if (error > 0)
			{
				h *= 0.95*pow(tolerance/error, 0.25); 
			}
			else 
			{
				h *= 2;
			}
		}
}
