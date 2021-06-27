void stepper_rk12 (
		int n, 
		void f (double t, double y[], double dydt[]),
		double t, double y[], double h, 
		double step[], double error[]
		);
void driver_rk12 (
		int n, 
		void f (double t, double y[], double dydt[]),
		double a, double b, 
		double y[], double h, 
		double acc, double eps
		);
