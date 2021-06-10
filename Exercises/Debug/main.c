#include"stdio.h"
#include"gsl/gsl_matrix.h"
/* THE ORIGINAL CODE CAN BE SEEN BELOW:
 * #include"gsl_matrix.h"
 * 
 * // To include GSL matrix, one shall use the h-file: gsl/gsl_matrix.h
 */

int print_half_00(gsl_matrix* m)
{
	double half = 1./2;
	/* THE ORGINAL CODE CAN BE SEEN BELOW:
	 * double half = 1/2;
	 * 
	 * THE ERROR WAS:
	 * // The error here is that 1/2 = 0 while 1./2 = 0.5 which was what the student ment to write.
	 */
	int status = printf( "half m_{00} = %g\n", gsl_matrix_get(m,0,0)*half );
	/* THE ORIGINAL CODE CAN BE SEEN BELOW:
	 * int *status = printf( "half m_{00} = %i\n", gsl_matrix_get(&m,0,0)*half );
	 * 
	 * THE ERRORS WERE:
	 * // Here the student used %i to represent the number, but the number is a double (thus %g or %f) thus the integer notation is wrong to use.
	 * // Furthermore one should neither have a pointer to the integer (*status).
	 * // There should not be a reference to the adress of m (&m).
	 */
	/* THE ORIGINAL CODE CAN BE SEEN BELOW:
	 * gsl_matrix_free(m);
	 * 
	 * THE ERROR was:
	 * // One should not free the GSL matrix here, since it is not defined in this scope.
	 */
	return status;
}

int main(void)
{
	gsl_matrix* m = gsl_matrix_alloc(1,1);
	/* THE ORIGINAL CODE CAN BE SEEN BELOW:
	 * gsl_matrix m = gsl_matrix_alloc(0,0);
	 * 
	 * THE ERRORS WERE:
	 * // One should have a pointer to the GSL matrix while allocating space.
	 * // Futhermore, one should allocate at least for a 1x1 matrix (gsl_matrix_alloc(1,1)).
	 */
	gsl_matrix_set(m,0,0,66);
	/* THE ORIGINAL CODE CAN BE SEEN BELOW:
	 * gsl_matrix_set(&m,0,0,66);
	 * 
	 * THE ERROR WAS:
	 * // Using the set function, one should just use the pointer (m) from above and not (&m).
	 */
	printf("half m_{00} (should be 33):\n");
	int status = print_half_00(m);
	/* THE ORIGINAL CODE CAN BE SEEN BELOW:
	 * int *status = print_half_00(&m);
	 * 
	 * THE ERRORS WERE:
	 * // Again, one should not make a pointer to the integer (*status).
	 * // One should though use the pointer to the matrix (m) printing, since the function requires a gsl_matrix*, i.e. GSL matix pointer.
	 */
	if(status<=0)
		printf("status=%d : SOMETHING WENT TERRIBLY WRONG (status<=0)\n",status);
	else
		printf("status=%d : everything went just fine (status>0)\n",status);
	/* THE ORIGINAL CODE CAN BE SEEN BELOW:
	 * if(status>0)
	 * 	printf("status=%g : SOMETHING WENT TERRIBLY WRONG (status>0)\n",*status);
	 * else
	 * 	printf("status=%g : everything went just fine (status=0)\n",*status);
	 * 
	 * THE ERRORS WERE:
	 * // Cases are inverted. printf returns an integer proportional to the number of characters printed, thus all is well for status>0, while something went wrong for status<0. Printing a signle digit also let the function return status>=1. An empty string printed is 0 characters, thus status==0, but that should not happen, when we are working with numbers; there should always be printed some number. Negative integers are also possible to get, but there clearly yields errors. So due to us wanting to always print a number, I would call it an werror to get status==0.
	 * // status is an integer, thus one should use %d (or %i) instead of %g, which is for doubles.
	 * // status>0 in the first print-statement and status=0 in the second print-statement shall change due to condition change in if-statement.
	 * // Furthermore, when refering to the value of status, one should not use the pointer (*status) but just the variable name (status).
	 */
	gsl_matrix_free(m);
	/* THE ORIGINAL CODE CAN BE SEEN BELOW:
	 * gsl_matrix_free(&m);
	 * 
	 * THE ERROR WAS:
	 * // When freeing the allocated space for a GSL matrix, one shall use the pointer to the matrix (m).
	 */
	return 0;
}
