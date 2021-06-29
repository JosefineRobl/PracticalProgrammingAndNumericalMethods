#include <math.h>
#include <stdio.h>
#include <gsl/gsl_integration.h>


#include "integration_tri.h"

void print_testResults(char* string, double integralVal, int numOfCalls){
    printf("\n%s  : %g\n", string, integralVal);
    printf("Function was called %i times...\n", numOfCalls);
}

void print_whichTest(char* string, double exactVal){
    printf("\n%s %g", string, exactVal);
    printf("\n---------------------------------");
}


int main(int argc, char* argv[]){
    printf("\n\nMy student number is 201708238, and  mod(38, 22) = 16\n");

    printf("\n###################################");
    printf("\n# ------------------------------- #");
    printf("\n# ------- EXAM PROJECT 16 ------- #");
    printf("\n# ------------------------------- #");
    printf("\n###################################\n");
    printf("Title:\n");
    printf("\n¤ Adaptive recursive integrator with subdivision into three subintervals.\n");
    printf("\n\nDescription:\n");
    printf("\n¤ Implement a (one-dimensional) adaptive recursive integrator (open or closed quadrature, at your choice) which at each iteration subdivides the interval not into two, but into three sub-intervals. Reuse points. Compare with the adaptive integrator from your homework.\n");

    double leftEndPt    =   0;
    double rightEndPt   =   1;
    double absAcc       =   1e-6;
    double relAcc       =   1e-6;

    int     numOfCalls          =   0;
    double  integrationError    =   0;


    printf("= Recursive adaptive tri-division integrator =");
    printf("\n============================================== \n");
    printf("\nTesting integrator on different integrals...\n");

    double exactVal     =   2.0/3.0;
    print_whichTest("∫_0^1 dx √(x) = 2/3 =", exactVal);


    double firstTestFunc ( double x ){
        numOfCalls++;
        return sqrt( x ) ;
    };

    double integralVal  =   integrate( firstTestFunc, leftEndPt, rightEndPt, absAcc, relAcc);

    print_testResults("Result of numerical integration", integralVal, numOfCalls);
    numOfCalls          =   0;
    double ptrDump = 0;
    printf("-- Comparing with the integrator from the homework exercise:");
    double integralVal_old  =   integrate_old( firstTestFunc, leftEndPt, rightEndPt, absAcc, relAcc, &ptrDump);
    print_testResults("Result of numerical integration", integralVal_old, numOfCalls);


    double secondTestFunc ( double x ){
        numOfCalls++;
        return 4 * sqrt( 1 - x*x );
    };


    numOfCalls = 0;
    integrationError = 0;
    integralVal = integrate( secondTestFunc, leftEndPt, rightEndPt, absAcc, relAcc);
    exactVal = M_PI;
    print_whichTest("∫_0^1 dx 4√(1-x²) = π =", exactVal);
    print_testResults("Result of numerical integration", integralVal, numOfCalls);
    numOfCalls          =   0;
    ptrDump = 0;
    printf("-- Comparing with the integrator from the homework exercise:");
    integralVal_old  =   integrate_old( secondTestFunc, leftEndPt, rightEndPt, absAcc, relAcc, &ptrDump);
    print_testResults("Result of numerical integration", integralVal_old, numOfCalls);


    double err = 0;

    double fifthTestFunc( double x ){
        numOfCalls++;
        return exp(-x*x);
    }
    exactVal = sqrt(M_PI);
    print_whichTest("∫_-inf^inf dx exp(-x²) = √π =", exactVal);
    numOfCalls = 0;
    integrationError = 0;
    integralVal = integrate(fifthTestFunc, -INFINITY, INFINITY, absAcc, relAcc);
    print_testResults("Result of numerical integration", integralVal, numOfCalls);
    numOfCalls  =   0;
    ptrDump = 0;
    printf("-- Comparing with the integrator from the homework exercise:");
    integralVal_old  =   integrate_old( fifthTestFunc, -INFINITY, INFINITY, absAcc, relAcc, &ptrDump);
    print_testResults("Result of numerical integration", integralVal_old, numOfCalls);



    double sixthTestFunc( double x ){
        numOfCalls++;
        return 1/(1 + x*x);
    }

    exactVal = M_PI/2;
    print_whichTest("∫_0^inf dx 1/(1+x²) = π/2 =", exactVal);
    numOfCalls = 0;
    integrationError = 0;
    integralVal = integrate( sixthTestFunc, 0, INFINITY, absAcc, relAcc);
    print_testResults("Result of numerical integration", integralVal, numOfCalls);
    numOfCalls =   0;
    ptrDump = 0;
    printf("-- Comparing with the integrator from the homework exercise:");
    integralVal_old  =   integrate_old( sixthTestFunc, 0, INFINITY, absAcc, relAcc, &ptrDump);
    print_testResults("Result of numerical integration", integralVal_old, numOfCalls);



    exactVal = M_PI/2;
    print_whichTest("∫-inf,0 dx  1/(1+x²) = π/2 =", exactVal);
    numOfCalls = 0;
    integrationError = 0;
    integralVal = integrate( sixthTestFunc, -INFINITY, 0, absAcc, relAcc);
    print_testResults("Result of numerical integration", integralVal, numOfCalls);
    numOfCalls =   0;
    ptrDump = 0;
    printf("-- Comparing with the integrator from the homework exercise:");
    integralVal_old  =   integrate_old( sixthTestFunc, -INFINITY, 0, absAcc, relAcc, &ptrDump);
    print_testResults("Result of numerical integration", integralVal_old, numOfCalls);


    printf("==============================================\n");
    printf("\n\nConclusion:\n");
    printf("\n¤ The tri-division variant uses less recursive calls in general, so performs better all in all.\n");


    return 0;
}
