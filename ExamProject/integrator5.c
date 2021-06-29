#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "integration_tri.h"

double adapt24 ( double func(double), double leftEndPt, double rightEndPt, double absAcc, double relAcc, double second_funcVal, int integrationLimit ){

    double first_funcVal    =   func( leftEndPt + 1*( rightEndPt - leftEndPt )/6 );     // The value of the function in the left,
    double third_funcVal    =   func( leftEndPt + 5*( rightEndPt - leftEndPt )/6 );     // and right end points, respectively


    double higherOrderEval  =   ( 3*first_funcVal +  2*second_funcVal + 3*third_funcVal) * ( rightEndPt - leftEndPt ) / 8;     // Q from Dmitris notes, the trapezoid quadrature
    double lowerOrderEval   =   (   first_funcVal +   second_funcVal +   third_funcVal) * ( rightEndPt - leftEndPt ) / 3;     // q from Dmitris notes, the rectangular quadrature

    double tolerance        =   absAcc + relAcc * fabs(higherOrderEval);
    double error            =   fabs( higherOrderEval - lowerOrderEval );



    //assert( integrationLimit  > 0);
    if ( error < tolerance || integrationLimit  == 0 ){
        return higherOrderEval;
    }
    else{
        //printf("%g %g %g\n", first_funcVal, second_funcVal, third_funcVal);
        double new_higherOrderEval_left   =  adapt24( func,          leftEndPt                                 , leftEndPt +   (rightEndPt - leftEndPt) / 3, absAcc / sqrt(3.), relAcc, first_funcVal,  integrationLimit - 1 );
        double new_higherOrderEval_mid    =  adapt24( func, leftEndPt +   (rightEndPt - leftEndPt) / 3, leftEndPt + 2*(rightEndPt - leftEndPt) / 3, absAcc / sqrt(3.), relAcc, second_funcVal, integrationLimit - 1 );
        double new_higherOrderEval_right  =  adapt24( func, leftEndPt + 2*(rightEndPt - leftEndPt) / 3,          rightEndPt                                , absAcc / sqrt(3.), relAcc, third_funcVal , integrationLimit - 1 );


        return new_higherOrderEval_left + new_higherOrderEval_mid + new_higherOrderEval_right ;
    }
}

double adapt ( double func(double), double leftEndPt, double rightEndPt, double absAcc, double relAcc ){
    double  second_funcVal      =   func( leftEndPt + 3*( rightEndPt - leftEndPt ) / 6 );

    int     integrationLimit    =   (int) 1e6;

    return adapt24 ( func , leftEndPt , rightEndPt, absAcc , relAcc , second_funcVal, integrationLimit ) ;
}
