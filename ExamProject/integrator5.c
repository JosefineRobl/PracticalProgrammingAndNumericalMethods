#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "integrator5.h"

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

double adapt24_old ( double func(double), double leftEndPt, double rightEndPt, double absAcc, double relAcc, double second_funcVal, double third_funcVal, int numOfRecursions, double* integrationError  ){
    assert( numOfRecursions < 1e6 );

    double first_funcVal    =   func( leftEndPt + 1*( rightEndPt - leftEndPt )/6 );     // The value of the function in the left,
    double fourth_funcVal   =   func( leftEndPt + 5*( rightEndPt - leftEndPt )/6 );     // and right end points, respectively

    double higherOrderEval  =   ( 2*first_funcVal + second_funcVal + third_funcVal + fourth_funcVal*2 ) * ( rightEndPt - leftEndPt ) / 6;     // Q from Dmitris notes, the trapezoid quadrature
    double lowerOrderEval   =   (   first_funcVal + second_funcVal + third_funcVal + fourth_funcVal   ) * ( rightEndPt - leftEndPt ) / 4;     // q from Dmitris notes, the rectangular quadrature

    double tolerance        =   absAcc + relAcc * fabs(higherOrderEval);
    double error            =   fabs( higherOrderEval - lowerOrderEval );

    if( error < tolerance ){
        *integrationError += error;
        return higherOrderEval;
    }
    else{
        double new_higherOrderEval_left  = adapt24_old( func,           leftEndPt                     , (leftEndPt + rightEndPt) / 2, absAcc / sqrt(2.), relAcc, first_funcVal, second_funcVal, numOfRecursions + 1, integrationError ) ;
        double new_higherOrderEval_right = adapt24_old( func, ( leftEndPt + rightEndPt ) / 2 ,           rightEndPt                  , absAcc / sqrt(2.), relAcc, third_funcVal, fourth_funcVal, numOfRecursions + 1, integrationError ) ;

        return new_higherOrderEval_left + new_higherOrderEval_right ;
    }
}

double adapt_old ( double func(double), double leftEndPt, double rightEndPt, double absAcc, double relAcc, double* integrationError ){
    double  second_funcVal      =   func( leftEndPt + 2*( rightEndPt - leftEndPt ) / 6 );
    double  third_funcVal       =   func( leftEndPt + 4*( rightEndPt - leftEndPt ) / 6 );

    int     numOfRecursions     =   0;

    return adapt24_old ( func , leftEndPt , rightEndPt ,absAcc , relAcc , second_funcVal , third_funcVal , numOfRecursions, integrationError ) ;
}

double integrate( double func(double), double leftEndPt, double rightEndPt, double absAcc, double relAcc ){
    if (isinf(-leftEndPt)){
        if (isinf(rightEndPt)){
            double transformInterval (double t) {
                return func(t / (1 - t * t)) * (1 + t * t) / pow(1 - t * t, 2);
            }
            return adapt(transformInterval, -1, 1, absAcc, relAcc);
        }
        else{
            double transformInterval (double t) {
                return func(rightEndPt + t / (1 + t)) / pow(1 + t, 2);
            }
            return adapt(transformInterval, -1, 0, absAcc, relAcc);
        }
    }
    else if (isinf(rightEndPt)) {
        double transformInterval (double t) {
            return func(leftEndPt + t / (1 - t)) / pow(1 - t, 2);
        }
        return adapt(transformInterval, 0, 1, absAcc, relAcc);
    }
    else{
        return adapt(func, leftEndPt, rightEndPt, absAcc, relAcc);
    }
}

double integrate_old( double func(double), double leftEndPt, double rightEndPt, double absAcc, double relAcc, double*  integrationError ){
    if (isinf(-leftEndPt)){
        if (isinf(rightEndPt)){
            double transformInterval (double t) {
                return func(t / (1 - t * t)) * (1 + t * t) / pow(1 - t * t, 2);
            }
            return adapt_old(transformInterval, -1, 1, absAcc, relAcc, integrationError);
        }
        else{
            double transformInterval (double t) {
                return func(rightEndPt + t / (1 + t)) / pow(1 + t, 2);
            }
            return adapt_old(transformInterval, -1, 0, absAcc, relAcc, integrationError);
        }
    }
    else if (isinf(rightEndPt)) {
        double transformInterval (double t) {
            return func(leftEndPt + t / (1 - t)) / pow(1 - t, 2);
        }
        return adapt_old(transformInterval, 0, 1, absAcc, relAcc, integrationError);
    }
    else{
        return adapt_old(func, leftEndPt, rightEndPt, absAcc, relAcc, integrationError);
    }
}
