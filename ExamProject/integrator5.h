#ifndef NUMERICAL_INTEGRATION_INTEGRATION_TRI_H
#define NUMERICAL_INTEGRATION_INTEGRATION_TRI_H

double adapt24 (  double func(double), double leftEndPt, double rightEndPt, double absAcc, double relAcc, double second_funcVal, int integrationLimit );
double adapt (    double func(double), double leftEndPt, double rightEndPt, double absAcc, double relAcc );
double adapt24_old ( double func(double), double leftEndPt, double rightEndPt, double absAcc, double relAcc, double second_funcVal, double third_funcVal, int numOfRecursions, double* integrationError  );
double adapt_old ( double func(double), double leftEndPt, double rightEndPt, double absAcc, double relAcc, double* integrationError );
double integrate( double func(double), double leftEndPt, double rightEndPt, double absAcc, double relAcc );
double integrate_old( double func(double), double leftEndPt, double rightEndPt, double absAcc, double relAcc, double*  integrationError );

#endif //NUMERICAL_INTEGRATION_INTEGRATION_TRI_H
