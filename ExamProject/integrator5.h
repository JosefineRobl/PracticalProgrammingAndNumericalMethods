#ifndef NUMERICAL_INTEGRATION_INTEGRATION_TRI_H
#define NUMERICAL_INTEGRATION_INTEGRATION_TRI_H

double adapt24 (  double func(double), double leftEndPt, double rightEndPt, double absAcc, double relAcc, double second_funcVal, int integrationLimit );
double adapt (    double func(double), double leftEndPt, double rightEndPt, double absAcc, double relAcc );

#endif //NUMERICAL_INTEGRATION_INTEGRATION_TRI_H
