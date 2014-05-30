#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#define C 1
#define C_Roe 1

#define PI 3.14159265359
#define MAXI(x, y) ((x) > (y) ? (x) : (y))
#define MINI(x, y) ((x) < (y) ? (x) : (y))
#define SQRT(X) (pow(X,0.5))

double** AllocateDynamicArray( int nRows, int nCols);

void FreeDynamicArray(double** dArray);

void Initialize_w(input* inputData, double** Array, double val);

void Initialize_f(input* inputData, double** Array, double val);

void update_f(double**f, double** w);

double alpha(double** w, int Row, int Col, int sign);

double flux(double** f, double** w, int Row, int Col, int sign, string FA);

double flux_LF(double** f, double** w, int Row, int Col, int sign);

double flux_simple(double** f, double** w, int Row, int Col, int sign);

double flux_Roe(double** f, double** w, int Row, int Col, int sign);

double flux_HLL(double** f, double** w, int Row, int Col, int sign);

double flux_HLLC(double** f, double** w, int Row, int Col, int sign);

void MVMult(double** M, double* V);

double flux_Roe(double** f, double** w, int Row, int Col, int sign);

double tau(double** w, double val);

void solve(input* inputData, double** w, double** f, double h);

void PrintFlowVariables(input* inputData, double** rho, double** u, double** E, double** e, double** p, double** Mach, int count, double h);

void PrintArray(double** Array, int nRows, int nCols, double h, string fileName);

#endif
