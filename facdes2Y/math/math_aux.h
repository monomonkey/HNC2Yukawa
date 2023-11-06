#ifndef FACDES2Y_MATH_DOT_H    /* This is an "include guard" */
#define FACDES2Y_MATH_DOT_H    /* prevents the file from being included twice. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

// We define some global parameters
extern int nrows, ncols;
extern double rho, dr;
extern double *r, *q, x[2];
extern double *U, *Up, *sigmaVec;

extern void closrel(double *gamma, int potentialID, int closureID, double *cFuncMatrix, double T, double alpha);

void pp(double dr, double *matrix1, double *matrix2, double *prod);
void ONg(double *gammaInput, double *gammaOutput, int potentialID, int closureID, double *cFuncMatrix, \
         double T, double Tfin, double alpha, double rmax);    
void Extrap(double *gammaOutput, double *gammaInput, double rho, double drho);
void interp(int m, int n, double *r1, double *gammaInput, double *gammaOutput);
void Pres(double *f, double dr, double *eta);
void FFTM(double *inputDataMatrix, double rmax, int isDirect);
double calint(double *f, double dr);
void intt(double *h, double dr, double *sft);
void FT(double *c, double *c1, double *rk, double dr);
void FFT(double *inputData, double rmax, int isDirect);
void sinft(double *y, int nmax);
void realft(double *data, int n, int isign, int nmax);
void four1(double *data, int nn, int isign, int nmax);

#endif /* FACDES2Y_MATH_DOT_H */
