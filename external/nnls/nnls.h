#ifndef NNLS_H
#define NNLS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "myMath.h"
typedef enum {FALSE = 0, TRUE = 1} Boolean;

void H12(int Mode, int LPivot, int L1, int M, double U[], int IUE, double* Up, double C[], int ICE, int ICV, int NCV);
void G1(double* A, double* B, double* COS, double* SIN, double* SIG);
void G2(double* COS, double* SIN, double* X, double* Y);
void nnls(double* A, int MDA, int M, int N, double* B, double* X, double* RNorm, double* W, double* ZZ, int* Index, int* Mode);
#endif
