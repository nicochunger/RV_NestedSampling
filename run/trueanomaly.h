#ifndef TRUEANOMALY_H
#define TREUANOMALY_H

float trueanomaly(float M, float ecc, int niterationmax, int tol);
int trueanomaly_array(float* M, int n, float ecc, float* nu, int niterationmax, int tol);

#endif