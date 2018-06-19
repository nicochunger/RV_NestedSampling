// This is a C implementation of the trueanomaly function to increase speed.

#include <math.h>
#include <stdlib.h>
#include "trueanomaly.h"

float trueanomaly(float M, float ecc, int niterationmax, int tol)
{
    float E = M;
    float E0 = M;
    float nu;
    int i,j;

    int niteration = 0;

    while(abs(E-E0)>tol || niteration==0)
    {
        for(j=0; j<50; j++)
        {
            E0 = E;

            float ff = E - ecc*sin(E) - M;
            float dff = 1 - ecc*cos(E);

            // Use Newton method
            E = E0 - ff / dff;
        }
        niteration+=50;
        if(niteration>=niterationmax)
            return -1;
    }

    // Compute true anomaly
    nu = 2. * atan2(sqrt(1.+ecc)*sin(E/2.), sqrt(1.-ecc)*cos(E/2.));

    return nu;
}