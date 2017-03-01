#include "functions.h"

double
Lerp(double start, double end, double proportion)
{

    return start+(end-start)*proportion;

}

double
GetNorm(double *vector)
{
 
    return sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);

}

double
Dot(double *vector1, double *vector2)
{

    double dotprod;

    dotprod = 0;

    for(int i = 0; i < 3; i++){
        dotprod   += vector1[i] * vector2[i];
    }

    return dotprod;

}

void
Normalize(double *vector)
{

    double norm = GetNorm(vector);

    for(int i = 0; i < 3; i++)
        vector[i] = vector[i] / norm;

}

void
CrossProduct(double *createdvector, double *vector1, double *vector2)
{
    createdvector[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1];
    createdvector[1] = vector1[2]*vector2[0] - vector1[0]*vector2[2];
    createdvector[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0];
}  

void
SubtractVectors(double *createdvector, double *vector1, double *vector2)
{
    createdvector[0] = vector1[0] - vector2[0];
    createdvector[1] = vector1[1] - vector2[1];
    createdvector[2] = vector1[2] - vector2[2];
}

double 
ceil441(double f)
{
    return ceil(f-0.00001);
}

double 
floor441(double f)
{
    return floor(f+0.00001);
}