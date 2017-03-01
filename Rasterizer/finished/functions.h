#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#include <cmath>

double Lerp(double start, double end, double proportion);

double GetNorm(double *vector);

double Dot(double *vector1, double *vector2);

void Normalize(double *vector);

void CrossProduct(double *createdvector, double *vector1, double *vector2);

void SubtractVectors(double *createdvector, double *vector1, double *vector2);

double ceil441(double f);

double floor441(double f);

#endif