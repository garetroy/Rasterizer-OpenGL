#ifndef _TRIANGLE_H_
#define _TRIANGLE_H_

#include <stdlib.h>
#include <string.h>

#include <iostream>

#include "matrix.h"
#include "functions.h"

typedef struct tritype
{
    bool flatbottom;
    bool flattop;
    bool irregular;

}tritype;

class Coordinates
{
    public:
        double x;
        double y;
        double z;
        double colors[3];
        double normals[3];
        double s;

        Coordinates(void);

};

class Triangle
{
  public:
    double         X[3];
    double         Y[3];
    double         Z[3];
    double         colors[3][3];
    double         normals[3][3];
    double         S[3];
    Coordinates    coord[3]; // This is arranged as [Left Coordinate, Tip Coordinate, Right Coordinate] (Irreg = Bottom, Tip, Top)
    tritype        triangletype;
    bool           coordset;

                   Triangle(void);
    void           Translate(Matrix totaltrans);
    void           Untranslate(void);
    void           Print(ostream &o);
    void           InitializeTriangle(void);
    Triangle**     SplitTriangle(void);
    void           Determinetype(void);
    void           SetCoordinates(void);

};

#endif
