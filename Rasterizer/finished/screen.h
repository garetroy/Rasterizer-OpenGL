#ifndef _SCREEN_H_
#define _SCREEN_H_

#include "matrix.h"
#include "camera.h"
#include "triangle.h"
#include "lighting.h"
#include "functions.h"


class Screen
{
    public:

        unsigned char           *buffer;
        int                     width, height;
        double                  *zbuffer;
        double                  shade;
        double                  bound_colors[2][3]; //This contains left and right x colors
        double                  pixel_colors[3];
        double                  bound_z[2];
        double                  viewdirection[3];
        double                  bound_shading[2];
        Matrix                  totaltrans;
        LightingParameters      lp;
        Camera                  c;

                                Screen(unsigned char *_buffer, int _width, int _height, double positionx, double positiony);
        void                    CalculatePhongShading(Triangle *t);
        void                    SetBoundColors(const Triangle * const t, const double proportion);
        void                    SetBoundZ(const Triangle * const t, const double proportion);
        void                    SetBoundShading(const Triangle * const t, const double proportion);
        void                    SetPixelColors(const double proportion);
        int                     SetZCoord(const int w, const int h, const double proportion);
        void                    SetPixel(const int w, const int h);
        void                    SetShading(const double proportion);
        void                    SetRow(double x1, double x2, const double y);
        void                    SetTriangle(Triangle *t);
                                ~Screen(void){delete [] zbuffer;}
};

#endif