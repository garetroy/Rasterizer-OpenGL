#include "camera.h"

void 
Camera::SetCameraFrame(void)
{
    SubtractVectors(frame[2],position,focus);
    Normalize(frame[2]);
    CrossProduct(frame[0],up,frame[2]);
    Normalize(frame[0]);
    CrossProduct(frame[1],frame[2],frame[0]);
    Normalize(frame[1]);
}


Matrix 
Camera::CameraTransform(void)
{  
    Matrix cameramatrix;
    int i;
    double t[3];
    double empty[3] = {0,0,0};

    SetCameraFrame();
    SubtractVectors(t,empty,position);

    for(i = 0; i < 3; i++){
        cameramatrix.A[i][3] = 0;
        cameramatrix.A[3][i] = Dot(frame[i],t);
    }

    for(i = 0; i < 3; i++)
        for(int k = 0; k < 3; k++)
            cameramatrix.A[i][k] = frame[k][i];
        
    cameramatrix.A[3][3] = 1;

    return cameramatrix;

}

Matrix
Camera::ViewTransform(void)
{
    Matrix transdmatrix;

    for(int i = 0; i < 4; i++)
        for(int k = 0; k < 4; k++)
            transdmatrix.A[i][k] = 0;

    transdmatrix.A[0][0] = (1/std::tan(angle/2));
    transdmatrix.A[1][1] = (1/std::tan(angle/2));
    transdmatrix.A[2][2] = (far+near) / (far-near);
    transdmatrix.A[2][3] = -1;
    transdmatrix.A[3][2] = 2*far*near / (far-near);

    return transdmatrix;

}

Matrix
Camera::DeviceTransform(double width, double height)
{
    Matrix transdmatrix;

    for(int i = 0; i < 4; i++)
       for(int k = 0; k < 4; k++)
           transdmatrix.A[i][k] = 0;

    transdmatrix.A[0][0] = width/2;
    transdmatrix.A[1][1] = height/2;
    transdmatrix.A[2][2] = 1;
    transdmatrix.A[3][0] = width/2;
    transdmatrix.A[3][1] = height/2;
    transdmatrix.A[3][3] = 1;

    return transdmatrix;

}


double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera
GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;

    return c;
}