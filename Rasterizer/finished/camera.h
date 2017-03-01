#ifndef _CAMERA_H_
#define _CAMERA_H_

#include "functions.h"
#include "matrix.h"

class Camera
{
  public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];
    double          frame[3][3]; //v1,v2,v3
    
    void            SetCameraFrame(void);
    Matrix          CameraTransform(void);
    Matrix          ViewTransform(void);
    Matrix          DeviceTransform(double width, double height);
};

Camera GetCamera(int frame, int nframes);

#endif