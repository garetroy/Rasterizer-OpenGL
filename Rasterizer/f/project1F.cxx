#include <iostream>
#include <stdlib.h>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <iostream>
#include <cstdlib>
#include <climits>
#include <cmath>

using std::cerr;
using std::endl;

double ceil441(double f)
{
    return ceil(f-0.00001);
}

double floor441(double f)
{
    return floor(f+0.00001);
}

typedef struct coordinates
{
    double x;
    double y;
    double z;
    double colors[3];
    double normals[3];
    double s;

}coordinates;

typedef struct tritype
{
    bool flatbottom;
    bool flattop;
    bool irregular;

}tritype;

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename, int k)
{
    std::string full_filename = filename;
    std::string number        = std::to_string(k);
    full_filename += number;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

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
        dotprod += vector1[i] * vector2[i];
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


struct LightingParameters
{
    LightingParameters(void)
    {
         lightDir[0] = -0.6;
         lightDir[1] = 0;
         lightDir[2] = -0.8;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 5.3;
         alpha = 7.5;
    };
  

    double lightDir[3];  // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
};

class Matrix
{
  public:
    double          A[4][4];

    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

void
Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix
Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
    Matrix rv;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }

    return rv;
}

void
Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
    ptOut[0] = ptIn[0]*A[0][0]
             + ptIn[1]*A[1][0]
             + ptIn[2]*A[2][0]
             + ptIn[3]*A[3][0];
    ptOut[1] = ptIn[0]*A[0][1]
             + ptIn[1]*A[1][1]
             + ptIn[2]*A[2][1]
             + ptIn[3]*A[3][1];
    ptOut[2] = ptIn[0]*A[0][2]
             + ptIn[1]*A[1][2]
             + ptIn[2]*A[2][2]
             + ptIn[3]*A[3][2];
    ptOut[3] = ptIn[0]*A[0][3]
             + ptIn[1]*A[1][3]
             + ptIn[2]*A[2][3]
             + ptIn[3]*A[3][3];
}

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

class Triangle
{
  public:
    double         X[3];
    double         Y[3];
    double         Z[3];
    double         colors[3][3];
    double         normals[3][3];
    double         S[3];
    coordinates    coord[3]; // This is arranged as [Left Coordinate, Tip Coordinate, Right Coordinate] (Irreg = Bottom, Tip, Top)
    tritype        triangletype;
    bool           coordset;

                   Triangle(void){triangletype.flatbottom = false; triangletype.flattop = false; triangletype.irregular = false; coordset = false;}
    void           Translate(Matrix totaltrans);
    void           Untranslate(void);
    void           Print(ostream &o);
    int            InitializeTriangle(void);
    Triangle**     SplitTriangle(void);
    void           Determinetype(void);
    void           SetCoordinates(void);

};

void
Triangle::Translate(Matrix totaltrans)
{

    double temp[4] = {0,0,0,0};

    for(int i = 0; i < 3; i++){
        double tcoordinates[4] = {X[i],Y[i],Z[i],1.};
        totaltrans.TransformPoint(tcoordinates, temp);
        X[i] = temp[0]/temp[3];
        Y[i] = temp[1]/temp[3];
        Z[i] = temp[2]/temp[3];
    }
}

void
Triangle::Print(ostream &o)
{
        char str[256];
        sprintf(str, "coord[0]:(%f,%f,%f) Shading:%f\n coord[1]:(%f,%f,%f) Shading:%f\n coord[2]:(%f,%f,%f) Shading:%f\n", 
                      coord[0].x,coord[0].y,coord[0].z,coord[0].s,
                      coord[1].x,coord[1].y,coord[1].z,coord[1].s,
                      coord[2].x,coord[2].y,coord[2].z,coord[2].s);
        o << str;
}


int
Triangle::InitializeTriangle(void)
{
    if((Y[0] == Y[1] && Y[1] == Y[2]) || (X[0] == X[1] && X[1] == X[2])){
        return 0;
    }

    if(Y[0] == Y[1] || Y[1] == Y[2] || Y[0] == Y[2]){
        Determinetype();
        SetCoordinates();
        return 1;
    }

    triangletype.irregular  = true;
    triangletype.flattop    = false;
    triangletype.flatbottom = false;

    SetCoordinates();

    return 1;
}

Triangle**
Triangle::SplitTriangle(void)
{
    if(!coordset){
        SetCoordinates();
    }

    Triangle** temparr  = (Triangle **)calloc(2, sizeof(Triangle *));
    if(temparr == NULL)
        return NULL;

    Triangle *temp = (Triangle *)malloc(sizeof(Triangle));
    if(temp == NULL)
        return NULL;

    Triangle *temp2 = (Triangle *)malloc(sizeof(Triangle));
    if(temp == NULL)
        return NULL;

    double desiredy, desiredx, desiredz, desiredshading;
    double proportion;
    double desiredcolors[3];
    double desirednormals[3];
    desiredy = coord[1].y;

    proportion = (desiredy - coord[0].y) / (coord[2].y - coord[0].y);

    desiredx       = Lerp(coord[0].x, coord[2].x, proportion);
    desiredz       = Lerp(coord[0].z,coord[2].z,proportion);
    desiredshading = Lerp(coord[0].s,coord[2].s,proportion);

    for(int i = 0; i < 3; i ++){
        desiredcolors[i]   = Lerp(coord[0].colors[i], coord[2].colors[i], proportion);
        desirednormals[i]  = Lerp(coord[0].normals[i], coord[2].normals[i], proportion);
    }
    
    Normalize(desirednormals); 
    Normalize(coord[0].normals);
    Normalize(coord[1].normals);
    Normalize(coord[2].normals);

    temp->X[0] = coord[2].x;
    temp->Y[0] = coord[2].y;
    temp->Z[0] = coord[2].z;
    temp->S[0] = coord[2].s;
    memcpy(temp->colors[0],coord[2].colors,sizeof(double) * 3);
    memcpy(temp->normals[0],coord[2].normals,sizeof(double) * 3);
    temp->X[1] = coord[1].x;
    temp->Y[1] = coord[1].y;
    temp->Z[1] = coord[1].z;
    temp->S[1] = coord[1].s;
    memcpy(temp->colors[1],coord[1].colors,sizeof(double) * 3);
    memcpy(temp->normals[1],coord[1].normals,sizeof(double) * 3);
    temp->X[2] = desiredx;
    temp->Y[2] = desiredy;
    temp->Z[2] = desiredz;
    temp->S[2] = desiredshading;
    memcpy(temp->colors[2],desiredcolors,sizeof(double) * 3);
    memcpy(temp->normals[2],desirednormals,sizeof(double) * 3);

    temp2->X[0] = coord[1].x;
    temp2->Y[0] = coord[1].y;
    temp2->Z[0] = coord[1].z;
    temp2->S[0] = coord[1].s;
    memcpy(temp2->colors[0],coord[1].colors,sizeof(double) * 3);
    memcpy(temp2->normals[0],coord[1].normals,sizeof(double) * 3);
    temp2->X[1] = coord[0].x;
    temp2->Y[1] = coord[0].y;
    temp2->Z[1] = coord[0].z;
    temp2->S[1] = coord[0].s;
    memcpy(temp2->colors[1],coord[0].colors,sizeof(double) * 3);
    memcpy(temp2->normals[1],coord[0].normals,sizeof(double) * 3);
    temp2->X[2] = desiredx;
    temp2->Y[2] = desiredy;
    temp2->Z[2] = desiredz;
    temp2->S[2] = desiredshading;
    memcpy(temp2->colors[2],desiredcolors,sizeof(double) * 3);
    memcpy(temp2->normals[2],desirednormals,sizeof(double) * 3);

    temparr[0] = temp;
    temparr[1] = temp2;

    return temparr;
}

void
Triangle::Determinetype()
{
    if(Y[0] == Y[1]){
        if(Y[2] > Y[0]){
            triangletype.flatbottom = true;
        }else{
            triangletype.flattop = true;
        }
    }
    else if(Y[1] == Y[2]){
        if(Y[0] > Y[1]){
            triangletype.flatbottom = true;
        }else{
            triangletype.flattop = true;
        }
    }
    else if(Y[0] == Y[2]){
        if(Y[1] > Y[0]){
            triangletype.flatbottom = true;
        }else{
            triangletype.flattop = true;
        }
    }
    triangletype.irregular = false;
    return; 
}

void
Triangle::SetCoordinates(void)
{
    coordinates left;
    coordinates tip;
    coordinates right;

    for(int i = 0; i < 3; i++){
        coord[i].x = X[i]; 
        coord[i].y = Y[i];
        coord[i].z = Z[i];
        coord[i].s = S[i];
        memcpy(coord[i].colors, colors[i], sizeof(double) * 3);
        memcpy(coord[i].normals, normals[i], sizeof(double) * 3);
    }

    if(triangletype.irregular){
        if((coord[0].y > coord[1].y && coord[0].y < coord[2].y) || (coord[0].y > coord[2].y && coord[0].y < coord[1].y)){ //if 0 is middle y
            tip.x = coord[0].x;
            tip.y = coord[0].y;
            tip.z = coord[0].z;
            tip.s = coord[0].s;
            memcpy(tip.colors,coord[0].colors, sizeof(double) * 3);
            memcpy(tip.normals,coord[0].normals, sizeof(double) * 3);
            if(coord[1].y < coord[2].y){
                left.x = coord[1].x;
                left.y = coord[1].y;
                left.z = coord[1].z;
                left.s = coord[1].s;
                memcpy(left.colors,coord[1].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[1].normals, sizeof(double) * 3);
                right.x = coord[2].x;
                right.y = coord[2].y;
                right.z = coord[2].z;
                right.s = coord[2].s;
                memcpy(right.colors,coord[2].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[2].normals, sizeof(double) * 3);
            } else {
                left.x = coord[2].x;
                left.y = coord[2].y;
                left.z = coord[2].z;
                left.s = coord[2].s;
                memcpy(left.colors,coord[2].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[2].normals, sizeof(double) * 3);
                right.x = coord[1].x;
                right.y = coord[1].y;
                right.z = coord[1].z;
                right.s = coord[1].s;
                memcpy(right.colors,coord[1].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[1].normals, sizeof(double) * 3);
            }
        }
        else if((coord[1].y > coord[0].y && coord[1].y < coord[2].y) || (coord[1].y > coord[2].y && coord[1].y < coord[0].y)){ //if 1 is middle y
            tip.x = coord[1].x;
            tip.y = coord[1].y;
            tip.z = coord[1].z;
            tip.s = coord[1].s;
            memcpy(tip.colors,coord[1].colors, sizeof(double) * 3);
            memcpy(tip.normals,coord[1].normals, sizeof(double) * 3);
            if(coord[0].y < coord[2].y){
                left.x = coord[0].x;
                left.y = coord[0].y;
                left.z = coord[0].z;
                left.s = coord[0].s;
                memcpy(left.colors,coord[0].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[0].normals, sizeof(double) * 3);
                right.x = coord[2].x;
                right.y = coord[2].y;
                right.z = coord[2].z;
                right.s = coord[2].s;
                memcpy(right.colors,coord[2].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[2].normals, sizeof(double) * 3);
            } else {
                left.x = coord[2].x;
                left.y = coord[2].y;
                left.z = coord[2].z;
                left.s = coord[2].s;
                memcpy(left.colors,coord[2].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[2].normals, sizeof(double) * 3);
                right.x = coord[0].x;
                right.y = coord[0].y;
                right.z = coord[0].z;
                right.s = coord[0].s;
                memcpy(right.colors,coord[0].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[0].normals, sizeof(double) * 3);
            }
        }       
        else if((coord[2].y > coord[1].y && coord[2].y < coord[0].y) || (coord[2].y > coord[0].y && coord[2].y < coord[1].y)){ //if 2 is middle y
            tip.x = coord[2].x;
            tip.y = coord[2].y;
            tip.z = coord[2].z;
            tip.s = coord[2].s;
            memcpy(tip.colors,coord[2].colors, sizeof(double) * 3);
            memcpy(tip.normals,coord[2].normals, sizeof(double) * 3);
            if(coord[1].y < coord[0].y){
                left.x  = coord[1].x;
                left.y  = coord[1].y;
                left.z  = coord[1].z;
                left.s = coord[1].s;
                memcpy(left.colors,coord[1].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[1].normals, sizeof(double) * 3);
                right.x = coord[0].x;
                right.y = coord[0].y;
                right.z = coord[0].z;
                right.s = coord[0].s;
                memcpy(right.colors,coord[0].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[0].normals, sizeof(double) * 3);
            } else {
                left.x  = coord[0].x;
                left.y  = coord[0].y;
                left.z  = coord[0].z;
                left.s = coord[0].s;
                memcpy(left.colors,coord[0].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[0].normals, sizeof(double) * 3);
                right.x = coord[1].x;
                right.y = coord[1].y;
                right.z = coord[1].z;
                right.s = coord[1].s;
                memcpy(right.colors,coord[1].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[1].normals, sizeof(double) * 3);
            }
        }
    }

    else if(triangletype.flattop || triangletype.flatbottom){
        if(coord[0].y == coord[1].y){
            tip.x = coord[2].x;
            tip.y = coord[2].y;
            tip.z = coord[2].z;
            tip.s = coord[2].s;
            memcpy(tip.colors,coord[2].colors, sizeof(double) * 3);
            memcpy(tip.normals,coord[2].normals, sizeof(double) * 3);
            if(coord[0].x > coord[1].x){
                left.x = coord[1].x;
                left.y = coord[1].y;
                left.z = coord[1].z;
                left.s = coord[1].s;
                memcpy(left.colors,coord[1].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[1].normals, sizeof(double) * 3);
                right.x = coord[0].x;
                right.y = coord[0].y;
                right.z = coord[0].z;
                right.s = coord[0].s;
                memcpy(right.colors,coord[0].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[0].normals, sizeof(double) * 3);
            } else {
                left.x = coord[0].x;
                left.y = coord[0].y;
                left.z = coord[0].z;
                left.s = coord[0].s;
                memcpy(left.colors,coord[0].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[0].normals, sizeof(double) * 3);
                right.x = coord[1].x;
                right.y = coord[1].y;
                right.z = coord[1].z;
                right.s = coord[1].s;
                memcpy(right.colors,coord[1].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[1].normals, sizeof(double) * 3);
            }
        }
        else if(coord[1].y == coord[2].y){
            tip.x = coord[0].x;
            tip.y = coord[0].y;
            tip.z = coord[0].z;
            tip.s = coord[0].s;
            memcpy(tip.colors,coord[0].colors, sizeof(double) * 3);
            memcpy(tip.normals,coord[0].normals, sizeof(double) * 3);
            if(coord[1].x > coord[2].x){
                left.x = coord[2].x;
                left.y = coord[2].y;
                left.z = coord[2].z;
                left.s = coord[2].s;
                memcpy(left.colors,coord[2].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[2].normals, sizeof(double) * 3);
                right.x = coord[1].x;
                right.y = coord[1].y;
                right.z = coord[1].z;
                right.s = coord[1].s;
                memcpy(right.colors,coord[1].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[1].normals, sizeof(double) * 3);
            } else {
                left.x = coord[1].x;
                left.y = coord[1].y;
                left.z = coord[1].z;
                left.s = coord[1].s;
                memcpy(left.colors,coord[1].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[1].normals, sizeof(double) * 3);
                right.x = coord[2].x;
                right.y = coord[2].y;
                right.z = coord[2].z;
                right.s = coord[2].s;
                memcpy(right.colors,coord[2].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[2].normals, sizeof(double) * 3);
            }
        }
        else if(coord[0].y == coord[2].y){
            tip.x = coord[1].x;
            tip.y = coord[1].y;
            tip.z = coord[1].z;
            tip.s = coord[1].s;
            memcpy(tip.colors,coord[1].colors, sizeof(double) * 3);
            memcpy(tip.normals,coord[1].normals, sizeof(double) * 3);
            if(coord[0].x > coord[2].x){
                left.x = coord[2].x;
                left.y = coord[2].y;
                left.z = coord[2].z;
                left.s = coord[2].s;
                memcpy(left.colors,coord[2].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[2].normals, sizeof(double) * 3);
                right.x = coord[0].x;
                right.y = coord[0].y;
                right.z = coord[0].z;
                right.s = coord[0].s;
                memcpy(right.colors,coord[0].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[0].normals, sizeof(double) * 3);
            } else {
                left.x = coord[0].x;
                left.y = coord[0].y;
                left.z = coord[0].z;
                left.s = coord[0].s;
                memcpy(left.colors,coord[0].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[0].normals, sizeof(double) * 3);
                right.x = coord[2].x;
                right.y = coord[2].y;
                right.z = coord[2].z;
                right.s = coord[2].s;
                memcpy(right.colors,coord[2].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[2].normals, sizeof(double) * 3);
            }
        }
    }
    coord[0] = left;
    coord[1] = tip;
    coord[2] = right;
    coordset = true;

    return;
}


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
        void                    SetShading(const double proportion);
        void                    SetPixel(const int w, const int h);
        void                    SetRow(double x1, double x2, const double y);
        void                    SetTriangle(Triangle *t);
                                ~Screen(void){if(zbuffer != NULL) delete [] zbuffer;}
};


Screen::Screen(unsigned char *_buffer, int _width, int _height, double positionx, double positiony)
{

    width   = _width;
    height  = _height;

    buffer  = _buffer;
    for (int i = 0 ; i < width*height*3 ; i++)
            buffer[i] = 0;

    zbuffer = new double [height*width]; 
    for(int i = 0; i < height*width; i++) 
        zbuffer[i]=-1;

    Matrix camtrans;
    Matrix viewtrans;
    Matrix devicetrans;

    c = GetCamera(positionx,positiony);

    camtrans    = c.CameraTransform();
    viewtrans   = c.ViewTransform();
    devicetrans = c.DeviceTransform((double)width,(double)height);
    totaltrans  = Matrix::ComposeMatrices(camtrans,viewtrans);
    totaltrans  = Matrix::ComposeMatrices(totaltrans,devicetrans);

}

void
Screen::CalculatePhongShading(Triangle *t)
{


    for(int i = 0; i < 3; i++){
        double R[3] = {0.,0.,0.};
        double vertex[3] = {t->X[i],t->Y[i],t->Z[i]};

        SubtractVectors(viewdirection, vertex, c.position);
        Normalize(viewdirection);

        double multiple = 2. * Dot(t->normals[i],lp.lightDir);
        for(int k = 0; k < 3; k++)
            R[k] = multiple * t->normals[i][k] - lp.lightDir[k];
 
        double diffuse  = fabs(Dot(t->normals[i],lp.lightDir)) * lp.Kd ;
        double power    = pow(Dot(R,viewdirection),lp.alpha);
        double specular = fmax(0.,lp.Ks*power);
        if(power != power)
            specular = 0.;

        t->S[i]= (diffuse + lp.Ka + specular);
   }
}

void
Screen::SetBoundColors(const Triangle* const t, const double proportion)
{

    for(int i = 0; i < 3; i++ ){
            bound_colors[0][i] = Lerp(t->coord[0].colors[i],t->coord[1].colors[i],proportion);
            bound_colors[1][i] = Lerp(t->coord[2].colors[i],t->coord[1].colors[i],proportion);
    }
}

void
Screen::SetBoundZ(const Triangle* const t, const double proportion)
{
    bound_z[0] = Lerp(t->coord[0].z,t->coord[1].z,proportion);
    bound_z[1] = Lerp(t->coord[2].z,t->coord[1].z,proportion);
}

void
Screen::SetBoundShading(const Triangle* const t, const double proportion)
{
    bound_shading[0] = Lerp(t->coord[0].s,t->coord[1].s,proportion);
    bound_shading[1] = Lerp(t->coord[2].s,t->coord[1].s,proportion);
}

void
Screen::SetPixelColors(const double proportion)
{
    for(int i = 0; i < 3; i++)
        pixel_colors[i] = Lerp(bound_colors[0][i], bound_colors[1][i], proportion);
}

int
Screen::SetZCoord(const int w, const int h, const double proportion)
{

    if(w >= width || w < 0 || h >= height || h < 0)
        return 0;

    double lerpval = Lerp(bound_z[0],bound_z[1],proportion);

    if(zbuffer[h*width+w] > lerpval){
        return 0; //Do not set pixel
    }

    if(zbuffer[h*width+w] <= lerpval){
        zbuffer[h*width+w] = lerpval;       
        return 1; //Set pixel
    }

    return 0;
}

void
Screen::SetShading(const double proportion)
{
    shade = Lerp(bound_shading[0], bound_shading[1], proportion);
}


void
Screen::SetPixel(const int w, const int h)
{

    if(w >= width || w < 0 || h >= height || h < 0)
        return;

    buffer[3*(h*width+w)    ] = fmin(255.,ceil441(shade * pixel_colors[0] * 255));
    buffer[3*(h*width+w) + 1] = fmin(255.,ceil441(shade * pixel_colors[1] * 255));
    buffer[3*(h*width+w) + 2] = fmin(255.,ceil441(shade * pixel_colors[2] * 255));

}

void
Screen::SetRow(double x1, double x2, const double y)
{
    double prex1,startx1,startx2,proportion = 0;

    startx1 = x1;
    startx2 = x2;
    x1 = ceil441(x1);
    x2 = floor441(x2);
    for(;x1 <=  x2; x1++){
        proportion = (x1 - startx1) / (startx2 - startx1);
        if(SetZCoord(x1,y,proportion) == 1){
            SetPixelColors(proportion);
            SetShading(proportion);
            SetPixel(x1,y);
        }
    }
}

void
Screen::SetTriangle(Triangle *t)
{

    double x1,x2,maxycoord,minycoord,proportion = 0;

    if(t->triangletype.flatbottom){
            minycoord = ceil441(t->coord[0].y);
            maxycoord = floor441(t->coord[1].y);
            for(; minycoord <= maxycoord; minycoord++){
                proportion = (minycoord - t->coord[0].y) / (t->coord[1].y - t->coord[0].y);
                x1 = Lerp(t->coord[0].x, t->coord[1].x, proportion);
                x2 = Lerp(t->coord[2].x, t->coord[1].x, proportion);
                SetBoundColors(t,proportion);
                SetBoundShading(t,proportion);
                SetBoundZ(t,proportion);
                SetRow(x1,x2,minycoord); 
                
            }
    }

    if(t->triangletype.flattop){
            maxycoord = floor441(t->coord[0].y);
            minycoord = ceil441(t->coord[1].y);
            for(; maxycoord >= minycoord; maxycoord--){
                proportion = (t->coord[0].y - maxycoord) / (t->coord[0].y - t->coord[1].y);
                x1 = Lerp(t->coord[0].x, t->coord[1].x, proportion);
                x2 = Lerp(t->coord[2].x, t->coord[1].x, proportion);
                SetBoundColors(t,proportion);
                SetBoundShading(t,proportion);
                SetBoundZ(t,proportion);
                SetRow(x1,x2,maxycoord);
            }
    }
}



std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1e_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();
/*
vtkDataSetWriter *writer = vtkDataSetWriter::New();
writer->SetInput(pd);
writer->SetFileName("hrc.vtk");
writer->Write();
 */

    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];

        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 },
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 }
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}

int main(void)
{
    vtkImageData *image = NewImage(1000, 1000);
    unsigned char *buffer = 
    (unsigned char *) image->GetScalarPointer(0,0,0);
    int npixels = 1000*1000;
    for (int i = 0 ; i < npixels*3 ; i++)
            buffer[i] = 0;

    std::vector<Triangle> triangles = GetTriangles();

    int frames = 1;
    Triangle **splittriangles;
    fprintf(stdout, "Starting Write\n");
    for(int k = 0; k < frames; k++){
        Screen screen(buffer,1000,1000,k,1000);
        for(int i = 0; i < triangles.size(); i++){  
            if(i % 15000 == 0){
                fprintf(stdout, ".");
                fflush(stdout);
            }

            screen.CalculatePhongShading(&triangles[i]);
            triangles[i].Translate(screen.totaltrans);
            if(triangles[i].InitializeTriangle() == 0)
                fprintf(stderr, "Found a flat line\n");

            if(triangles[i].triangletype.irregular){
              splittriangles = triangles[i].SplitTriangle();

              if(splittriangles == NULL){
                    fprintf(stderr, "SplitTriangles() failed\n");
                    return -1;
              }

                splittriangles[0]->InitializeTriangle();
                if(splittriangles[0]->triangletype.irregular){
                    fprintf(stderr, "Double irregular");
                    return -1;
                }

                screen.SetTriangle(splittriangles[0]);
                free(splittriangles[0]);

                splittriangles[1]->InitializeTriangle();
                if(splittriangles[1]->triangletype.irregular){
                    fprintf(stderr, "Double irregular");
                    return -1;
                }

                screen.SetTriangle(splittriangles[1]);
                free(splittriangles[1]);
                free(splittriangles);

            } else {
                screen.SetTriangle(&triangles[i]);
            }
        }
        fprintf(stdout,"\n");
        WriteImage(image, "allTriangles",k);
    }
}
