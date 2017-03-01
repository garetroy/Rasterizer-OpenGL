#include "screen.h"

Screen::Screen(unsigned char* _buffer, int _width, int _height, double positionx, double positiony)
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
Screen::SetBoundShading(const Triangle* const t, const double proportion)
{
    bound_shading[0] = Lerp(t->coord[0].s,t->coord[1].s,proportion);
    bound_shading[1] = Lerp(t->coord[2].s,t->coord[1].s,proportion);
}

void
Screen::SetBoundZ(const Triangle* const t, const double proportion)
{
    bound_z[0] = Lerp(t->coord[0].z,t->coord[1].z,proportion);
    bound_z[1] = Lerp(t->coord[2].z,t->coord[1].z,proportion);
}

void
Screen::SetPixelColors(const double proportion)
{
    for(int i = 0; i < 3; i++)
        pixel_colors[i] = Lerp(bound_colors[0][i], bound_colors[1][i], proportion);
}


void
Screen::SetShading(const double proportion)
{
    shade = Lerp(bound_shading[0], bound_shading[1], proportion);
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
    double prex1 = 0, startx1 = 0, startx2 = 0, proportion = 0;

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

    double x1 = 0, x2 = 0, maxycoord = 0, minycoord = 0, proportion = 0;

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
