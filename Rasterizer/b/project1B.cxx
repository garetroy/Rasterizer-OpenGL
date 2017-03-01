#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

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

typedef struct coordinates{
  double x;
  double y;
}coordinates;

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

class Triangle
{
  public:
      double         X[3];
      double         Y[3];
      unsigned char  color[3];
      coordinates    coord[3]; // This is arranged as [Left Coordinate, Top Coordinate, Right Coordinate]
      double         slopes[2]; //Arranged as [Left Slope, Right Slope]
      bool           coordset;
                     Triangle(){coordset = false;}
      void           SetCoordinates();
      void           SetSlopes();
  // would some methods for transforming the triangle in place be helpful?
};

void
Triangle::SetCoordinates(){
  coordinates left;
  coordinates top;
  coordinates right;
  int tempindex, tempindex2, tempindex3 = 0;
  double topy = 0;

  for(int i = 0; i < 3; i++){
    if(Y[i] > topy){
      tempindex = i;
      topy = Y[i];
    }
  }

  top.y = Y[tempindex];
  top.x = X[tempindex];

  switch(tempindex){
    case 0:
      tempindex2 = (X[1] < X[2] ? 1 : 2);
      left.x = X[tempindex2];
      left.y = Y[tempindex2];
      tempindex3 = (X[1] > X[2] ? 1 : 2);
      right.x = X[tempindex3];
      right.y = Y[tempindex3];
      break;
    case 1:
      tempindex2 = (X[0] < X[2] ? 0 : 2);
      left.x = X[tempindex2];
      left.y = Y[tempindex2];
      tempindex3 = (X[0] > X[2] ? 0 : 2);
      right.x = X[tempindex3];
      right.y = Y[tempindex3];
      break;
    case 2:
      tempindex2 = (X[0] < X[1] ? 0 : 1);
      left.x = X[tempindex2];
      left.y = Y[tempindex2];
      tempindex3 = (X[0] > X[1] ? 0 : 1);
      right.x = X[tempindex3];
      right.y = Y[tempindex3];
      break;
  } 
  coord[0] = left;
  coord[1] = top;
  coord[2] = right;
  return;
}

void
Triangle::SetSlopes(){
  if(!coordset)
    SetCoordinates();
  slopes[0] = ((coord[1].y - coord[0].y) != 0 ? ((coord[1].x - coord[0].x) / (coord[1].y - coord[0].y)) : 0);
  slopes[1] = ((coord[1].y - coord[2].y) != 0 ? ((coord[1].x - coord[2].x) / (coord[1].y - coord[2].y)) : 0);
}

class Screen
{
  public:
      unsigned char   *buffer;
      int             width, height;
      unsigned char   GetPixel(int w, int h) {return buffer[3*(h*width+w)];}
      void            SetPixel(int w, int h, unsigned char p[3]);
      void            SetRow(double x1, double x2, double y, unsigned char p[3]);

  // would some methods for accessing and setting pixels be helpful?
};

void
Screen::SetPixel(int w, int h, unsigned char p[3])
{
  if(w > width || w < 0 || h >= height || h < 0)
    return;
  buffer[3*(h*width+w)    ] = p[0];
  buffer[3*(h*width+w) + 1] = p[1];
  buffer[3*(h*width+w) + 2] = p[2];
}

void
Screen::SetRow(double x1, double x2, double y, unsigned char p[3]){
  x1 = ceil441(x1);
  x2 = floor441(x2);
  for(;x1 <= x2; x1++){
    SetPixel(x1,y,p);
  }
}

std::vector<Triangle>
GetTriangles(void)
{
   std::vector<Triangle> rv(100);

   unsigned char colors[6][3] = { {255,128,0}, {255, 0, 127}, {0,204,204}, 
                                  {76,153,0}, {255, 204, 204}, {204, 204, 0}};
   for (int i = 0 ; i < 100 ; i++)
   {
       int idxI = i%10;
       int posI = idxI*100;
       int idxJ = i/10;
       int posJ = idxJ*100;
       int firstPt = (i%3);
       rv[i].X[firstPt] = posI;
       if (i == 50)
           rv[i].X[firstPt] = -10;
       rv[i].Y[firstPt] = posJ;
       rv[i].X[(firstPt+1)%3] = posI+99;
       rv[i].Y[(firstPt+1)%3] = posJ;
       rv[i].X[(firstPt+2)%3] = posI+i;
       rv[i].Y[(firstPt+2)%3] = posJ+10*(idxJ+1);
       if (i == 95)
          rv[i].Y[(firstPt+2)%3] = 1050;
       rv[i].color[0] = colors[i%6][0];
       rv[i].color[1] = colors[i%6][1];
       rv[i].color[2] = colors[i%6][2];
   }

   return rv;
}



int main()
{
   vtkImageData *image = NewImage(1000, 1000);
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   int npixels = 1000*1000;
   for (int i = 0 ; i < npixels*3 ; i++)
       buffer[i] = 0;

   std::vector<Triangle> triangles = GetTriangles();

   Screen screen;
   screen.buffer = buffer;
   screen.width = 1000;
   screen.height = 1000;

   double x1, x2;
   unsigned char temp[3];
   temp[0] = 255;
   temp[1] = 0;
   temp[2] = 255;


   //YOUR CODE GOES HERE TO DEPOSIT TRIANGLES INTO PIXELS USING THE SCANLINE ALGORITHM
   //Find top right pixel
   //from top pixel to bottom pixel from coord, fill in color, also calculate hypotenuse, fill thoes pixels inbetween

  /*for(int i = 0; i < 100; i++){
    printf("\nTrangle #%d", i+1);
    printf("\nUnsorted Coord:");
    for(int k = 0; k < 3; k++)
      printf("(%f,%f)",triangles[i].X[k],triangles[i].Y[k]);
      triangles[i].SetCoordinates();
      triangles[i].SetSlopes();
    printf("\nLeftCoord:(%f,%f), TopCoord:(%f,%f), RightCoord:(%f,%f)\n", triangles[i].coord[0].x,triangles[i].coord[0].y,triangles[i].coord[1].x,triangles[i].coord[1].y,triangles[i].coord[2].x,triangles[i].coord[2].y);
    printf("\nLeft Slope:%f, Right Slope:%f\n", triangles[i].slopes[0], triangles[i].slopes[1]);
  }*/

  for(int i = 0; i < 100; i++){
    triangles[i].SetCoordinates();
    triangles[i].SetSlopes();
    x1 = triangles[i].coord[0].x;
    x2 = triangles[i].coord[2].x;
    for(double ycoord = ceil441(triangles[i].coord[0].y); ycoord <= floor441(triangles[i].coord[1].y); ycoord++){
      screen.SetRow(x1,x2,ycoord,triangles[i].color);
      x1 += triangles[i].slopes[0];
      x2 += triangles[i].slopes[1];
   }
  }


   WriteImage(image, "allTriangles");
}
