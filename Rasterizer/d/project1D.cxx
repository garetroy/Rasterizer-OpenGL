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

double
Lerp(double start, double end, double proportion)
{

	return start+(end-start)*proportion;

}

class Triangle
{
  public:
	double         X[3];
	double         Y[3];
	double 		   Z[3];
	double  	   colors[3][3];
	coordinates    coord[3]; // This is arranged as [Left Coordinate, Tip Coordinate, Right Coordinate] (irreg = Bottom, Tip Top)
	tritype   	   triangletype;
	bool           coordset;
				   Triangle(){triangletype.flatbottom, triangletype.flattop, triangletype.irregular, coordset = false;}
	void           InitializeTriangle();
	Triangle**	   SplitTriangle();
	void           Determinetype();
	void           SetCoordinates();
};

void
Triangle::InitializeTriangle()
{
	if((Y[0] == Y[1] && Y[1] == Y[2]) || (X[0] == X[1] && X[1] == X[2])){
		return;
	}
	if(Y[0] == Y[1] || Y[1] == Y[2] || Y[0] == Y[2]){
		Determinetype();
		SetCoordinates();

		return;
	}

	triangletype.irregular  = true;
	triangletype.flattop    = false;
	triangletype.flatbottom = false;
	SetCoordinates();

	return;
}

Triangle**
Triangle::SplitTriangle()
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

	double desiredy, desiredx, desiredz;
	double proportion;
	double desiredcolors[3];
	desiredy = coord[1].y;

	proportion = (desiredy - coord[0].y) / (coord[2].y - coord[0].y);
	desiredx = Lerp(coord[0].x, coord[2].x, proportion);

	desiredz = Lerp(coord[0].z,coord[2].z,proportion);

	for(int i = 0; i < 3; i ++){
		desiredcolors[i] = Lerp(coord[0].colors[i], coord[2].colors[i], proportion);
	}

	temp->X[0] = coord[2].x;
	temp->Y[0] = coord[2].y;
	temp->Z[0] = coord[2].z;
	memcpy(temp->colors[0],coord[2].colors,sizeof(double) * 3);
	temp->X[1] = coord[1].x;
	temp->Y[1] = coord[1].y;
	temp->Z[1] = coord[1].z;
	memcpy(temp->colors[1],coord[1].colors,sizeof(double) * 3);
	temp->X[2] = desiredx;
	temp->Y[2] = desiredy;
	temp->Z[2] = desiredz;
	memcpy(temp->colors[2],desiredcolors,sizeof(double) * 3);

	temp2->X[0] = coord[1].x;
	temp2->Y[0] = coord[1].y;
	temp2->Z[0] = coord[1].z;
	memcpy(temp2->colors[0],coord[1].colors,sizeof(double) * 3);
	temp2->X[1] = coord[0].x;
	temp2->Y[1] = coord[0].y;
	temp2->Z[1] = coord[0].z;
	memcpy(temp2->colors[1],coord[0].colors,sizeof(double) * 3);
	temp2->X[2] = desiredx;
	temp2->Y[2] = desiredy;
	temp2->Z[2] = desiredz;
	memcpy(temp2->colors[2],desiredcolors,sizeof(double) * 3);

	temparr[0] = temp;
	temparr[1] = temp2;

	return       temparr;
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
Triangle::SetCoordinates()
{
	coordinates left;
	coordinates tip;
	coordinates right;

	for(int i = 0; i <= 2; i++){
		coord[i].x = X[i];
		coord[i].y = Y[i];
		coord[i].z = Z[i];
		memcpy(coord[i].colors, colors[i], sizeof(double) * 3);
	}

	if(triangletype.irregular){
		if((coord[0].y > coord[1].y && coord[0].y < coord[2].y) || (coord[0].y > coord[2].y && coord[0].y < coord[1].y)){ //if 0 is middle y
			tip.x = coord[0].x;
			tip.y = coord[0].y;
			tip.z = coord[0].z;
			memcpy(tip.colors, coord[0].colors, sizeof(double) * 3);
			if(coord[1].y < coord[2].y){
				left.x  = coord[1].x;
				left.y 	= coord[1].y;
				left.z  = coord[1].z;
				memcpy(left.colors, coord[1].colors, sizeof(double) * 3);
				right.x = coord[2].x;
				right.y = coord[2].y;
				right.z = coord[2].z;
				memcpy(right.colors, coord[2].colors, sizeof(double) * 3);
			} else {
				left.x  = coord[2].x;
				left.y  = coord[2].y;
				left.z  = coord[2].z;
				memcpy(left.colors, coord[2].colors, sizeof(double) * 3);
				right.x = coord[1].x;
				right.y = coord[1].y;
				right.z = coord[1].z;
				memcpy(right.colors, coord[1].colors, sizeof(double) * 3);
			}
		}
		else if((coord[1].y > coord[0].y && coord[1].y < coord[2].y) || (coord[1].y > coord[2].y && coord[1].y < coord[0].y)){ //if 1 is middle y
			tip.x = coord[1].x;
			tip.y = coord[1].y;
			tip.z = coord[1].z;
			memcpy(tip.colors, coord[1].colors, sizeof(double) * 3);
			if(coord[0].y < coord[2].y){
				left.x  = coord[0].x;
				left.y 	= coord[0].y;
				left.z  = coord[0].z;
				memcpy(left.colors, coord[0].colors, sizeof(double) * 3);
				right.x = coord[2].x;
				right.y = coord[2].y;
				right.z = coord[2].z;
				memcpy(right.colors, coord[2].colors, sizeof(double) * 3);
			} else {
				left.x  = coord[2].x;
				left.y  = coord[2].y;
				left.z  = coord[2].z;
				memcpy(left.colors, coord[2].colors, sizeof(double) * 3);
				right.x = coord[0].x;
				right.y = coord[0].y;
				right.z = coord[0].z;
				memcpy(right.colors, coord[0].colors, sizeof(double) * 3);
			}
		}		
		else if((coord[2].y > coord[1].y && coord[2].y < coord[0].y) || (coord[2].y > coord[0].y && coord[2].y < coord[1].y)){ //if 2 is middle y
			tip.x = coord[2].x;
			tip.y = coord[2].y;
			tip.z = coord[2].z;
			memcpy(tip.colors, coord[2].colors, sizeof(double) * 3);
			if(coord[1].y < coord[0].y){
				left.x  = coord[1].x;
				left.y 	= coord[1].y;
				left.z  = coord[1].z;
				memcpy(left.colors, coord[1].colors, sizeof(double) * 3);
				right.x = coord[0].x;
				right.y = coord[0].y;
				right.z = coord[0].z;
				memcpy(right.colors, coord[0].colors, sizeof(double) * 3);
			} else {
				left.x  = coord[0].x;
				left.y  = coord[0].y;
				left.z  = coord[0].z;
				memcpy(left.colors, coord[0].colors, sizeof(double) * 3);
				right.x = coord[1].x;
				right.y = coord[1].y;
				right.z = coord[1].z;
				memcpy(right.colors, coord[1].colors, sizeof(double) * 3);
			}
		}
	}

	if(triangletype.flattop || triangletype.flatbottom){
		if(coord[0].y == coord[1].y){
			tip.x = coord[2].x;
			tip.y = coord[2].y;
			tip.z = coord[2].z;
			memcpy(tip.colors, coord[2].colors, sizeof(double) * 3);
			if(coord[0].x > coord[1].x){
				left.x  = coord[1].x;
				left.y  = coord[1].y;
				left.z  = coord[1].z;
				memcpy(left.colors, coord[1].colors, sizeof(double) * 3);
				right.x = coord[0].x;
				right.y = coord[0].y;
				right.z = coord[0].z;
				memcpy(right.colors, coord[0].colors, sizeof(double) * 3);
			} else {
				left.x  = coord[0].x;
				left.y  = coord[0].y;
				left.z  = coord[0].z;
				memcpy(left.colors, coord[0].colors, sizeof(double) * 3);
				right.x = coord[1].x;
				right.y = coord[1].y;
				right.z = coord[1].z;
				memcpy(right.colors, coord[1].colors, sizeof(double) * 3);
			}
		}

		if(coord[1].y == coord[2].y){
			tip.x = coord[0].x;
			tip.y = coord[0].y;
			tip.z = coord[0].z;
			memcpy(tip.colors, coord[0].colors, sizeof(double) * 3);
			if(coord[1].x > coord[2].x){
				left.x  = coord[2].x;
				left.y  = coord[2].y;
				left.z  = coord[2].z;
				memcpy(left.colors, coord[2].colors, sizeof(double) * 3);
				right.x = coord[1].x;
				right.y = coord[1].y;
				right.z = coord[1].z;
				memcpy(right.colors, coord[1].colors, sizeof(double) * 3);
			} else {
				left.x  = coord[1].x;
				left.y  = coord[1].y;
				left.z  = coord[1].z;
				memcpy(left.colors, coord[1].colors, sizeof(double) * 3);
				right.x = coord[2].x;
				right.y = coord[2].y;
				right.z = coord[2].z;
				memcpy(right.colors, coord[2].colors, sizeof(double) * 3);
			}
		}

		if(coord[0].y == coord[2].y){
			tip.x = coord[1].x;
			tip.y = coord[1].y;
			tip.z = coord[1].z;
			memcpy(tip.colors, coord[1].colors, sizeof(double) * 3);
			if(coord[0].x > coord[2].x){
				left.x  = coord[2].x;
				left.y  = coord[2].y;
				left.z  = coord[2].z;
				memcpy(left.colors, coord[2].colors, sizeof(double) * 3);
				right.x = coord[0].x;
				right.y = coord[0].y;
				right.z = coord[0].z;
				memcpy(right.colors, coord[0].colors, sizeof(double) * 3);
			} else {
				left.x  = coord[0].x;
				left.y  = coord[0].y;
				left.z = coord[0].z;
				memcpy(left.colors, coord[0].colors, sizeof(double) * 3);
				right.x = coord[2].x;
				right.y = coord[2].y;
				right.z = coord[2].z;
				memcpy(right.colors, coord[2].colors, sizeof(double) * 3);
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

		unsigned char   *buffer;
		int             width, height;
		double          *zbuffer;
		double 			bound_colors[2][3]; //This contains left and right x colors
		double 			pixel_colors[3];
		double 			bound_z[2];
		void			AllocateBuffer(){zbuffer = new double [height*width]; for(int i = 0; i < height*width; i++) zbuffer[i]=-1;}
		void            SetPixel(int w, int h);
		void            SetRow(double x1, double x2, double y);
		void   			SetTriangle(Triangle *t);
		void 			SetBoundColors(Triangle *t, double proportion);
		void  			SetBoundZ(Triangle *t, double proportion);
		int 			SetZCoord(int w, int h, double proportion);
		void 			SetPixelColors(double proportion);
						~Screen(){if(buffer != NULL) delete [] zbuffer;}
};

void
Screen::SetBoundColors(Triangle* t, double proportion)
{

	for(int i = 0; i < 3; i++ ){
			bound_colors[0][i] = Lerp(t->coord[0].colors[i],t->coord[1].colors[i],proportion);
			bound_colors[1][i] = Lerp(t->coord[2].colors[i],t->coord[1].colors[i],proportion);
	}

}

void
Screen::SetBoundZ(Triangle* t, double proportion)
{
	bound_z[0] = Lerp(t->coord[0].z,t->coord[1].z,proportion);
	bound_z[1] = Lerp(t->coord[2].z,t->coord[1].z,proportion);

}

void
Screen::SetPixelColors(double proportion)
{

	for(int i = 0; i < 3; i++)
		pixel_colors[i] = Lerp(bound_colors[0][i], bound_colors[1][i], proportion);

}

int
Screen::SetZCoord(int w, int h, double proportion)
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
Screen::SetPixel(int w, int h)
{

	if(w >= width || w < 0 || h >= height || h < 0)
		return;
	
	buffer[3*(h*width+w)    ] = ceil441(pixel_colors[0] * 255);
	buffer[3*(h*width+w) + 1] = ceil441(pixel_colors[1] * 255);
	buffer[3*(h*width+w) + 2] = ceil441(pixel_colors[2] * 255);

}

void
Screen::SetRow(double x1, double x2, double y)
{
	double prex1,prex2,startx,proportion= 0;

	prex2 = x2;
	startx = x1;
	x1 = ceil441(x1);
	x2 = floor441(x2);
	for(;x1 <=  x2; x1++){
		proportion = (x1 - startx) / (prex2 - startx);
		if(SetZCoord(x1,y,proportion) == 1){
			SetPixelColors(proportion);
			SetPixel(x1,y);
		}
	}
}

void
Screen::SetTriangle(Triangle *t)
{

	double x1,x2,maxycoord,minycoord,proportion,proportion2 = 0;

	if(t->triangletype.flatbottom){
    		minycoord = ceil441(t->coord[0].y);
    		maxycoord = floor441(t->coord[1].y);
   			for(; minycoord <= maxycoord; minycoord++){
    			proportion = (minycoord - t->coord[0].y) / (t->coord[1].y - t->coord[0].y);
   				x1 = Lerp(t->coord[0].x, t->coord[1].x, proportion);
   				x2 = Lerp(t->coord[2].x, t->coord[1].x, proportion);
   				SetBoundColors(t,proportion);
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
				SetBoundColors(t, proportion);
				SetBoundZ(t,proportion);
     			SetRow(x1,x2,maxycoord);

			}
	}
}



std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1d_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();
    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    float *colors_ptr = var->GetPointer(0);
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
        tris[idx].X[0] = pts->GetPoint(ptIds[0])[0];
        tris[idx].X[1] = pts->GetPoint(ptIds[1])[0];
        tris[idx].X[2] = pts->GetPoint(ptIds[2])[0];
        tris[idx].Y[0] = pts->GetPoint(ptIds[0])[1];
        tris[idx].Y[1] = pts->GetPoint(ptIds[1])[1];
        tris[idx].Y[2] = pts->GetPoint(ptIds[2])[1];
        tris[idx].Z[0] = pts->GetPoint(ptIds[0])[2];
        tris[idx].Z[1] = pts->GetPoint(ptIds[1])[2];
        tris[idx].Z[2] = pts->GetPoint(ptIds[2])[2];
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
            float val = colors_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate colors for " << val << endl;
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


int main()
{
	vtkImageData *image = NewImage(1000, 1000);
	unsigned char *buffer = 
	(unsigned char *) image->GetScalarPointer(0,0,0);
	int npixels = 1000*1000;
	for (int i = 0 ; i < npixels*3 ; i++)
			buffer[i] = 0;

	std::vector<Triangle> triangles = GetTriangles();

	/*std::vector<Triangle> triangles;
	Triangle t;
	Triangle t2;

  t.X[0] = 200;
  t.X[1] = 500;
  t.X[2] = 700;
  t.Y[0] = 600;
  t.Y[1] = 800;
  t.Y[2] = 200;
  t.Z[0] = 1;
  t.Z[1] = 1;
  t.Z[2] = 1;
  t.colors[0][0] = 1;
  t.colors[0][1] = 0;
  t.colors[0][2] = 0;
  t.colors[1][0] = 0;
  t.colors[1][1] = 0;
  t.colors[1][2] = 1;
  t.colors[2][0] = 0;
  t.colors[2][1] = 1;
  t.colors[2][2] = 0;
	triangles.push_back(t);*/
	
	double x1,x2 = 0;

	Screen screen;
	screen.buffer = buffer;
	screen.width = 1000;
	screen.height = 1000;
	screen.AllocateBuffer();

	Triangle **splittriangles;
	fprintf(stdout, "Starting Write\n");
	for(int i = 0; i < triangles.size(); i++){
	//for(int i = 2; i < 3; i++){
		if(i % 150000 == 0){
			fprintf(stdout, ".");
			fflush(stdout);
		}
		triangles[i].InitializeTriangle();
		if(triangles[i].triangletype.irregular){
			splittriangles = triangles[i].SplitTriangle();

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
   	fprintf(stdout,"\nDone Writting\n");
	WriteImage(image, "allTriangles");
}
