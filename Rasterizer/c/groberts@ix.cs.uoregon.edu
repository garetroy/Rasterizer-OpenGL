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
Lerp(double start, double end, double progress)
{

	return start+(end-start)*progress;

}

class Triangle
{
  public:
	double         X[3];
	double         Y[3];
	unsigned char  color[3];
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

	double desiredy, desiredx, m, b;
	double progress;
	desiredy = coord[1].y;

	progress = (desiredy - coord[0].y) / (coord[2].y - coord[0].y);
	desiredx = Lerp(coord[0].x, coord[2].x, progress);

	temp->X[0]  = coord[2].x;
	temp->Y[0]  = coord[2].y;
	temp->X[1]  = coord[1].x;
	temp->Y[1]  = coord[1].y;
	temp->X[2]  = desiredx;
	temp->Y[2]  = desiredy;

	temp2->X[0] = coord[1].x;
	temp2->Y[0] = coord[1].y;
	temp2->X[1] = coord[0].x;
	temp2->Y[1] = coord[0].y;
	temp2->X[2] = desiredx;
	temp2->Y[2] = desiredy;

	temp->color[0]  = color[0];
	temp->color[1]  = color[1];
	temp->color[2]  = color[2];
	temp2->color[0] = color[0];
	temp2->color[1] = color[1];
	temp2->color[2] = color[2];

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
	}

	if(triangletype.irregular){
		if((coord[0].y > coord[1].y && coord[0].y < coord[2].y) || (coord[0].y > coord[2].y && coord[0].y < coord[1].y)){ //if 0 is middle y
			tip.x = coord[0].x;
			tip.y = coord[0].y;
			if(coord[1].y < coord[2].y){
				left.x  = coord[1].x;
				left.y 	= coord[1].y;
				right.x = coord[2].x;
				right.y = coord[2].y;
			} else {
				left.x  = coord[2].x;
				left.y  = coord[2].y;
				right.x = coord[1].x;
				right.y = coord[1].y;
			}
		}
		else if((coord[1].y > coord[0].y && coord[1].y < coord[2].y) || (coord[1].y > coord[2].y && coord[1].y < coord[0].y)){ //if 1 is middle y
			tip.x = coord[1].x;
			tip.y = coord[1].y;
			if(coord[0].y < coord[2].y){
				left.x  = coord[0].x;
				left.y 	= coord[0].y;
				right.x = coord[2].x;
				right.y = coord[2].y;
			} else {
				left.x  = coord[2].x;
				left.y  = coord[2].y;
				right.x = coord[0].x;
				right.y = coord[0].y;
			}
		}		
		else if((coord[2].y > coord[1].y && coord[2].y < coord[0].y) || (coord[2].y > coord[0].y && coord[2].y < coord[1].y)){ //if 2 is middle y
			tip.x = coord[2].x;
			tip.y = coord[2].y;
			if(coord[1].y < coord[0].y){
				left.x  = coord[1].x;
				left.y 	= coord[1].y;
				right.x = coord[0].x;
				right.y = coord[0].y;
			} else {
				left.x  = coord[0].x;
				left.y  = coord[0].y;
				right.x = coord[1].x;
				right.y = coord[1].y;
			}
		}
	}

	if(triangletype.flattop || triangletype.flatbottom){
		if(coord[0].y == coord[1].y){
			tip.x = coord[2].x;
			tip.y = coord[2].y;
			if(coord[0].x > coord[1].x){
				left.x  = coord[1].x;
				left.y  = coord[1].y;
				right.x = coord[0].x;
				right.y = coord[0].y;
			} else {
				left.x  = coord[0].x;
				left.y  = coord[0].y;
				right.x = coord[1].x;
				right.y = coord[1].y;
			}
		}

		if(coord[1].y == coord[2].y){
			tip.x = coord[0].x;
			tip.y = coord[0].y;
			if(coord[1].x > coord[2].x){
				left.x  = coord[2].x;
				left.y  = coord[2].y;
				right.x = coord[1].x;
				right.y = coord[1].y;
			} else {
				left.x  = coord[1].x;
				left.y  = coord[1].y;
				right.x = coord[2].x;
				right.y = coord[2].y;
			}
		}
		if(coord[0].y == coord[2].y){
			tip.x = coord[1].x;
			tip.y = coord[1].y;
			if(coord[0].x > coord[2].x){
				left.x  = coord[2].x;
				left.y  = coord[2].y;
				right.x = coord[0].x;
				right.y = coord[0].y;
			} else {
				left.x  = coord[0].x;
				left.y  = coord[0].y;
				right.x = coord[2].x;
				right.y = coord[2].y;
			}
		}
	}


  	coord[0] = left;
  	coord[1] = tip;
  	coord[2] = right;
#ifdef DEBUG
 	fprintf(stderr, "Coordinates now arranged as (%f,%f)(%f,%f)(%f,%f)\n", coord[0].x, coord[0].y, coord[1].x, coord[1].y, coord[2].x, coord[2].y);
#endif
  	coordset = true;
  	return;
}


class Screen
{
	public:
		unsigned char   *buffer;
		int             width, height;
		//unsigned char   GetPixel(int w, int h) {return buffer[3*(h*width+w)];}
		void            SetPixel(int w, int h, unsigned char p[3]);
		void            SetRow(double x1, double x2, double y, unsigned char p[3]);
		void   			SetTriangle(Triangle *t);
};

void
Screen::SetPixel(int w, int h, unsigned char p[3])
{
	if(w >= width || w < 0 || h >= height || h < 0)
		return;
	buffer[3*(h*width+w)    ] = p[0];
	buffer[3*(h*width+w) + 1] = p[1];
	buffer[3*(h*width+w) + 2] = p[2];
}

void
Screen::SetRow(double x1, double x2, double y, unsigned char p[3])
{

	x1 = ceil441(x1);
	x2 = floor441(x2);
	for(;x1 <= x2; x1++){
		SetPixel(x1,y,p);
  	}
}

void
Screen::SetTriangle(Triangle *t)
{

	double x1,x2,maxycoord,minycoord,progress = 0;
	if(t->triangletype.flatbottom){
    		minycoord = ceil441(t->coord[0].y);
    		maxycoord = floor441(t->coord[1].y);

   			for(; minycoord <= maxycoord; minycoord++){
    			progress = (minycoord - t->coord[0].y) / (t->coord[1].y - t->coord[0].y);
   				x1 = Lerp(t->coord[0].x, t->coord[1].x, progress);
   				x2 = Lerp(t->coord[2].x, t->coord[1].x, progress);
      			SetRow(x1,x2,minycoord,t->color); // Half the jump, need to proprtion from over the scanline
			}

	}

	if(t->triangletype.flattop){
			x1 = t->coord[0].x;
    		x2 = t->coord[2].x;
    		maxycoord = floor441(t->coord[0].y);
    		minycoord = ceil441(t->coord[1].y);

   			for(; maxycoord >= minycoord; maxycoord--){
    			progress = (t->coord[0].y - maxycoord) / (t->coord[0].y - t->coord[1].y);
   				x1 = Lerp(t->coord[0].x, t->coord[1].x, progress);
   				x2 = Lerp(t->coord[2].x, t->coord[1].x, progress);
     			SetRow(x1,x2,maxycoord,t->color);

			}
	}

}



std::vector<Triangle>
GetTriangles(void)
{
	vtkPolyDataReader *rdr = vtkPolyDataReader::New();
	rdr->SetFileName("proj1c_geometry.vtk");
	cerr << "Reading" << endl;
	rdr->Update();
	if (rdr->GetOutput()->GetNumberOfCells() == 0)
	{
		cerr << "Unable to open file!!" << endl;
		exit(EXIT_FAILURE);
	}
	vtkPolyData *pd = rdr->GetOutput();
	int numTris = pd->GetNumberOfCells();
	vtkPoints *pts = pd->GetPoints();
	vtkCellArray *cells = pd->GetPolys();
	vtkFloatArray *colors = (vtkFloatArray *) pd->GetPointData()->GetArray("color_nodal");
	float *color_ptr = colors->GetPointer(0);
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
		tris[idx].color[0] = (unsigned char) color_ptr[4*ptIds[0]+0];
		tris[idx].color[1] = (unsigned char) color_ptr[4*ptIds[0]+1];
		tris[idx].color[2] = (unsigned char) color_ptr[4*ptIds[0]+2];
	}
	cerr << "Done reading" << endl;

	return tris;
}

int main()
{
	vtkImageData *image = NewImage(1786, 1344);
	unsigned char *buffer = 
	(unsigned char *) image->GetScalarPointer(0,0,0);
	int npixels = 1786*1344;
	for (int i = 0 ; i < npixels*3 ; i++)
		buffer[i] = 0;

	std::vector<Triangle> triangles = GetTriangles();
	
	double x1,x2 = 0;

	Screen screen;
	screen.buffer = buffer;
	screen.width = 1786;
	screen.height = 1344;

	Triangle **splittriangles;
	fprintf(stdout, "Starting Write\n");
	for(int i = 0; i <= triangles.size(); i++){
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