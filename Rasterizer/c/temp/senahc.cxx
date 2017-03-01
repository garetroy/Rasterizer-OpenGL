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

//#define DEBUG;

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
	double    	   coord[5]; // for flat leftx, rightx, top/botx, miny, maxy. for irregular botx, topx, midy, miny, maxy 
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
	fprintf(stderr, "\n\n-------Init Triangle--------\n");
	if((Y[0] == Y[1] && Y[1] == Y[2]) || (X[0] == X[1] && X[1] == X[2])){
		fprintf(stderr, "Non-traingle Input\n");
		return;
	}
	if(Y[0] == Y[1] || Y[1] == Y[2] || Y[0] == Y[2]){
		fprintf(stderr, "It's a top or Bottom\n");
		Determinetype();
		SetCoordinates();

		return;
	}
	fprintf(stderr, "Irreular\n");
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
	desiredy = coord[3];

	//progress = (coord[4] - desiredy) / (coord[4] - desiredy);
	//desiredx = Lerp(coord[0], coord[4], progress);
	desiredy = coord[2];
	m = (coord[0] - coord[4])/(coord[0] - coord[4]);
	b = coord[4] - m*coord[2];
	desiredx = (desiredy - b)/m;


	temp->X[0]  = X[2];
	temp->Y[0]  = Y[2];
	temp->X[1]  = X[1];
	temp->Y[1]  = Y[1];
	temp->X[2]  = desiredx;
	temp->Y[2]  = desiredy;

	temp2->X[0] = X[1];
	temp2->Y[0] = Y[1];
	temp2->X[1] = X[0];
	temp2->Y[1] = Y[0];
	temp2->X[2] = desiredx;
	temp2->Y[2] = desiredy;

	temp->color[0]  = color[0];
	temp->color[1]  = color[1];
	temp->color[2]  = color[2];
	temp2->color[0] = color[0];
	temp2->color[1] = color[1];
	temp2->color[2] = color[2];

	//fprintf(stderr, "\n\n-----------------Split Triangle------------\n" );
	//fprintf(stderr, "Desired x: %f and y: %f\n", desiredx, desiredy );
	//fprintf(stderr, "Coordinates for Irreular Triangle:(%f,%f)(%f,%f)(%f,%f)\n", coord[0].x, coord[0].y, coord[1].x, coord[1].y, coord[2].x, coord[2].y);
	//fprintf(stderr, "Coordinates for Top Triangle:(%f,%f)(%f,%f)(%f,%f)\n", temp->X[0], temp->Y[0], temp->X[1], temp->Y[1], temp->X[2], temp->Y[2]);
	//fprintf(stderr, "Coordinates for Bottom Triangle:(%f,%f)(%f,%f)(%f,%f)\n\n", temp->X[0], temp->Y[0], temp->X[1], temp->Y[1], temp->X[2], temp->Y[2]);
	temparr[0] = temp;
	temparr[1] = temp2;

	return       temparr;
}

void
Triangle::Determinetype()
{
	fprintf(stderr, "\n----------------Determine type-----------\n");
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
	fprintf(stderr, "\n-----------Set Coordinates------------\n");fprintf(stderr, "irregular:%d, flatbottom:%d flattop:%d\n", triangletype.irregular, triangletype.flatbottom, triangletype.flattop);
	int i;
	double temp = 0;

	for(i = 0; i < 3; i ++){
		fprintf(stderr, "(%f,%f)", X[i], Y[i]);
	}

	if(triangletype.flattop || triangletype.flatbottom){
		temp = X[1];
		for(i = 0; i < 3; i++){
			if(X[i] < temp)
				temp = X[i];
		}
		coord[0] = temp;

		temp = X[1];
		for(i = 0; i < 3; i++){
			if(X[i] > temp)
				temp = X[i];
		}
		coord[1] = temp;

		temp = Y[1];
		for(i = 0; i < 3; i++){
			if(Y[i] < temp)
				temp = Y[i];
		}
		coord[3] = temp;

		temp = Y[1];
		for(i = 0; i < 3; i++){
			if(Y[i] > temp)
				temp = Y[i];
		}
		coord[4] = temp;

		if(triangletype.flattop){
			temp = 0;
			for(i = 0; i < 3; i++){
				if(Y[i] > coord[3] && Y[i] > coord[4])
					temp = i;
			}
			coord[2] = X[i];
		}

		if(triangletype.flatbottom){
			temp = 0;
			for(i = 0; i < 3; i++){
				if(Y[i] < coord[3] && Y[i] < coord[4])
					temp = i;
			}
			coord[2] = X[i];
		}
	}

	if(triangletype.irregular){
		temp = X[1];
		for(i = 0; i < 3; i++){
			if(X[i] < temp)
				temp = X[i];
		}
		coord[0] = temp;

		temp = X[1];
		for(i = 0; i < 3; i++){
			if(X[i] > temp)
				temp = X[i];
		}
		coord[1] = temp;

		temp = Y[1];
		for(i = 0; i < 3; i++){
			if(Y[i] < temp)
				temp = Y[i];
		}
		coord[3] = temp;

		temp = Y[1];
		for(i = 0; i < 3; i++){
			if(Y[i] > temp)
				temp = Y[i];
		}
		coord[4] = temp;

		for(i = 0; i < 3; i++){
			if(Y[i] < coord[4] && Y[i] > coord[3])
				coord[2] = Y[i];
		}
	}

	//fprintf(stderr, "Left X:%f, Right X:%f, Top/Bot(left/right) X:%f, miny:%f, minx:%f\n", coord[0],coord[1],coord[2],coord[3],coord[4] );

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
	//fprintf(stderr, "\n------------------------Set Row----------------\n");
	//fprintf(stderr, "Before Celing x1:%f, before flooringx2:%f\n",x1,x2 );
	x1 = ceil441(x1);
	x2 = floor441(x2);
	//fprintf(stderr, "After Celing x1:%f, After flooring x2:%f\n", x1,x2);
	for(;x1 <= x2; x1++){
		SetPixel(x1,y,p);
  	}
}

void
Screen::SetTriangle(Triangle *t)
{
	//fprintf(stderr, "\n------------------------Set Triangle----------------\n");
	//fprintf(stderr, "flatbottom:%d flattop:%d\n", t->triangletype.flatbottom, t->triangletype.flattop);
	//fprintf(stderr, "Coordinates for Triangle:(%f,%f)(%f,%f)(%f,%f)\n\n", t->coord[0].x, t->coord[0].y, t->coord[1].x, t->coord[1].y, t->coord[2].x, t->coord[2].y);
	double x1,x2,maxycoord,minycoord,progress = 0;
	if(t->triangletype.flatbottom){
    		minycoord = ceil441(t->coord[3]);
    		maxycoord = floor441(t->coord[4]);
    		//fprintf(stderr, "HRC: x1:%f, x2:%f\n", x1, x2);
    		//fprintf(stderr, "We are going from (ceil)minycoord:%f to (floor)maxycoord:%f, using ++\n", minycoord,maxycoord);
   			for(; minycoord <= maxycoord; minycoord++){
    			progress = (t->coord[4] - minycoord) / (t->coord[4] - t->coord[3]);
   				x1 = Lerp(t->coord[0], t->coord[3], progress);
   				x2 = Lerp(t->coord[1], t->coord[3], progress);
      			SetRow(x1,x2,minycoord,t->color); // Half the jump, need to proprtion from over the scanline
			}
		//fprintf(stderr, "HRC2 X1:%f X2:%f\n", x1,x2);
	}

	if(t->triangletype.flattop){
    		minycoord = ceil441(t->coord[3]);
    		maxycoord = floor441(t->coord[4]);
    		//fprintf(stderr, "We are going from (ceil)maxycoord:%f to (floor)minycoord:%f, using --\n", maxycoord,minycoord);
   			for(; maxycoord >= minycoord; maxycoord--){
    			progress = (t->coord[4] - minycoord) / (t->coord[4] - t->coord[3]);
   				x1 = Lerp(t->coord[0], t->coord[3], progress);
   				x2 = Lerp(t->coord[1], t->coord[3], progress);
     			SetRow(x1,x2,maxycoord,t->color);

			}
		//fprintf(stderr, "X1:%f X2:%f\n", x1,x2);
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

	/*Triangle temp3;

	temp3.X[0] = 400;
	temp3.Y[0] = 400;
	temp3.X[1] = 250;
	temp3.Y[1] = 600;
	temp3.X[2] = 400;
	temp3.Y[2] = 800;
	temp3.color[0] = 255;


	Triangle temp2;

	temp2.X[0] = 400;
	temp2.Y[0] = 400;
	temp2.X[1] = 600;
	temp2.Y[1] = 560;
	temp2.X[2] = 400;
	temp2.Y[2] = 800;
	temp2.color[1] = 255;

	Triangle temp[3];
	temp[0] = temp2;
	temp[1] = temp3;*/

	Triangle **splittriangles;
	fprintf(stdout, "Starting Write\n");
	//2515113
	for(int i = 0; i <= triangles.size(); i++){
		if(i % 150000 == 0){
			//fprintf(stdout, ".");
			//fflush(stdout);
		}
		triangles[i].InitializeTriangle();
		if(triangles[i].triangletype.irregular){
			splittriangles = triangles[i].SplitTriangle();

			splittriangles[0]->InitializeTriangle();
			if(splittriangles[0]->triangletype.irregular){
				//fprintf(stderr, "Double irregular");
				return -1;
			}

			screen.SetTriangle(splittriangles[0]);
			free(splittriangles[0]);

			splittriangles[1]->InitializeTriangle();
			if(splittriangles[1]->triangletype.irregular){
				//fprintf(stderr, "Double irregular");
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