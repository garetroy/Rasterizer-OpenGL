#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

using std::cerr;
using std::endl;

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

int main()
{
    std::cerr << "In main!" << endl;
    vtkImageData *image = NewImage(1024, 1350);
    unsigned char *buffer = 
    (unsigned char *) image->GetScalarPointer(0,0,0);
    int x = -1;
    int r,g,b = 0;
    for(int i = 0; i<1350; i++){
        if(i%50 == 0){x++;r = x/9;g = (x/3)%3;b = x%3;}
        for(int j=0;j<1024;j++){
            if(r == 2){buffer[3*(i*1024+j)]=255;}else if(r == 1){buffer[3*(i*1024+j)]=128;}else if(r == 0){buffer[3*(i*1024+j)]=0;};
            if(g == 2){buffer[3*(i*1024+j)+1]=255;}else if(g == 1){buffer[3*(i*1024+j)+1]=128;}else if(g == 0){buffer[3*(i*1024+j)+1]=0;};
            if(b == 2){buffer[3*(i*1024+j)+2]=255;}else if(b == 1){buffer[3*(i*1024+j)+2]=128;}else if(b == 0){buffer[3*(i*1024+j)+2]=0;};
        }
    }
    WriteImage(image, "strips");
}
