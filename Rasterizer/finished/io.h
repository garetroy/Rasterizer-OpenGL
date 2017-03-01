#ifndef _IO_H_
#define _IO_H_

#include <vtkImageData.h>
#include <vtkPNGWriter.h>

vtkImageData * NewImage(int width, int height);

void WriteImage(vtkImageData *img, const char *filename);

#endif