#ifndef _MAIN_H
#define _MAIN_H

#include <iostream>

#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>

#include "matrix.h"
#include "camera.h"
#include "functions.h"
#include "io.h"
#include "triangle.h"
#include "screen.h"
#include "lighting.h"

using std::cerr;
using std::endl;

std::vector<Triangle>
GetTriangles(void);

#endif