#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"

#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkImageActor.h>
#include <vtkCallbackCommand.h>
#include <GLUT/glut.h>
#include <vector>
#include <cmath>

int earthday = 0;
int moonday = 0;
bool display = true;

/*
`
  For debugging and finding camera positions.

static void CameraModifiedCallback(vtkObject* caller,
                                   long unsigned int vtkNotUsed(eventId),
                                   void* vtkNotUsed(clientData),
                                   void* vtkNotUsed(callData) )
{
  std::cout << caller->GetClassName() << " modified" << std::endl;
 
  vtkCamera* camera = static_cast<vtkCamera*>(caller);
  // print the interesting stuff
  std::cout << "\tPosition: "
            << camera->GetPosition()[0] << ", "
            << camera->GetPosition()[1] << ", "
            << camera->GetPosition()[2] << std::endl;
  std::cout << "\tFocal point: "
            << camera->GetFocalPoint()[0] << ", "
            << camera->GetFocalPoint()[1] << ", "
            << camera->GetFocalPoint()[2] << std::endl;
}
*/

/*
  
  For debugging, this helps fine tune rotation

class customMouseInteractorStyle : public vtkInteractorStyleTrackballCamera
{
  public:
    static customMouseInteractorStyle* New();
    vtkTypeMacro(customMouseInteractorStyle, vtkInteractorStyleTrackballCamera);
 
    virtual void OnLeftButtonDown() 
    {
      earthday += 10;
      std::cout << "Pressed middle mouse button." << std::endl;
      // Forward events
      vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    }
 
    virtual void OnMiddleButtonDown() 
    {
      std::cout << "Pressed middle mouse button." << std::endl;
      // Forward events
      vtkInteractorStyleTrackballCamera::OnMiddleButtonDown();
    }
 
    virtual void OnRightButtonDown() 
    {
      std::cout << "Pressed right mouse button." << std::endl;
      earthday += 365;
      // Forward events
      vtkInteractorStyleTrackballCamera::OnRightButtonDown();
    }
 
};
 
vtkStandardNewMacro(customMouseInteractorStyle);
*/

class Triangle
{
  public:
      double         X[3];
      double         Y[3];
      double         Z[3];
};

class vtk441 : public vtkOpenGLPolyDataMapper
{
  /* This class is for setting up the scene for the mapper */
  protected:
   GLuint displayList;
   bool   initialized;
   float  size;

  public:
   vtk441(void)
   {
     initialized = false;
   }

  void IncrementSize(void)
   {
       size += 0.01;
       if (size > 2.0)
           size = 1.0;
   }
    
   void
   RemoveVTKOpenGLStateSideEffects(void)
   {
     float Info[4] = { 0, 0, 0, 1 };
     glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Info);
     float ambient[4] = { 1,1, 1, 1.0 };
     glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
     float diffuse[4] = { 1, 1, 1, 1.0 };
     glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
     float specular[4] = { 1, 1, 1, 1.0 };
     glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
   }


   void SetupLight(void)
   {
       glTranslatef(-1.0f, -1.0f, -1.0f);
       glEnable(GL_LIGHTING);
       glEnable(GL_LIGHT0);
       glEnable(GL_COLOR_MATERIAL);

       GLfloat diffuse0[4]  = { 0.6, 0.6, 0.6, 1 };
       GLfloat ambient0[4]  = { 0.2, 0.2, 0.2, 1 };
       GLfloat specular0[4] = { 0.3, 0.3, 0.3, 1 };
       GLfloat pos0[4]      = { 0, 0, 0, 1};
       glLightfv(GL_LIGHT0, GL_POSITION, pos0);
       glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse0);
       glLightfv(GL_LIGHT0, GL_AMBIENT, ambient0);
       glLightfv(GL_LIGHT0, GL_SPECULAR, specular0);
       float ambient[] = {0,0,0,0};
       glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, ambient);

       glDisable(GL_LIGHT2);
       glDisable(GL_LIGHT3);
       glDisable(GL_LIGHT5);
       glDisable(GL_LIGHT6);
       glDisable(GL_LIGHT7);
   }
};

class Planet
{ 
  /*This class is for each "planet" or spherical object in the sky
    All the variables not mentioned are self explanitory 
    Variables:
      textloc is the location of the texture file 
      locx and locy is the location (x,y) of the planet
      yearratio is the ratio of the planet's year in respect to earth year
      rotation is the ratio of the planets day in respect to earths day
  */

public:
  double        radius;
  double        distance;
  int           dims[3];
  unsigned char *buffer;
  char          textloc[1024];
  double        locx;
  double        locy;
  double        yearratio;
  double        rotation;
  vtkImageData  *img;
  vtkJPEGReader *rdr;
  GLuint        texture;

  Planet()
  {
    memset(textloc, 0, 1024 * sizeof(char));
    rdr = vtkJPEGReader::New();
  }

  void SetBuffer(void);
};

void
Planet::SetBuffer(void)
{
  rdr->SetFileName(textloc);
  rdr->Update();
  img = rdr->GetOutput();
  img->GetDimensions(dims);
  buffer = (unsigned char *) img->GetScalarPointer(0,0,0);
}

class vtk441Mapper : public vtk441
{
  /*
    This class is the mapper for the whole solar system
  */

 public:
   static vtk441Mapper *New();
   
   GLuint        displayList;
   bool          initialized;

   Planet solarsystem[11]; //Sun, Mercury, Venus, Earth, Moon, Mars ... etc

   vtk441Mapper()
   {
     initialized = false;

     solarsystem[0].radius  = 432.288;
     solarsystem[1].radius  = 1.516;
     solarsystem[2].radius  = 3.760;
     solarsystem[3].radius  = 3.959;
     solarsystem[4].radius  = 1.079;
     solarsystem[5].radius  = 2.106;
     solarsystem[6].radius  = 43.441;
     solarsystem[7].radius  = 36.184;
     solarsystem[8].radius  = 15.759;
     solarsystem[9].radius  = 15.299;
     solarsystem[10].radius = 0.7376;

     solarsystem[0].distance  = 0;
     solarsystem[1].distance  = 35.98;
     solarsystem[2].distance  = 67.24;
     solarsystem[3].distance  = 92.96;
     solarsystem[4].distance  = 0.2389;
     solarsystem[5].distance  = 141.6;
     solarsystem[6].distance  = 483.8;
     solarsystem[7].distance  = 888.2;
     solarsystem[8].distance  = 1784;
     solarsystem[9].distance  = 2795;
     solarsystem[10].distance = 3670;

     solarsystem[0].yearratio  = 0;
     solarsystem[1].yearratio  = 0.2409308693;
     solarsystem[2].yearratio  = 0.6160164271;
     solarsystem[3].yearratio  = 0.07392197125;
     solarsystem[4].yearratio  = 1;
     solarsystem[5].yearratio  = 1.8809034908;
     solarsystem[6].yearratio  = 11.9917864476;
     solarsystem[7].yearratio  = 30.1163586585;
     solarsystem[8].yearratio  = 84.8733744011;
     solarsystem[9].yearratio  = 164.8186173854;
     solarsystem[10].yearratio = 247.7;

     solarsystem[0].rotation  = 1.4166666667;
     solarsystem[1].rotation  = 58.6458333333;
     solarsystem[2].rotation  = 116.75;
     solarsystem[3].rotation  = 1;
     solarsystem[4].rotation  = 1;
     solarsystem[5].rotation  = 1.0277777778;
     solarsystem[6].rotation  = 0.02777777778;
     solarsystem[7].rotation  = 0.4458333333;
     solarsystem[8].rotation  = 0.7180555556;
     solarsystem[9].rotation  = 0.6708333333;
     solarsystem[10].rotation = 6.4;

     strcpy(solarsystem[0].textloc, "./textures/sun.jpg");
     strcpy(solarsystem[1].textloc, "./textures/mercury.jpg");
     strcpy(solarsystem[2].textloc, "./textures/venus.jpg");
     strcpy(solarsystem[3].textloc, "./textures/earth.jpg");
     strcpy(solarsystem[4].textloc, "./textures/moon.jpg");
     strcpy(solarsystem[5].textloc, "./textures/mars.jpg");
     strcpy(solarsystem[6].textloc, "./textures/jupiter.jpg");
     strcpy(solarsystem[7].textloc, "./textures/saturn.jpg");
     strcpy(solarsystem[8].textloc, "./textures/uranus.jpg");
     strcpy(solarsystem[9].textloc, "./textures/neptune.jpg");
     strcpy(solarsystem[10].textloc, "./textures/pluto.jpg");

     for(int i = 0; i < 11; i++)
      solarsystem[i].SetBuffer();

   }

   void DrawCylinder(void)
   {
       int nfacets = 30;
       glBegin(GL_TRIANGLES);
       for (int i = 0 ; i < nfacets ; i++)
       {
           double angle = 3.14159*2.0*i/nfacets;
           double nextAngle = (i == nfacets-1 ? 0 : 3.14159*2.0*(i+1)/nfacets);
           glNormal3f(0,0,1);
           glVertex3f(0, 0, 1);
           glVertex3f(cos(angle), sin(angle), 1);
           glVertex3f(cos(nextAngle), sin(nextAngle), 1);
           glNormal3f(0,0,-1);
           glVertex3f(0, 0, 0);
           glVertex3f(cos(angle), sin(angle), 0);
           glVertex3f(cos(nextAngle), sin(nextAngle), 0);
       }
       glEnd();
       glBegin(GL_QUADS);
       for (int i = 0 ; i < nfacets ; i++)
       {
           double angle = 3.14159*2.0*i/nfacets;
           double nextAngle = (i == nfacets-1 ? 0 : 3.14159*2.0*(i+1)/nfacets);
           glNormal3f(cos(angle), sin(angle), 0);
           glVertex3f(cos(angle), sin(angle), 1);
           glVertex3f(cos(angle), sin(angle), 0);
           glNormal3f(cos(nextAngle), sin(nextAngle), 0);
           glVertex3f(cos(nextAngle), sin(nextAngle), 0);
           glVertex3f(cos(nextAngle), sin(nextAngle), 1);
       }
       glEnd();
   }

   void SetUpTextures(Planet *p)
   {
    /*This sends each Planets texture down to the graphics card so that we don't
      have to keep sending them down (slowing the whole process down)
    */
    glGenTextures(1, &p->texture);
    glBindTexture(GL_TEXTURE_2D, p->texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, p->dims[0], p->dims[1], 0, GL_RGB,
                     GL_UNSIGNED_BYTE, p->buffer);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
   }

   void DrawSphere(Planet *p, bool sun = false, bool pluto = false)
   {
    /*
      This function set's up the textures for each planet, draws the shpeere, and applys the textures
      you can see that locx and locy are applying geometry to transform for orbits
      rotation around axis is done by Rotatef's
    */

    if(!initialized)
      SetUpTextures(p);

    GLUquadric *qobj = gluNewQuadric(); 

    gluQuadricTexture(qobj,GL_TRUE); 

    if(sun){
      gluQuadricNormals(qobj,GLU_SMOOTH);
      gluQuadricOrientation(qobj,GLU_OUTSIDE);
      gluQuadricDrawStyle(qobj, GLU_FILL);
    }

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,p->texture);

    if(!pluto)
      gluSphere(qobj,p->radius/10,50,50); 
    if(pluto)
      gluSphere(qobj,p->radius,50,50); 

    gluDeleteQuadric(qobj); 
   }

   void DrawSun(void)
   {
    glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    glRotatef((earthday/solarsystem[0].rotation),0,0,1);
    DrawSphere(&solarsystem[0],true);
    glPopMatrix();
   }

   void DrawMercury(void)
   {
    glPushMatrix();
    solarsystem[1].locx = (solarsystem[0].radius/10 + solarsystem[1].distance/10) * cos((earthday/solarsystem[1].yearratio) / 180.0 * 3.1415);
    solarsystem[1].locy = (solarsystem[0].radius/10 + solarsystem[1].distance/10) * sin((earthday/solarsystem[1].yearratio) / 180.0 * 3.1415);
    glTranslatef(solarsystem[1].locx, solarsystem[1].locy, 0);
    glRotatef((earthday/solarsystem[1].rotation),0,0,1);
    DrawSphere(&solarsystem[1]);
    glPopMatrix();
   }   

   void DrawVenus(void)
   {
    glPushMatrix();
    solarsystem[2].locx = (solarsystem[0].radius/10 + solarsystem[2].distance/10) * cos((earthday/solarsystem[2].yearratio) / 180.0 * 3.1415);
    solarsystem[2].locy = (solarsystem[0].radius/10 + solarsystem[2].distance/10) * sin((earthday/solarsystem[2].yearratio) / 180.0 * 3.1415);
    glTranslatef(solarsystem[2].locx, solarsystem[2].locy, 0);
    glRotatef((earthday/solarsystem[2].rotation),0,0,1);
    DrawSphere(&solarsystem[2]);
    glPopMatrix();
   }  

   void DrawMoon(void)
   {
    glPushMatrix();
    solarsystem[4].locx = (solarsystem[4].distance + solarsystem[3].radius - 2) * cos((moonday/solarsystem[3].yearratio) / 180.0 * 3.1415);
    solarsystem[4].locy = (solarsystem[4].distance + solarsystem[3].radius - 2) * sin((moonday/solarsystem[3].yearratio) / 180.0 * 3.1415);
    glTranslatef(solarsystem[4].locx, solarsystem[4].locy, 0);
    DrawSphere(&solarsystem[4],true);
    glPopMatrix();
   }

   void DrawEarth(void)
   {
    glPushMatrix();
    solarsystem[3].locx = (solarsystem[0].radius/10 + solarsystem[3].distance/10) * cos((earthday/solarsystem[4].yearratio) / 180.0 * 3.1415);
    solarsystem[3].locy = (solarsystem[0].radius/10 + solarsystem[3].distance/10) * sin((earthday/solarsystem[4].yearratio) / 180.0 * 3.1415);
    glTranslatef(solarsystem[3].locx, solarsystem[3].locy, 0);
    if(display)
      glRotatef((moonday),0,0,1);
    if(!display)
      glRotatef((earthday),0,0,1);
    DrawSphere(&solarsystem[3]);
    DrawMoon();
    glPopMatrix();
   }

   void DrawMars(void)
   {
    glPushMatrix();
    solarsystem[5].locx = (solarsystem[0].radius/10 + solarsystem[5].distance/10) * cos((earthday/solarsystem[5].yearratio) / 180.0 * 3.1415);
    solarsystem[5].locy = (solarsystem[0].radius/10 + solarsystem[5].distance/10) * sin((earthday/solarsystem[5].yearratio) / 180.0 * 3.1415);
    glTranslatef(solarsystem[5].locx, solarsystem[5].locy, 0);
    glRotatef((earthday/solarsystem[5].rotation),0,0,1);
    DrawSphere(&solarsystem[5]);
    glPopMatrix();
   }

   void DrawJupiter(void)
   {
    glPushMatrix();
    solarsystem[6].locx = (solarsystem[0].radius/10 + solarsystem[6].distance/10) * cos((earthday/solarsystem[6].yearratio) / 180.0 * 3.1415);
    solarsystem[6].locy = (solarsystem[0].radius/10 + solarsystem[6].distance/10) * sin((earthday/solarsystem[6].yearratio) / 180.0 * 3.1415);
    glTranslatef(solarsystem[6].locx, solarsystem[6].locy, 0);
    glRotatef((earthday/solarsystem[6].rotation),0,0,1);
    DrawSphere(&solarsystem[6]);
    glPopMatrix();
   }

   void DrawSaturn(void)
   {
    glPushMatrix();
    
    solarsystem[7].locx = (solarsystem[0].radius/10 + solarsystem[7].distance/10) * cos((earthday/solarsystem[7].yearratio) / 180.0 * 3.1415);
    solarsystem[7].locy = (solarsystem[0].radius/10 + solarsystem[7].distance/10) * sin((earthday/solarsystem[7].yearratio) / 180.0 * 3.1415);
    glTranslatef(solarsystem[7].locx, solarsystem[7].locy, 0);
    glRotatef((earthday/solarsystem[7].rotation),0,0,1);
    glRotatef(3,1,0,0);

      glPushMatrix();
      DrawSphere(&solarsystem[7]);
      glBindTexture(GL_TEXTURE_2D, solarsystem[7].texture);
      glPopMatrix();

     //Ring
     glPushMatrix();
     glScalef(solarsystem[7].radius/10 + 2,solarsystem[7].radius/10 + 2,0);
     DrawCylinder();
     glPopMatrix();

    glPopMatrix();
   }

   void DrawUranus(void)
   {
    glPushMatrix();
    solarsystem[8].locx = (solarsystem[0].radius/10 + solarsystem[8].distance/10) * cos((earthday/solarsystem[8].yearratio) / 180.0 * 3.1415);
    solarsystem[8].locy = (solarsystem[0].radius/10 + solarsystem[8].distance/10) * sin((earthday/solarsystem[8].yearratio) / 180.0 * 3.1415);
    glTranslatef(solarsystem[8].locx, solarsystem[8].locy, 0);
    glRotatef((earthday/solarsystem[8].rotation),0,0,1);
    DrawSphere(&solarsystem[8]);
    glPopMatrix();
   }

   void DrawNeptune(void)
   {
    glPushMatrix();
    solarsystem[9].locx = (solarsystem[0].radius/10 + solarsystem[9].distance/10) * cos((earthday/solarsystem[9].yearratio) / 180.0 * 3.1415);
    solarsystem[9].locy = (solarsystem[0].radius/10 + solarsystem[9].distance/10) * sin((earthday/solarsystem[9].yearratio) / 180.0 * 3.1415);
    glTranslatef(solarsystem[9].locx, solarsystem[9].locy, 0);
    glRotatef((earthday/solarsystem[9].rotation),0,0,1);
    DrawSphere(&solarsystem[9]);
    glPopMatrix();
   }

   void DrawPluto(void)
   {
    glPushMatrix();
    solarsystem[10].locx = (solarsystem[0].radius/10 + solarsystem[10].distance/10) * cos((earthday/solarsystem[10].yearratio) / 180.0 * 3.1415);
    solarsystem[10].locy = (solarsystem[0].radius/10 + solarsystem[10].distance/10) * sin((earthday/solarsystem[10].yearratio) / 180.0 * 3.1415);
    glTranslatef(solarsystem[10].locx, solarsystem[10].locy, 0);
    glRotatef((earthday/solarsystem[10].rotation),0,0,1);
    DrawSphere(&solarsystem[10],false,true);
    glPopMatrix();
   }

   void DrawSolarSystem(void)
   {
    if(!display)
      glFrontFace(GL_CW); //Makes the sun glow
    DrawSun();
    if(!display)
      glFrontFace(GL_CCW); //Dosen't invert the lighting for the other planets
    DrawMercury();
    DrawVenus();
    DrawEarth();
    DrawMars();
    DrawJupiter();
    DrawSaturn();
    DrawUranus();
    DrawNeptune();
    DrawPluto();
    initialized = true; //We won't have redundent textures being sent to the graphics card
   }

   virtual void RenderPiece(vtkRenderer *ren, vtkActor *act)
   {
       RemoveVTKOpenGLStateSideEffects();
       SetupLight();
       glEnable(GL_COLOR_MATERIAL);
       glEnable(GL_TEXTURE_2D);
       glMatrixMode(GL_MODELVIEW);

       DrawSolarSystem();

   }
};

vtkStandardNewMacro(vtk441Mapper);


class vtkTimerCallback : public vtkCommand
{
  /*
    This timer creates the magic of rotation, orbiting, and the position of the camera
  */
  public:
    static vtkTimerCallback *New()
    {
      vtkTimerCallback *cb = new vtkTimerCallback;
      cb->TimerCount = 0;
      cb->mapper = NULL;
      cb->renWin = NULL;
      cb->cam    = NULL;
      return cb;
    }

    void   SetMapper(vtk441Mapper *m) { mapper = m; exe = 0;};
    void   SetRenderWindow(vtkRenderWindow *rw) { renWin = rw; };
    void   SetCamera(vtkCamera *c) { cam = c; };
    virtual void Execute(vtkObject *vtkNotUsed(caller), unsigned long eventId,
                         void *vtkNotUsed(callData))
    {

      if (vtkCommand::TimerEvent == eventId)
        {
        ++this->TimerCount;
        }

      if (mapper != NULL)
            mapper->IncrementSize();

      if(exe == 300)
        moonday = 0;

      if(exe <= 150){
        moonday += 1; //Rotate the moon, needed individual to display rotation
        exe++;
        glEnable(GL_CULL_FACE);
        glFrontFace(GL_CW); // To show the textures better in first scene
      } else if(exe > 150 && exe < 300){
        moonday += 1;
        exe++;
        display = false;
      } else if (cam != NULL && exe > 300 && exe <= 900){
        cam->SetPosition(-6.48569, -71.0278, 320.251);
        cam->SetFocalPoint( 4.95772, 1.88858e-07, 6.75249);
        moonday += 1;
        earthday += 1;
        exe++;
      } else if(cam != NULL && exe > 900 && exe <= 1500) {
        cam->SetPosition(-230.864, 186.501, -107.54);
        moonday += 10;
        earthday += 10;
        exe++;
      } else if(cam != NULL && exe < 1600){ 
        cam->SetPosition(-216.655, 178.394, -1443.57);
        moonday += 50;
        earthday += 50;
        exe++;
      } else {
        moonday += 50;
        earthday += 50;
      }

      if (renWin != NULL)
        renWin->Render();

      }
 
  private:
    int TimerCount;
    vtk441Mapper *mapper;
    vtkRenderWindow *renWin;
    vtkCamera *cam;
    int   exe;
};

int main()
{
  // Dummy input so VTK pipeline mojo is happy.
  //
  vtkSmartPointer<vtkSphereSource> sphere =
    vtkSmartPointer<vtkSphereSource>::New();
  sphere->SetThetaResolution(100);
  sphere->SetPhiResolution(50);

  //Setting up background (decided to not enable for now, visibility)
  /*vtkSmartPointer<vtkImageData> imageData;
  vtkSmartPointer<vtkJPEGReader> jpegReader =
      vtkSmartPointer<vtkJPEGReader>::New();

  jpegReader->SetFileName("./textures/stars.jpg");
  jpegReader->Update();
  imageData = jpegReader->GetOutput();

  vtkSmartPointer<vtkImageActor> imageActor =
    vtkSmartPointer<vtkImageActor>::New();

  imageActor->SetInputData(imageData);

  vtkSmartPointer<vtkRenderer> backgroundRenderer = 
    vtkSmartPointer<vtkRenderer>::New();
  //*/

  vtkSmartPointer<vtk441Mapper> win3Mapper =
    vtkSmartPointer<vtk441Mapper>::New();
  win3Mapper->SetInputConnection(sphere->GetOutputPort());

  vtkSmartPointer<vtkActor> win3Actor =
    vtkSmartPointer<vtkActor>::New();
  win3Actor->SetMapper(win3Mapper);

  vtkSmartPointer<vtkRenderer> ren3 =
    vtkSmartPointer<vtkRenderer>::New();

  vtkSmartPointer<vtkRenderWindow> renWin =
    vtkSmartPointer<vtkRenderWindow>::New();

  //backgroundRenderer->SetLayer(0);
  //backgroundRenderer->InteractiveOff();
  ren3->SetLayer(0);
  renWin->SetNumberOfLayers(1);
  //renWin->AddRenderer(backgroundRenderer);
  renWin->AddRenderer(ren3);

  ren3->SetViewport(0, 0, 1, 1);

  vtkSmartPointer<vtkRenderWindowInteractor> iren =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  iren->SetRenderWindow(renWin);

  // Add the actors to the renderer, set the background and size.
  //
  ren3->AddActor(win3Actor);
  //backgroundRenderer->AddActor(imageActor);
  renWin->SetSize(2000,2000);
  renWin->Render();

  //vtkSmartPointer<vtkCallbackCommand> modifiedCallback =
  //  vtkSmartPointer<vtkCallbackCommand>::New();
  //modifiedCallback->SetCallback (CameraModifiedCallback);

  ren3->GetActiveCamera()->SetFocalPoint(4.95772, 0, 6.75249);
  ren3->GetActiveCamera()->SetPosition(73.4703, 7.26642, -5.62595);
  ren3->GetActiveCamera()->SetViewUp(1,0,0);
  ren3->GetActiveCamera()->SetDistance(70);
  ren3->GetActiveCamera()->SetClippingRange(0.1,10000);
  //ren3->GetActiveCamera()->AddObserver(vtkCommand::ModifiedEvent,modifiedCallback);

  
  ((vtkInteractorStyle *)iren->GetInteractorStyle())->SetAutoAdjustCameraClippingRange(0);
  iren->Initialize();

  vtkSmartPointer<vtkTimerCallback> cb = 
  vtkSmartPointer<vtkTimerCallback>::New();
  iren->AddObserver(vtkCommand::TimerEvent, cb);
  cb->SetMapper(win3Mapper);
  cb->SetRenderWindow(renWin);
  cb->SetCamera(ren3->GetActiveCamera());

  int timerId = iren->CreateRepeatingTimer(50);  // repeats every 10 microseconds <--> 0.01 seconds
  renWin->Render();
  iren->Start();

  return EXIT_SUCCESS;
}