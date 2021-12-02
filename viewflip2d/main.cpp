/**
 * Main class for the viewer for simulation results
 *
 * @author Ante Qu, 2017
 * Based on Bridson's simple_flip2d starter code at http://www.cs.ubc.ca/~rbridson/
 */

#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include "gluvi.h"
#include "vec2.h"

using namespace std;

char frame_number[100]="frame 0";
const char *file_format;
unsigned int frame=0;
vector<Vec2f> x;

bool read_frame(int newframe)
{
   if(newframe<0) return false;

   char filename[100];
   sprintf(filename, file_format, newframe);
   ifstream in(filename);
   if(!in.good())
      return false;
   int n;
   in>>n;
   x.resize(n);
   for(int i=0; i<n; ++i)
      in>>x[i];
   frame=newframe;
   sprintf(frame_number, "frame %d", frame);
   return true;
}

void set_view(float &bottom, float &left, float &height)
{
   bottom=x[0][1];
   left=x[0][0];
   float top=bottom, right=left;

   for(unsigned int i=1; i<x.size(); ++i){
      if(x[i][0]<left) left=x[i][0];
      else if(x[i][0]>right) right=x[i][0];
      if(x[i][1]<bottom) bottom=x[i][1];
      else if(x[i][1]>top) top=x[i][1];
   }
   if (right < 1) right = 1;
   if (top < 1) top = 1;
   if(right-left > top-bottom)
      height=1.5*(right-left);
   else
      height=1.5*(top-bottom);
   

   left=(left+right)/2-height/2;
   bottom=(top+bottom)/2-height/2;
}

void display(void)
{
   glDisable(GL_LIGHTING);
   glColor3f(1, 1, 1);
   glBegin(GL_POINTS);
   for(unsigned int i=0; i<x.size(); ++i)
      glVertex2fv(x[i].v);
   glEnd();
}

struct ScreenShotButton : public Gluvi::Button{
   const char *filename_format;
   ScreenShotButton(const char *label, const char *filename_format_) : Gluvi::Button(label), filename_format(filename_format_) {}
   void action()
   { Gluvi::ppm_screenshot(filename_format, frame); }
};

void special_key_handler(int key, int x, int y)
{
   switch(key){
      case GLUT_KEY_LEFT:
         if(read_frame(frame-1))
            glutPostRedisplay();
         break;
      case GLUT_KEY_RIGHT:
         if(read_frame(frame+1))
            glutPostRedisplay();
         break;
      default:
         ;
   }
}

int main(int argc, char **argv)
{
   Gluvi::init("viewflip2d", &argc, argv);

   if(argc!=2){
      cerr<<"Expecting one argument: format for particle filenames"<<endl;
      return 1;
   }

   file_format=argv[1];   
   if(!read_frame(0))
      return 1;

   glutSpecialFunc(special_key_handler);

   float bottom, left, height;
   set_view(bottom, left, height);
   Gluvi::PanZoom2D cam(bottom, left, height);
   Gluvi::camera=&cam;

   Gluvi::userDisplayFunc=display;

   Gluvi::StaticText frametext(frame_number);
   Gluvi::root.list.push_back(&frametext);

   char ppmfileformat[strlen(file_format)+5];
   sprintf(ppmfileformat, "%s.ppm", file_format);
   ScreenShotButton screenshot("Screenshot", ppmfileformat);
   Gluvi::root.list.push_back(&screenshot);

   Gluvi::run();
   return 0;
}

