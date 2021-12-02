/**
 * File for setting up the main simulator with viewer. Hit spacebar to advance forward one frame.
 *
 * @author Ante Qu, 
 * Based on Bridson's simple_flip2d starter code at http://www.cs.ubc.ca/~rbridson/
 */

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include "particles.h"
#include "util.h"
#include "viewflip2d/gluvi.h"
#include "shared_main.h"

using namespace std;

char frame_number[100]="frame 0";
Grid* pGrid;
Particles* pParticles;
std::string outputpath;
double frametime = 1./30;
int stepCount;



void set_view(float &bottom, float &left, float &height)
{
   bottom=0;
   left=0;
   float top=bottom, right=left;
   right = 1;
   top = 1;
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
   if( pParticles )
      for(unsigned int i=0; i<pParticles->x.size(); ++i)
         glVertex2fv(pParticles->x[i].v);
   glEnd();
}

struct ScreenShotButton : public Gluvi::Button{
   
   ScreenShotButton(const char *label) : Gluvi::Button(label) {}
   void action()
   { ; }
};

void key_handler(unsigned char key, int x, int y)
{
   if( !pGrid || !pParticles)
      return;
   switch(key){
      case ' ':
         advance_one_frame(*pGrid, *pParticles, 1./30);
         ++stepCount;
         printf("===================================================> step %d...\n", stepCount);
         pParticles->write_to_file("%s/frameparticles%04d", outputpath.c_str(), stepCount);
         sprintf(frame_number, "frame %d", stepCount);
         glutPostRedisplay();
         break;
      default:
         ;
   }
}


int main(int argc, char **argv)
{
   float gravity = 9.8;
   if( USE_SPHERICAL_GRAV )
      gravity *= GRAV_FACTOR;
   pGrid = new Grid(gravity, 50, 50, 1);
   SimulationType sType = SIMULATION_TYPE;
   
   outputpath=".";

   if(argc>1) outputpath=argv[1];
   else printf("using default output path...\n");
   printf("Output sent to %s\n", outputpath.c_str() );

   if(argc>2){
      std::string  simType = argv[2];
      std::transform(simType.begin(), simType.end(), simType.begin(), ::tolower);
      if (!simType.compare("apic"))
         sType = APIC;
      else if(!simType.compare("flip"))
         sType = FLIP;
      else if(!simType.compare("pic"))
         sType = PIC;
   }
   pParticles = new Particles(*pGrid, sType);

   Gluvi::init("fluid simulation viewer woohoo", &argc, argv);
   init_water_drop(*pGrid, *pParticles, 2, 2);
   pParticles->write_to_file("%s/frameparticles%04d", outputpath.c_str(), 0);
   stepCount = 0;

   glutKeyboardFunc(key_handler);

   float bottom, left, height;
   set_view(bottom, left, height);
   Gluvi::PanZoom2D cam(bottom, left, height);
   Gluvi::camera=&cam;
   
   Gluvi::userDisplayFunc=display;

   Gluvi::StaticText frametext(frame_number);
   Gluvi::root.list.push_back(&frametext);
/*
   char ppmfileformat[strlen(file_format)+5];
   sprintf(ppmfileformat, "%s.ppm", file_format);
   ScreenShotButton screenshot("Screenshot", ppmfileformat);
   Gluvi::root.list.push_back(&screenshot);
*/
   Gluvi::run();
   return 0;
}


