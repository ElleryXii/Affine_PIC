/**
 * Implementation of the Gluvi commands to enable the viewr to run.
 *
 * @author Ante Qu, 2017
 * Based on Bridson's simple_flip2d starter code at http://www.cs.ubc.ca/~rbridson/
 */

#include <cmath>
#include <cstdarg>
#include <cstdlib>
#include <fstream>
#include "gluvi.h"

using namespace std;

namespace Gluvi{

Target3D::
Target3D(float target_[3], float dist_, float heading_, float pitch_, float fovy_, float near_clip_factor_, float far_clip_factor_)
   : dist(dist_), heading(heading_), pitch(pitch_), fovy(fovy_), 
     near_clip_factor(near_clip_factor_), far_clip_factor(far_clip_factor_), action_mode(INACTIVE)
{
   if(target_){
      target[0]=target_[0];
      target[1]=target_[1];
      target[2]=target_[2];
   }else{
      target[0]=0;
      target[1]=0;
      target[2]=0;
   }
   default_target[0]=target[0];
   default_target[1]=target[1];
   default_target[2]=target[2];
   default_dist=dist;
   default_heading=heading;
   default_pitch=pitch;
}

void Target3D::
click(int button, int state, int x, int y)
{
   if(state==GLUT_UP)
      action_mode=INACTIVE;
   else if(button==GLUT_LEFT_BUTTON)
      action_mode=ROTATE;
   else if(button==GLUT_MIDDLE_BUTTON)
      action_mode=TRUCK;
   else if(button==GLUT_RIGHT_BUTTON)
      action_mode=DOLLY;
   oldmousex=x;
   oldmousey=y;
}

void Target3D::
drag(int x, int y)
{
   switch(action_mode){
      case INACTIVE:
         return; // nothing to do
      case ROTATE:
         heading+=0.007*(oldmousex-x);
         if(heading<-M_PI) heading+=2*M_PI;
         else if(heading>M_PI) heading-=2*M_PI;
         pitch+=0.007*(oldmousey-y);
         if(pitch<-0.5*M_PI) pitch=-0.5*M_PI;
         else if(pitch>0.5*M_PI) pitch=0.5*M_PI;
         break;
      case TRUCK:
         target[0]+=(0.002*dist)*cos(heading)*(oldmousex-x);
         target[1]-=(0.002*dist)*(oldmousey-y);
         target[2]-=(0.002*dist)*sin(heading)*(oldmousex-x);
         break;
      case DOLLY:
         dist*=pow(1.01, oldmousey-y + x-oldmousex);
         break;
   }
   oldmousex=x;
   oldmousey=y;
   glutPostRedisplay();
}

void Target3D::
return_to_default(void)
{
   target[0]=default_target[0];
   target[1]=default_target[1];
   target[2]=default_target[2];
   dist=default_dist;
   heading=default_heading;
   pitch=default_pitch;
}

void Target3D::
transform_mouse(int x, int y, float ray_origin[3], float ray_direction[3])
{
   float ch=cos(heading), sh=sin(heading);
   float cp=cos(pitch), sp=sin(pitch);

   ray_origin[0]=target[0]+dist*sh*cp;
   ray_origin[1]=target[1]-dist*sp;
   ray_origin[2]=target[2]+dist*ch*cp;

   float scale=0.5*tan(fovy)/winheight;
   float camx=(x-0.5*winwidth)*scale, camy=(0.5*winheight-y)*scale, camz=-1.0; // in camera coordinates, this is ray_direction (but not normalized)
   // now need to rotate into world space from camera space
   float px=camx, py=camy*cp-camz*sp, pz=camy*sp+camz*cp;
   ray_direction[0]=px*ch+pz*sh;
   ray_direction[1]=py;
   ray_direction[2]=-px*sh+pz*ch;

   float mag=sqrt(ray_direction[0]*ray_direction[0]
                  + ray_direction[1]*ray_direction[1]
                  + ray_direction[2]*ray_direction[2]);
   ray_direction[0]/=mag;
   ray_direction[1]/=mag;
   ray_direction[2]/=mag;

   ray_origin[0]+=near_clip_factor*dist*ray_direction[0];
   ray_origin[1]+=near_clip_factor*dist*ray_direction[1];
   ray_origin[2]+=near_clip_factor*dist*ray_direction[2];
}

void Target3D::
gl_transform(void)
{
   glViewport(0, 0, (GLsizei)winwidth, (GLsizei)winheight);

   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   gluPerspective(fovy, winwidth/(float)winheight, near_clip_factor*dist, far_clip_factor*dist);

   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   GLfloat pos[3];
   pos[0]=target[0]-dist*sin(heading)*cos(pitch);
   pos[1]=target[1]-dist*sin(pitch);
   pos[2]=target[2]-dist*cos(heading)*cos(pitch);
   glTranslatef(0, 0, -dist); // translate target dist away in the z direction
   glRotatef(-180/M_PI*pitch, 1, 0, 0); // rotate pitch in the yz plane
   glRotatef(-180/M_PI*heading, 0, 1, 0); // rotate heading in the xz plane
   glTranslatef(-target[0], -target[1], -target[2]); // translate target to origin
}

void Target3D::
export_rib(ostream &output)
{
   output<<"Clipping "<<near_clip_factor*dist<<" "<<far_clip_factor*dist<<endl; // could be more generous here!
   output<<"Projection \"perspective\" \"fov\" "<<fovy<<endl;
   output<<"ReverseOrientation"<<endl;  // RenderMan has a different handedness from OpenGL's default
   output<<"Scale 1 1 -1"<<endl;        // so we need to correct for that here
   output<<"Translate 0 0 "<<-dist<<endl;
   output<<"Rotate "<<-180/M_PI*pitch<<" 1 0 0"<<endl;
   output<<"Rotate "<<-180/M_PI*heading<<" 0 1 0"<<endl;
   output<<"Translate "<<-target[0]<<" "<<-target[1]<<" "<<-target[2]<<endl;
}

//=================================================================================

PanZoom2D::
PanZoom2D(float bottom_, float left_, float height_)
   : bottom(bottom_), left(left_), height(height_), action_mode(INACTIVE)
{
   default_bottom=bottom;
   default_left=left;
   default_height=height;
}

void PanZoom2D::
click(int button, int state, int x, int y)
{
   if(state==GLUT_UP){
      float r=height/winheight;
      switch(action_mode){
         case PAN:
            if(!moved_since_mouse_down){
               // make mouse click the centre of the window
               left+=r*(x-0.5*winwidth);
               bottom+=r*(0.5*winheight-y);
               glutPostRedisplay();
            }
            break;
         case ZOOM_IN:
            if(moved_since_mouse_down){
               // zoom in to selection
               float desired_width=fabs((x-clickx)*height/winheight);
               float desired_height=fabs((y-clicky)*height/winheight);
               if(desired_height==0) desired_height=height/winheight;
               if(desired_width*winheight > desired_height*winwidth)
                  desired_height=winheight*desired_width/winwidth;
               else
                  desired_width=winwidth*desired_height/winheight;
               left+=0.5*(x+clickx)*height/winheight-0.5*desired_width;
               bottom+=(winheight-0.5*(y+clicky))*height/winheight-0.5*desired_height;
               height=desired_height;
            }else{
               // zoom in by some constant factor on the mouse click
               float factor=0.70710678118654752440084;
               left+=(1-factor)*height*(x/(float)winheight);
               bottom+=(1-factor)*height*(1-y/(float)winheight);
               height*=factor;
            }
            glutPostRedisplay();
            break;
         case ZOOM_OUT:
            // zoom out by some constant factor
            {
               float factor=1.41421356237309504880168;
               left-=0.5*(factor-1)*height;
               bottom-=0.5*(factor-1)*winwidth*height/winheight;
               height*=factor;
            }
            glutPostRedisplay();
            break;
         default:
            ;// nothing to do
      }
      action_mode=INACTIVE;

   }else if(button==GLUT_LEFT_BUTTON)
      action_mode=PAN;
   else if(button==GLUT_MIDDLE_BUTTON){
      clickx=x;
      clicky=y;
      action_mode=ZOOM_IN;
   }else if(button==GLUT_RIGHT_BUTTON)
      action_mode=ZOOM_OUT;
   moved_since_mouse_down=false;
   oldmousex=x;
   oldmousey=y;
}

void PanZoom2D::
drag(int x, int y)
{
   if(x!=oldmousex || y!=oldmousey){
      moved_since_mouse_down=true;
      if(action_mode==PAN){
         float r=height/winheight;
         left-=r*(x-oldmousex);
         bottom+=r*(y-oldmousey);
         glutPostRedisplay();
      }
      oldmousex=x;
      oldmousey=y;
   }
}

void PanZoom2D::
return_to_default(void)
{
   bottom=default_bottom;
   left=default_left;
   height=default_height;
}

void PanZoom2D::
transform_mouse(int x, int y, float coords[2])
{
   float r=height/winheight;
   coords[0]=x*r+left;
   coords[1]=(winheight-y)*r+bottom;
}

void PanZoom2D::
gl_transform(void)
{
   glViewport(0, 0, (GLsizei)winwidth, (GLsizei)winheight);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glOrtho(left, left+(height*winwidth)/winheight, bottom, bottom+height, 0, 1);
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
}

void PanZoom2D::
export_rib(ostream &output)
{
   // no projection matrix
   output<<"Clipping 1 2000"<<endl; // somewhat arbitrary - hopefully this is plenty of space
   output<<"ReverseOrientation"<<endl;  // RenderMan has a different handedness from OpenGL's default
   output<<"Scale 1 1 -1"<<endl;        // so we need to correct for that here
   // scale so that smaller dimension gets scaled to size 2
   float scalefactor;
   if(winwidth>winheight) scalefactor=2.0/height;
   else scalefactor=2.0/(winwidth*height/winheight);
   output<<"Scale "<<scalefactor<<" "<<scalefactor<<" 1"<<endl;
   // translate so centre of view gets mapped to (0,0,1000)
   output<<"Translate "<<-(left+0.5*winwidth*height/winheight)<<" "<<-(bottom+0.5*height)<< " 1000"<<endl;
}

//=================================================================================

StaticText::
StaticText(const char *text_)
   : text(text_)
{}

void StaticText::
display(int x, int y)
{
   dispx=x;
   dispy=y;
   width=glutBitmapLength(GLUT_BITMAP_HELVETICA_12, (const unsigned char*)text)+1;
   height=15;
   glColor3f(0.3, 0.3, 0.3);
   glRasterPos2i(x, y-height+2);
   for(int i=0; text[i]!=0; ++i)
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, text[i]);
   glColor3f(1, 1, 1);
   glRasterPos2i(x+1, y-height+3);
   for(int i=0; text[i]!=0; ++i)
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, text[i]);
}

//=================================================================================

Button::
Button(const char *text_, int minwidth_)
   : status(UNINVOLVED), text(text_), minwidth(minwidth_)
{}

void Button::
display(int x, int y)
{
   dispx=x;
   dispy=y;
   int textwidth=glutBitmapLength(GLUT_BITMAP_HELVETICA_12, (const unsigned char*)text);
   if(textwidth<minwidth) width=minwidth+24;
   else width=textwidth+24;
   height=17;
   if(status==UNINVOLVED){
      glColor3f(0.7, 0.7, 0.7);
      glBegin(GL_QUADS);
      glVertex2i(x+1, y-1);
      glVertex2i(x+width, y-1);
      glVertex2i(x+width, y-height+1);
      glVertex2i(x+1, y-height+1);
      glEnd();
      glColor3f(0.3, 0.3, 0.3);
      glLineWidth(1);
      glBegin(GL_LINE_STRIP);
      glVertex2i(x, y-2);
      glVertex2i(x, y-height);
      glVertex2i(x+width-1, y-height);
      glEnd();
      glColor3f(0.3, 0.3, 0.3);
   }else{
      if(status==SELECTED) glColor3f(0.8, 0.8, 0.8);
      else                 glColor3f(1, 1, 1);
      glBegin(GL_QUADS);
      glVertex2i(x, y-1);
      glVertex2i(x+width, y-1);
      glVertex2i(x+width, y-height);
      glVertex2i(x, y-height);
      glEnd();
      glColor3f(0, 0, 0);
   }
   glRasterPos2i(x+(width-textwidth)/2, y-height+5);
   for(int i=0; text[i]!=0; ++i)
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, text[i]);
}

bool Button::
click(int state, int x, int y)
{
   if(state==GLUT_DOWN && x>dispx && x<=dispx+width && y<dispy-2 && y>=dispy-height){
      status=HIGHLIGHTED;
      glutPostRedisplay();
      return true;
   }else if(state==GLUT_UP && status!=UNINVOLVED){
      status=UNINVOLVED;
      glutPostRedisplay();
      if(x>=dispx && x<dispx+width && y<dispy-2 && y>=dispy-height)
         action();
      return true;
   }else
      return false;
}

void Button::
drag(int x, int y)
{
   // needs to control highlighting (SELECTED vs. HIGHLIGHTED)
   if(status==SELECTED && x>=dispx && x<dispx+width && y<dispy-2 && y>=dispy-height){
      status=HIGHLIGHTED;
      glutPostRedisplay();
   }else if(status==HIGHLIGHTED && !(x>=dispx && x<dispx+width && y<dispy-2 && y>=dispy-height)){
      status=SELECTED;
      glutPostRedisplay();
   }
}

//=================================================================================

Slider::
Slider(const char *text_, int length_, int position_, int justify_)
   : status(UNINVOLVED), text(text_), length(length_), justify(justify_), position(position_)
{}

void Slider::
display(int x, int y)
{
   dispx=x;
   dispy=y;
   width=glutBitmapLength(GLUT_BITMAP_HELVETICA_12, (const unsigned char*)text);
   if(width<justify) width=justify;
   width+=11+6+length+1;
   height=15;
   glColor3f(0.3, 0.3, 0.3);
   glRasterPos2i(x, y-height+2);
   for(int i=0; text[i]!=0; ++i)
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, text[i]);
   glColor3f(1, 1, 1);
   glRasterPos2i(x+1, y-height+3);
   for(int i=0; text[i]!=0; ++i)
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, text[i]);
   scrollxmin=x+width-length-12;
   scrollxmax=x+width;
   scrollymin=y-height+1;
   scrollymax=y-2;
   glColor3f(0.3, 0.3, 0.3);
   glLineWidth(1);
   glBegin(GL_LINE_STRIP);
   glVertex2i(scrollxmin, scrollymax-1);
   glVertex2i(scrollxmin, scrollymin);
   glVertex2i(scrollxmax-1, scrollymin);
   glVertex2i(scrollxmax-1, scrollymax-1);
   glEnd();
   glColor3f(0.7, 0.7, 0.7);
   glBegin(GL_LINE_STRIP);
   glVertex2i(scrollxmin+1, scrollymax);
   glVertex2i(scrollxmin+1, scrollymin+1);
   glVertex2i(scrollxmax, scrollymin+1);
   glVertex2i(scrollxmax, scrollymax);
   glEnd();
   if(status==UNINVOLVED){
      glColor3f(0.3, 0.3, 0.3);
      glBegin(GL_LINE_STRIP);
      glVertex2i(scrollxmin+position+2, scrollymax-2);
      glVertex2i(scrollxmin+position+2, scrollymin+2);
      glVertex2i(scrollxmin+position+10, scrollymin+2);
      glEnd();
      glColor3f(0.7, 0.7, 0.7);
      glBegin(GL_QUADS);
      glVertex2i(scrollxmin+position+3, scrollymin+3);
      glVertex2i(scrollxmin+position+11, scrollymin+3);
      glVertex2i(scrollxmin+position+11, scrollymax);
      glVertex2i(scrollxmin+position+3, scrollymax);
      glEnd();
   }else{ // SELECTED
      glColor3f(1, 1, 1);
      glBegin(GL_QUADS);
      glVertex2i(scrollxmin+position+2, scrollymin+2);
      glVertex2i(scrollxmin+position+11, scrollymin+2);
      glVertex2i(scrollxmin+position+11, scrollymax);
      glVertex2i(scrollxmin+position+2, scrollymax);
      glEnd();
   }
}

bool Slider::
click(int state, int x, int y)
{
   if(state==GLUT_DOWN && x>scrollxmin+position+2 && x<=scrollxmin+position+11 && y<scrollymax-1 && y>=scrollymin+2){
      status=SELECTED;
      clickx=x;
      glutPostRedisplay();
      return true;
   }else if(status!=UNINVOLVED && state==GLUT_UP){
      status=UNINVOLVED;
      glutPostRedisplay();
      return true;
   }else
      return false;
}

void Slider::
drag(int x, int y)
{
   if(status==SELECTED){
      glutPostRedisplay();
      int newposition=position+(x-clickx);
      clickx=x;
      if(newposition<0){
         clickx+=(0-newposition);
         newposition=0;
      }else if(newposition>length){
         clickx+=(length-newposition);
         newposition=length;
      }
      if(newposition!=position){
         position=newposition;
         action();
         glutPostRedisplay();
      }
   }
}

//=================================================================================

WidgetList::
WidgetList(int indent_, bool hidden_)
   : indent(indent_), hidden(hidden_), downclicked_member(-1)
{
}

void WidgetList::
display(int x, int y)
{
   dispx=x;
   dispy=y;
   if(hidden){
      width=height=0;
   }else{
      height=0;
      for(unsigned int i=0; i<list.size(); ++i){
         list[i]->display(x+indent, y-height);
         height+=list[i]->height;
         width=(width<indent+list[i]->width) ? indent+list[i]->width : width;
      }
   }
}

bool WidgetList::
click(int state, int x, int y)
{
   //if(hidden || x<dispx || x>=dispx+width || y>=dispy || y<dispy-height) return false; // early exit
   if(state==GLUT_DOWN){ // search for correct widget
      for(unsigned int i=0; i<list.size(); ++i){
         if(list[i]->click(state, x, y)){
            downclicked_member=i;
            return true;
         }
      }
   }else if(state==GLUT_UP && downclicked_member>=0){
      list[downclicked_member]->click(state, x, y);
      downclicked_member=-1;
   }
   return false;
}

void WidgetList::
drag(int x, int y)
{
   if(downclicked_member>=0)
      list[downclicked_member]->drag(x, y);
}

//=================================================================================

static void gluviReshape(int w, int h)
{
   winwidth=w;
   winheight=h;
   glutPostRedisplay(); // triggers the camera to adjust itself to the new dimensions
}

//=================================================================================

static void gluviDisplay()
{
   glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

   // draw the scene
   if(camera) camera->gl_transform();
   if(userDisplayFunc) userDisplayFunc();

   // now draw widgets on top
   glPushAttrib(GL_CURRENT_BIT|GL_ENABLE_BIT|GL_LINE_BIT);
   glDisable(GL_DEPTH_TEST);
   glDisable(GL_LIGHTING);
   glLineWidth(1);
   // and probably more needs setting before widgets

   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   gluOrtho2D(0, winwidth, 0, winheight);

   root.display(0, winheight);

   glPopAttrib();

   glutSwapBuffers();
}

//=================================================================================

static enum {NOBODY, CAMERA, WIDGETS, USER} mouse_owner=NOBODY;

static void gluviMouse(int button, int state, int x, int y)
{
   if(state==GLUT_DOWN){
      int mods=glutGetModifiers();
      if(camera && mods==GLUT_ACTIVE_SHIFT){
         camera->click(button, state, x, y);
         mouse_owner=CAMERA;
      }else if(button==GLUT_LEFT_BUTTON && root.click(state, x, winheight-y)){
         mouse_owner=WIDGETS;
      }else if(userMouseFunc){
         userMouseFunc(button, state, x, y);
         mouse_owner=USER;
      }
   }else{ // mouse up - send event to whoever got the mouse down
      switch(mouse_owner){
         case CAMERA:
            camera->click(button, state, x, y);
            break;
         case WIDGETS:
            root.click(state, x, winheight-y);
            break;
         case USER:
            if(userMouseFunc) userMouseFunc(button, state, x, y);
            break;
         default:
           ;// nothing to do
      }
      mouse_owner=NOBODY;
   }
}

//=================================================================================

static void gluviDrag(int x, int y)
{
   switch(mouse_owner){
      case CAMERA:
         camera->drag(x, y);
         break;
      case WIDGETS:
         root.drag(x, winheight-y);
         break;
      case USER:
         if(userDragFunc) userDragFunc(x, y);
         break;
      default:
         ;// nothing to do
   }
}

//=================================================================================

void ppm_screenshot(const char *filename_format, ...)
{
   va_list ap;
   va_start(ap, filename_format);
#ifdef _MSC_VER
#define FILENAMELENGTH 256
   char filename[FILENAMELENGTH];
   _vsnprintf(filename, FILENAMELENGTH, filename_format, ap);
   ofstream out(filename, ofstream::binary);
#else
   char *filename;
   vasprintf(&filename, filename_format, ap);
   ofstream out(filename, ofstream::binary);
   free(filename);
#endif
   if(!out) return;
   GLubyte *image_buffer=new GLubyte[3*winwidth*winheight];
   glReadBuffer(GL_FRONT);
   glReadPixels(0, 0, winwidth, winheight, GL_RGB, GL_UNSIGNED_BYTE, image_buffer);
   out<<"P6\n"<<winwidth<<' '<<winheight<<" 255\n";
   for(int i=1; i<=winheight; ++i)
      out.write((const char*)image_buffer+3*winwidth*(winheight-i), 3*winwidth);
   delete[] image_buffer;
}

void set_generic_lights(void)
{
   glEnable(GL_LIGHTING);
   {
      GLfloat ambient[4] = {.3, .3, .3, 1};
      glLightModelfv (GL_LIGHT_MODEL_AMBIENT,ambient);
   }
   {
      GLfloat color[4] = {.8, .8, .8, 1};
      glLightfv (GL_LIGHT0, GL_DIFFUSE, color);
      glLightfv (GL_LIGHT0, GL_SPECULAR, color);
      glEnable (GL_LIGHT0);
   }
   {
      GLfloat color[4] = {.4, .4, .4, 1};
      glLightfv (GL_LIGHT1, GL_DIFFUSE, color);
      glLightfv (GL_LIGHT1, GL_SPECULAR, color);
      glEnable (GL_LIGHT1);
   }
   {
      GLfloat color[4] = {.2, .2, .2, 1};
      glLightfv (GL_LIGHT2, GL_DIFFUSE, color);
      glLightfv (GL_LIGHT2, GL_SPECULAR, color);
      glEnable (GL_LIGHT2);
   }
}

void set_generic_material(float r, float g, float b)
{
   GLfloat ambient[4], diffuse[4], specular[4];
   ambient[0]=r*r+0.1; ambient[1]=g*g+0.1; ambient[2]=b*b+0.1; ambient[3]=0;
   diffuse[0]=r; diffuse[1]=g; diffuse[2]=b; diffuse[3]=0;
   specular[0]=sqrt(r)+0.1; specular[1]=sqrt(g)+0.1; specular[2]=sqrt(b)+0.1; specular[3]=0;
   glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
   glMaterialfv (GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
   glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, specular);
   glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, 32);
}

//=================================================================================

void init(const char *windowtitle, int *argc, char **argv)
{
   glutInit(argc, argv);
   glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
   glutInitWindowSize(winwidth, winheight);
   glutCreateWindow(windowtitle);
   glutReshapeFunc(gluviReshape);
   glutDisplayFunc(gluviDisplay);
   glutMouseFunc(gluviMouse);
   glutMotionFunc(gluviDrag);
   glEnable(GL_DEPTH_TEST);
   glClearColor(0, 0, 0, 0);
   glClearDepth(1);
   glPixelStorei(GL_PACK_ALIGNMENT, 1);
   glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
}

//=================================================================================

void (*userDisplayFunc)(void)=0; 
void (*userMouseFunc)(int button, int state, int x, int y)=0;
void (*userDragFunc)(int x, int y)=0;
Camera *camera=0;
WidgetList root(0);
int winwidth=640, winheight=480;

//=================================================================================

void run(void)
{
   glutMainLoop();
}

};

