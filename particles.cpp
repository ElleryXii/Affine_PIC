/**
 * Implementation of the Particles functions. Most of your edits should be in here.
 *
 * @author Ante Qu, 2017
 * Based on Bridson's simple_flip2d starter code at http://www.cs.ubc.ca/~rbridson/
 */

#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include "particles.h"
#include "util.h"

using namespace std;

void Particles::
add_particle(const Vec2f &px, const Vec2f &pu)
{
   x.push_back(px);
   u.push_back(pu);
   /* TODO: initialize the variables you created in particles.h */
   //cx.push_back(Vec2f(0.f,0.f));
   //cy.push_back(Vec2f(0.f,0.f));
   ++np;
}

template<class T>
void Particles::
accumulate(T &accum, float q, int i, int j, float fx, float fy)
{
   float weight;

   weight=(1-fx)*(1-fy);
   accum(i,j)+=weight*q;
   sum(i,j)+=weight;

   weight=fx*(1-fy);
   accum(i+1,j)+=weight*q;
   sum(i+1,j)+=weight;

   weight=(1-fx)*fy;
   accum(i,j+1)+=weight*q;
   sum(i,j+1)+=weight;

   weight=fx*fy;
   accum(i+1,j+1)+=weight*q;
   sum(i+1,j+1)+=weight;
}

/* call this function to incorporate c[] when transfering particles to grid */
/* This function should take the c_pa^n values from c, and update them, with proper weighting, */
/*  into the correct grid velocity values in accum */
template<class T>
void Particles::
affineFix(T &accum, Vec2f c, int i, int j, float fx, float fy)
{

   /* TODO: fill this in */

   accum(i,j)+= 0; 
 
   accum(i+1,j)+= 0;

   accum(i,j+1)+= 0;

   accum(i+1,j+1)+= 0;
}

void Particles::
transfer_to_grid(void)
{
   int p, i, ui, j, vj;
   float fx, ufx, fy, vfy;

   grid.u.zero();
   sum.zero();
   for(p=0; p<np; ++p){
      grid.bary_x(x[p][0], ui, ufx);
      grid.bary_y_centre(x[p][1], j, fy);
      accumulate(grid.u, u[p][0], ui, j, ufx, fy);
      /* TODO: call affineFix to incorporate c_px^n into the grid.u update */
      // affineFix([FILL THIS IN]);
   }
   for(j=0; j<grid.u.ny; ++j) for(i=0; i<grid.u.nx; ++i){
      if(sum(i,j)!=0) grid.u(i,j)/=sum(i,j);
   }

   grid.v.zero();
   sum.zero();
   for(p=0; p<np; ++p){
      grid.bary_x_centre(x[p][0], i, fx);
      grid.bary_y(x[p][1], vj, vfy);
      accumulate(grid.v, u[p][1] , i, vj, fx, vfy);
      /* TODO: call affineFix to incorporate c_py^n into the grid.v update */
      // affineFix([FILL THIS IN]);
   }
   for(j=0; j<grid.v.ny; ++j) for(i=0; i<grid.v.nx; ++i){
      if(sum(i,j)!=0) grid.v(i,j)/=sum(i,j);
   }

   // identify where particles are in grid
   grid.marker.zero();
   for(p=0; p<np; ++p){
      grid.bary_x(x[p][0], i, fx);
      grid.bary_y(x[p][1], j, fy);
      grid.marker(i,j)=FLUIDCELL;
   }
}

/* this function computes c from the gradient of w and the velocity field from the grid. */
Vec2f Particles::
computeC(Array2f &ufield, int i, int j, float fx, float fy)
{
   /* TODO: fill this in */
   return Vec2f(0.f,0.f);
}

void Particles::
update_from_grid(void)
{
   int p;
   int i, ui, j, vj;
   float fx, ufx, fy, vfy;
   for(p=0; p<np; ++p){
      grid.bary_x(x[p][0], ui, ufx);
      grid.bary_x_centre(x[p][0], i, fx);
      grid.bary_y(x[p][1], vj, vfy);
      grid.bary_y_centre(x[p][1], j, fy);
      if( simType == FLIP )
      {
         u[p]+=Vec2f(grid.du.bilerp(ui, j, ufx, fy), grid.dv.bilerp(i, vj, fx, vfy)); // FLIP
      }
      else
      {
         u[p]=Vec2f(grid.u.bilerp(ui, j, ufx, fy), grid.v.bilerp(i, vj, fx, vfy)); // PIC and APIC

         if( simType == APIC )
         {
            /* TODO: call computeC with the right indices to compute c_px^n and c_py^n */
            //cx[p] = [FILL THIS IN]; // APIC
            //cy[p] = [FILL THIS IN]; // APIC
         }
      }
      
   }
}

void Particles::
move_particles_in_grid(float dt)
{
   Vec2f midx, gu;
   float xmin=1.001*grid.h, xmax=grid.lx-1.001*grid.h;
   float ymin=1.001*grid.h, ymax=grid.ly-1.001*grid.h;
   for(int p=0; p<np; ++p){
      // first stage of Runge-Kutta 2 (do a half Euler step)
      grid.bilerp_uv(x[p][0], x[p][1], gu[0], gu[1]);
      midx=x[p]+0.5*dt*gu;
      clamp(midx[0], xmin, xmax);
      clamp(midx[1], ymin, ymax);
      // second stage of Runge-Kutta 2
      grid.bilerp_uv(midx[0], midx[1], gu[0], gu[1]);
      x[p]+=dt*gu;
      clamp(x[p][0], xmin, xmax);
      clamp(x[p][1], ymin, ymax);
   }
}

void Particles::
write_to_file(const char *filename_format, ...)
{
   va_list ap;
   va_start(ap, filename_format);
   char *filename;
   vasprintf(&filename, filename_format, ap);
   FILE *fp=fopen(filename, "wt");
   free(filename);
   va_end(ap);

   fprintf(fp, "%d\n", np);
   for(int p=0; p<np; ++p)
      fprintf(fp, "%.5g %.5g\n", x[p][0], x[p][1]);
   fclose(fp);
}

