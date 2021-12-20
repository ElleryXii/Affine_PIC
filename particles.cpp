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
#include <iostream>
using namespace std;

void Particles::
add_particle(const Vec3f &px, const Vec3f &pu)
{
   x.push_back(px);
   u.push_back(pu);
   cx.push_back(Vec3f(0.f,0.f,0.f));
   cy.push_back(Vec3f(0.f,0.f,0.f));
   cz.push_back(Vec3f(0.f,0.f,0.f));
   ++np;
}

template<class T>
void Particles::
accumulate(T &accum, float q, int i, int j, int k, float fx, float fy, float fz)
{
   float weight;

   weight=(1-fx)*(1-fy)*fz;
   accum(i,j,k+1)+=weight*q;
   sum(i,j,k+1)+=weight;

   weight=fx*(1-fy)*fz;
   accum(i+1,j,k+1)+=weight*q;
   sum(i+1,j,k+1)+=weight;

   weight=(1-fx)*fy*fz;
   accum(i,j+1,k+1)+=weight*q;
   sum(i,j+1,k+1)+=weight;

   weight=fx*fy*fz;
   accum(i+1,j+1,k+1)+=weight*q;
   sum(i+1,j+1,k+1)+=weight;

   weight=(1-fx)*(1-fy)*(1-fz);
   accum(i,j,k)+=weight*q;
   sum(i,j,k)+=weight;

   weight=fx*(1-fy)*(1-fz);
   accum(i+1,j,k)+=weight*q;
   sum(i+1,j,k)+=weight;

   weight=(1-fx)*fy*(1-fz);
   accum(i,j+1,k)+=weight*q;
   sum(i,j+1,k)+=weight;

   weight=fx*fy*(1-fz);
   accum(i+1,j+1,k)+=weight*q;
   sum(i+1,j+1,k)+=weight;
}

/* call this function to incorporate c[] when transfering particles to grid */
/* This function should take the c_pa^n values from c, and update them, with proper weighting, */
/*  into the correct grid velocity values in accum */
template<class T>
void Particles::
affineFix(T &accum, Vec3f c, int i, int j, int k, float fx, float fy, float fz)
{
    float weight;
    weight = (1-fx)*(1-fy)*(1-fz);
    accum(i,j,k) +=  weight * (c[0] * fx + c[1] * fy + c[2] * fz);

    weight = fx*(1-fy)*(1-fz);
    accum(i+1,j,k) += weight * (c[0] * (1-fx) + c[1] * fy + c[2] * fz);

    weight = (1-fx)*fy*(1-fz);
    accum(i,j+1,k) += weight * (c[0] * fx + c[1] * (1-fy) + c[2] * fz);

    weight = fx*fy*(1-fz);
    accum(i+1, j+1,k) += weight * (c[0]*(1-fx) + c[1] * (1-fy) + c[2] * fz);

    weight = (1-fx)*(1-fy)*fz;
    accum(i,j,k+1) +=  weight * (c[0] * fx + c[1] * fy + c[2] * (1-fz));

    weight = fx*(1-fy)*fz;
    accum(i+1,j,k+1) += weight * (c[0] * (1-fx) + c[1] * fy + c[2] * (1-fz));

    weight = (1-fx)*fy*fz;
    accum(i,j+1,k+1) += weight * (c[0] * fx + c[1] * (1-fy) + c[2] * (1-fz));

    weight = fx*fy*fz;
    accum(i+1, j+1, k+1) += weight * (c[0]*(1-fx) + c[1] * (1-fy) + c[2] * (1-fz));
}

void Particles::
transfer_to_grid(void)
{
   int p, i, ui, j, vj, k, wk;
   float fx, ufx, fy, vfy, fz, wfz;

   grid.u.zero();
   sum.zero();
   for(p=0; p<np; ++p){
      grid.bary_x(x[p][0], ui, ufx);
      grid.bary_y_centre(x[p][1], j, fy);
      grid.bary_z_centre(x[p][2], k, fz);
      accumulate(grid.u, u[p][0], ui, j, k, ufx, fy, fz);
      affineFix(grid.u, cx[p], ui, j, k, ufx, fy, fz);
   }

   for(k=0; k<grid.u.nz; ++k){
      for(j=0; j<grid.u.ny; ++j){ 
         for(i=0; i<grid.u.nx; ++i){
            if(sum(i,j,k)!=0) grid.u(i,j,k)/=sum(i,j,k);
         }
      }
   }

   grid.v.zero();
   sum.zero();
   for(p=0; p<np; ++p){
      grid.bary_x_centre(x[p][0], i, fx);
      grid.bary_y(x[p][1], vj, vfy);
      grid.bary_z_centre(x[p][2], k, fz);
      accumulate(grid.v, u[p][1] , i, vj, k, fx, vfy, fz);
      affineFix(grid.v, cy[p], i, vj, k, fx, vfy, fz);
   }
   for(k=0; k<grid.v.nz; ++k){
      for(j=0; j<grid.v.ny; ++j){ 
         for(i=0; i<grid.v.nx; ++i){
            if(sum(i,j,k)!=0) grid.v(i,j,k)/=sum(i,j,k);
         }
      }
   }

   grid.w.zero();
   sum.zero();
   for(p=0; p<np; ++p){
      grid.bary_x_centre(x[p][0], i, fx);
      grid.bary_y_centre(x[p][1], j, fy);
      grid.bary_z(x[p][2], wk, wfz);
      accumulate(grid.w, u[p][2] , i, j, wk, fx, fy, wfz);
      affineFix(grid.w, cz[p], i, j, wk, fx, fy, wfz);
   }
   for(k=0; k<grid.w.nz; ++k){
      for(j=0; j<grid.w.ny; ++j){ 
         for(i=0; i<grid.w.nx; ++i){
            if(sum(i,j,k)!=0) grid.w(i,j,k)/=sum(i,j,k);
         }
      }
   }
   // identify where particles are in grid
   grid.marker.zero();
   for(p=0; p<np; ++p){
      grid.bary_x(x[p][0], i, fx);
      grid.bary_y(x[p][1], j, fy);
      grid.bary_z(x[p][2], k, fz);
      grid.marker(i,j,k)=FLUIDCELL;
   }
}

/* this function computes c from the gradient of w and the velocity field from the grid. */
Vec3f Particles::
computeC(Array3f &ufield, int i, int j, int k, float fx, float fy, float fz) //ufield: grid.u, grid.v, grid.v
{
   Vec3f c = Vec3f(-1.0*(1.0 - fy)*(1.0 - fz), -1.0*(1.0 - fx)*(1.0 - fz), -1.0*(1.0 - fx)*(1.0 - fy)) * ufield(i,j,k)
           + Vec3f((1.0 - fy)*(1.0 - fz), -1.0*fx*(1.0 - fz), -1.0*fx*(1.0 - fy)) * ufield(i+1,j,k)
           + Vec3f(-1.0* fy*(1.0 - fz), (1.0 - fx)*(1.0 - fz), -1.0*fx*(1.0 - fy)) * ufield(i,j+1,k)
           + Vec3f(fy*(1.0 - fz), fx*(1.0 - fz), -1.0*fx*fy) * ufield(i+1,j+1,k)
           + Vec3f(-1.0*(1.0 - fy)*fz, -1.0*(1.0 - fx)*fz, (1.0 -fx)*(1.0 - fy)) * ufield(i,j,k+1)
           + Vec3f((1.0 - fy)*fz, -1.0*fx*fz, fx*(1.0 - fy)) * ufield(i+1, j, k+1)
           + Vec3f(-1.0*fy*fz, (1.0-fx)*fz, (1.0-fx)*fy*fz) * ufield(i, j+1, k+1)
           + Vec3f(fy*fz, fx*fz, fx*fy) * ufield(i+1, j+1, k+1);
   return c;
}

void Particles::
update_from_grid(void)
{
   int p;
   int i, ui, j, vj, k, wk;
   float fx, ufx, fy, vfy, fz, wfz;
   for(p=0; p<np; ++p){
      grid.bary_x(x[p][0], ui, ufx);
      grid.bary_x_centre(x[p][0], i, fx);
      grid.bary_y(x[p][1], vj, vfy);
      grid.bary_y_centre(x[p][1], j, fy);
      grid.bary_z(x[p][2], wk, wfz);
      grid.bary_z_centre(x[p][2], k, fz);
      if( simType == FLIP )
      {
         /*working on*/
         u[p]+=Vec3f(grid.du.trilerp(ui, j, k, ufx, fy, fz), grid.dv.trilerp(i, vj, k, fx, vfy, fz), grid.dw.trilerp(i, j, wk, fx, fy, wfz)); // FLIP
      }
      else
      {
         /*working on*/
         u[p]=Vec3f(grid.u.trilerp(ui, j, k, ufx, fy, fz), grid.v.trilerp(i, vj, k, fx, vfy, fz), grid.w.trilerp(i, j, wk, fx, fy, wfz)); // PIC and APIC

         if( simType == APIC )
         {
             cx[p] = computeC(grid.u, ui, j, k, ufx, fy, fz);
             cy[p] = computeC(grid.v, i, vj, k, fx, vfy, fz);
             cz[p] = computeC(grid.w, i, j, wk, fx, fy, wfz);
         }
      }
      
   }
}

void Particles::
move_particles_in_grid(float dt)
{
   Vec3f midx, gu;
   float xmin=1.001*grid.h, xmax=grid.lx-1.001*grid.h;
   float ymin=1.001*grid.h, ymax=grid.ly-1.001*grid.h;
   float zmin=1.001*grid.h, zmax=grid.lz-1.001*grid.h;
   for(int p=0; p<np; ++p){
       // first stage of Runge-Kutta 2 (do a half Euler step)
      grid.trilerp_uvw(x[p][0], x[p][1], x[p][2], gu[0], gu[1], gu[2]); // working
      midx=x[p]+0.5*dt*gu;
      clamp(midx[0], xmin, xmax);
      clamp(midx[1], ymin, ymax);
      clamp(midx[2], zmin, zmax);
      // second stage of Runge-Kutta 2
      grid.trilerp_uvw(midx[0], midx[1], midx[2], gu[0], gu[1], gu[2]); // wokring
      x[p]+=dt*gu;
      clamp(x[p][0], xmin, xmax);
      clamp(x[p][1], ymin, ymax);
      clamp(x[p][2], zmin, zmax);
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

   fprintf(fp, "%d %lf\n", np, 0.002);
   for(int p=0; p<np; ++p)
      fprintf(fp, "%.5g %.5g %.5g\n", x[p][0], x[p][1], x[p][2]);
   fclose(fp);
}

