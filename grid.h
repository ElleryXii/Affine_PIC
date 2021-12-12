/**
 * Definition of the Grid struct, which holds the MAC grid.
 *
 * @author Ante Qu, 2017
 * Based on Bridson's simple_flip2d starter code at http://www.cs.ubc.ca/~rbridson/
 */

#ifndef GRID_H
#define GRID_H

#include "array3.h"
#include "array3.h"
#include "util.h"

#define AIRCELL 0
#define FLUIDCELL 1
#define SOLIDCELL 2

struct Grid{
   float gravity;
   float lx, ly, lz; // width，height, depth
   float h, overh; //cell side length, 1/cell side length

   // active variables
   Array3f u, v, w; // staggered MAC grid of velocities
   Array3f du, dv, dw; // saved velocities and differences for particle update
   Array3c marker; // identifies what sort of cell we have
   Array3f phi; // decays away from water into air (used for extrapolating velocity)
   Array3d pressure;
   // stuff for the pressure solve
   Array3x3f poisson;
   Array3d preconditioner;
   Array3d m;
   Array3d r, z, s;

   Grid(void)
   {}

   Grid(float gravity_, int cell_nx, int cell_ny, int cell_nz, float lx_)
   { init(gravity_, cell_nx, cell_ny, cell_nz, lx_); }

   void init(float gravity_, int cell_nx, int cell_ny, int cell_nz, float lx_);
   float CFL(void);
   void save_velocities(void);
   void add_gravity(float dt, bool centered, float cx, float cy, float cz);
   void compute_distance_to_fluid(void);
   void extend_velocity(void);
   void apply_boundary_conditions(void);
   void make_incompressible(void);
   void get_velocity_update(void);

   //TODO: unify bary functions


   void bary_x(float x, int &i, float &fx)
   {
      float sx=x*overh;
      i=(int)sx;
      fx=sx-floor(sx);
   }

   void bary_x_centre(float x, int &i, float &fx)
   {
      float sx=x*overh-0.5;
      i=(int)sx;
      if(i<0){ i=0; fx=0.0; }
      else if(i>pressure.nx-2){ i=pressure.nx-2; fx=1.0; }
      else{ fx=sx-floor(sx); }
   }

   void bary_y(float y, int &j, float &fy)
   {
      float sy=y*overh;
      j=(int)sy;
      fy=sy-floor(sy);
   }

   void bary_y_centre(float y, int &j, float &fy)
   {
      float sy=y*overh-0.5;
      j=(int)sy;
      if(j<0){ j=0; fy=0.0; }
      else if(j>pressure.ny-2){ j=pressure.ny-2; fy=1.0; }
      else{ fy=sy-floor(sy); }
   }

   void bilerp_uv(float px, float py, float &pu, float &pv)
   {
      int i, j;
      float fx, fy;
      bary_x(px, i, fx);
      bary_y_centre(py, j, fy);
      pu=u.bilerp(i, j, fx, fy);
      bary_x_centre(px, i, fx);
      bary_y(py, j, fy);
      pv=v.bilerp(i, j, fx, fy);
   }

   //TODO
   void trilerp_uvw(float px, float py, float pz, float &pu, float &pv, float &pw){

   }

   private:
   void init_phi(void);
   void sweep_phi(void);
   void sweep_u(int i0, int i1, int j0, int j1);
   void sweep_v(int i0, int i1, int j0, int j1);
   void sweep_w(int i0, int i1, int j0, int j1);
   void sweep_velocity(void);
   void find_divergence(void);
   void form_poisson(void);
   void form_preconditioner(void);
   void apply_poisson(const Array3d &x, Array3d &y);
   void apply_preconditioner(const Array3d &x, Array3d &y, Array3d &temp);
   void solve_pressure(int maxits, double tolerance);
   void add_gradient(void);
};

#endif

