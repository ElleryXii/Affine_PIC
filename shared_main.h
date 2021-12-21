/**
 * Edit this file to set up the initial conditions for the fluid simulation 
 * (what the water looks like, etc). 
 *
 * @author Ante Qu, 2017
 * Based on Bridson's simple_flip2d starter code at http://www.cs.ubc.ca/~rbridson/
 */

#ifndef SHARED_MAIN_H
#define SHARED_MAIN_H


#define SIMULATION_TYPE (APIC) // default simtype: APIC, FLIP, or PIC
#define INIT_DROP_RADIUS (0.05)
#define INIT_FLOOR_SIZE (0.05)
#define INIT_NA (4)
#define INIT_NB (4)
#define INIT_NC (4)
#define USE_SPHERICAL_GRAV (true) //not implemented for 3d
// the following only matter when USE_SPHERICAL_GRAV is true
#define INIT_VEL_MAGNITUDE (0.55)
#define GRAV_CENTER_X (0.5)
#define GRAV_CENTER_Y (0.5)
#define GRAV_CENTER_Z (0.5)
#define GRAV_FACTOR (0.01)

#define N_ITER (100)

using namespace std;
/* This sets the signed distance function phi of the fluid. You can selectively uncomment a line to decide */
/* which example to run. This is used in initializing the water. */
float fluidphi(Grid &grid, float x, float y, float z)
{
//   return y-0.05*grid.ly; // no drop
//   return min(sqrt(sqr(x-0.5*grid.lx)+sqr(y-0.625*grid.ly)+sqr(z-0.5*grid.lz))-0.02*grid.ly, y-0.6*grid.ly); // tiny drop
//   return min(sqrt(sqr(x-0.3333*grid.lx)+sqr(y-0.71*grid.ly)+sqr(z-0.3333*grid.lz))-0.3*grid.ly, y-0.2*grid.ly); // large drop
//   return max(y-0.8*grid.ly, -sqrt(sqr(x-0.5*grid.lx)+sqr(y-0.2*grid.ly)+sqr(z-0.5*grid.lz))+0.1*grid.lx); // bubble
   //return sqrt(sqr(x-0.5*grid.lx)+sqr(y-0.75*grid.ly))-0.15*grid.lx; // large drop without bottom
   return min(y-INIT_FLOOR_SIZE*grid.ly, sqrt(sqr(x-0.5*grid.lx)+sqr(y-0.7*grid.ly)+sqr(z-0.5*grid.lz))-INIT_DROP_RADIUS*grid.lx); // medium drop
//    return 0.95*grid.ly-y;; // horizontal
//    return 0.95*grid.lx-x; // vertical column
//   return max( max(x-0.75*grid.lx, 0.25*grid.lx-x), max(y-0.75*grid.ly, 0.25*grid.ly-y), max(z-0.75*grid.lz, 0.25*grid.lz-z)); // small box
}

/* Helper function for initializing the water */
void project(Grid &grid, float &x, float &y, float &z, float current, float target)
{
   float dpdx=(fluidphi(grid, x+1e-4, y,z)-fluidphi(grid, x-1e-4, y,z))/2e-4;
   float dpdy=(fluidphi(grid, x, y+1e-4, z)-fluidphi(grid, x, y-1e-4,z))/2e-4;
   float dpdz=(fluidphi(grid, x, y, z+1e-4)-fluidphi(grid, x, y, z-1e-4))/2e-4;
   float scale=(target-current)/sqrt(dpdx*dpdx+dpdy*dpdy+dpdz*dpdz);
   x+=scale*dpdx;
   y+=scale*dpdy;
   z+=scale*dpdz;
}

/* This function allocates particles based on the phi function. */
void init_water_drop(Grid &grid, Particles &particles, int na, int nb, int nc)
{
   int i, j, k, a, b, c;
   float x, y, z, phi;
   float vx, vy, vz;
   vx = 0;
   vy = 0;
   vz = 0;

   for(i=1; i<grid.marker.nx-1; ++i){
      for(j=1; j<grid.marker.ny-1; ++j){
          for(k=1; k<grid.marker.nz-1; ++k){
              for(a=0; a<na; ++a){
                  for(b=0; b<nb; ++b){
                      for (c=0; c<nc; ++c){
                          x=(i+(a+0.1+0.8*rand()/(double)RAND_MAX)/na)*grid.h;
                          y=(j+(b+0.1+0.8*rand()/(double)RAND_MAX)/nb)*grid.h;
                          z=(k+(c+0.1+0.8*rand()/(double)RAND_MAX)/nc)*grid.h;
                          if( USE_SPHERICAL_GRAV )
                          {
                              vx = INIT_VEL_MAGNITUDE * grid.lx * (rand() * 2.0 / (double)RAND_MAX + 0.5 );
                              vy = INIT_VEL_MAGNITUDE * grid.ly * (rand() * 2.0 / (double)RAND_MAX - 1.0 );
                              vz = INIT_VEL_MAGNITUDE * grid.lz * (rand() * 2.0 / (double)RAND_MAX - 1.0 );
                          }
                          phi=fluidphi(grid, x, y, z);
                          if(phi>-0.25*grid.h/na)
                              continue;
                          else if(phi>-1.5*grid.h/na){
                              project(grid, x, y, z, phi, -0.75*grid.h/na);
                              phi=fluidphi(grid, x, y,z);
                              project(grid, x, y, z, phi, -0.75*grid.h/na);
                              phi=fluidphi(grid, x, y,z);
                          }
                          particles.add_particle(Vec3f(x,y,z), Vec3f(vx,vy,vz));
                      }
                  }
              }
          }
      }
   }

}

void advance_one_step(Grid &grid, Particles &particles, double dt)
{

   for(int i=0; i<5; ++i)
      particles.move_particles_in_grid(0.2*dt);
   particles.transfer_to_grid();
   grid.save_velocities();
   grid.add_gravity(dt, USE_SPHERICAL_GRAV, GRAV_CENTER_X * grid.lx, GRAV_CENTER_Y  * grid.ly, GRAV_CENTER_Z *grid.lz);
   grid.compute_distance_to_fluid();
   grid.extend_velocity();
   grid.apply_boundary_conditions();
   grid.make_incompressible();
   grid.extend_velocity();
   grid.get_velocity_update();
   particles.update_from_grid();
}

void advance_one_frame(Grid &grid, Particles &particles, double frametime)
{
   double t=0;
   double dt;
   bool finished=false;
   while(!finished){
      dt=2*grid.CFL();
      if(t+dt>=frametime){
         dt=frametime-t;
         finished=true;
      }else if(t+1.5*dt>=frametime)
         dt=0.5*(frametime-t);
      printf("advancing %g (to %f%% of frame)\n", dt, 100.0*(t+dt)/frametime);
      advance_one_step(grid, particles, dt);
      t+=dt;
   }
}


#endif
