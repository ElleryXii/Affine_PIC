/**
 * Declaration for the Particles struct, which holds the particles.
 *
 * @author Ante Qu, 2017
 * Based on Bridson's simple_flip2d starter code at http://www.cs.ubc.ca/~rbridson/
 */

#ifndef PARTICLES_H
#define PARTICLES_H

#include <vector>
#include "grid.h"
#include "vec3.h"

typedef enum SimulationTypeEnum { PIC = 0, FLIP = 1, APIC = 2 } SimulationType;

struct Particles{
   Grid &grid;
   int np; // number of particles
   std::vector<Vec3f> x, u; // positions and velocities
   /* TODO: add helper variables */
   std::vector<Vec3f> cx, cy, cz; // c vectors stored, times h

   // transfer stuff
   Array3f sum;
   SimulationType simType;

   Particles(Grid &grid_, SimulationType simType_)
      :grid(grid_), np(0),
       sum(grid_.pressure.nx+1, grid_.pressure.ny+1, grid_.pressure.nz+1), simType( simType_ )
   {}

   void add_particle(const Vec3f &px, const Vec3f &pu, const Vec3f &pw);
   void transfer_to_grid(void);
   void update_from_grid(void);
   void move_particles_in_grid(float dt);
   void write_to_file(const char *filename_format, ...);

   private:
   template<class T> void accumulate(T &accum, float q, int i, int j, int k, float fx, float fy, float fz);
   template<class T> void affineFix(T &accum, Vec3f c, int i, int j, int k, float fx, float fy, float fz);
   Vec3f computeC(Array3f &ufield, int i, int j, int k, float fx, float fy, float fz);
};

#endif
