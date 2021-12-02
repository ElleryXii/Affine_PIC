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
#include "vec2.h"

typedef enum SimulationTypeEnum { PIC = 0, FLIP = 1, APIC = 2 } SimulationType;

struct Particles{
   Grid &grid;
   int np; // number of particles
   std::vector<Vec2f> x, u; // positions and velocities
   /* TODO: add helper variables */
   //std::vector<Vec2f> cx, cy; // c vectors stored, times h

   // transfer stuff
   Array2f sum;
   SimulationType simType;

   Particles(Grid &grid_, SimulationType simType_)
      :grid(grid_), np(0),
       sum(grid_.pressure.nx+1, grid_.pressure.ny+1), simType( simType_ )
   {}

   void add_particle(const Vec2f &px, const Vec2f &pu);
   void transfer_to_grid(void);
   void update_from_grid(void);
   void move_particles_in_grid(float dt);
   void write_to_file(const char *filename_format, ...);

   private:
   template<class T> void accumulate(T &accum, float q, int i, int j, float fx, float fy);
   template<class T> void affineFix(T &accum, Vec2f c, int i, int j, float fx, float fy);
   Vec2f computeC(Array2f &ufield, int i, int j, float fx, float fy);
};

#endif
