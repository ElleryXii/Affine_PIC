/**
 * File for setting up the main simulator
 *
 * @author Ante Qu, 2017
 * Based on Bridson's simple_flip2d starter code at http://www.cs.ubc.ca/~rbridson/
 */

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include "particles.h"
#include "util.h"
#include "shared_main.h"

#define N_ITER (100)
using namespace std;


int main(int argc, char **argv)
{
   float gravity = 9.8;
   if( USE_SPHERICAL_GRAV )
      gravity *= GRAV_FACTOR;
   Grid grid(gravity, 50, 50, 1);
   SimulationType sType = SIMULATION_TYPE;
   
   std::string outputpath=".";

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
   Particles particles(grid, sType);

   init_water_drop(grid, particles, 2, 2);
   particles.write_to_file("%s/frameparticles%04d", outputpath.c_str(), 0);

   for(int i=1; i<N_ITER + 1; ++i){
      printf("===================================================> step %d...\n", i);
      advance_one_frame(grid, particles, 1./30);
      particles.write_to_file("%s/frameparticles%04d", outputpath.c_str(), i);
   }

   return 0;
}
