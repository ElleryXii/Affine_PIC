/**
 * Implementation of the MAC grid. Most of your edits should be in particles.cpp instead.
 *
 * @author Ante Qu, 2017
 * Based on Bridson's simple_flip2d starter code at http://www.cs.ubc.ca/~rbridson/
 */

#include <cmath>
#include <vector>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include "grid.h"

using namespace std;

void Grid::
init(float gravity_, int cell_nx, int cell_ny, int cell_nz, float lx_)
{
   gravity=gravity_;
   lx = lx_;
   ly = cell_ny*lx/cell_nx;
   lz = cell_nz*lx/cell_nx;
   h=lx/cell_nx;
   overh=cell_nx/lx;
   // allocate all the grid variables
   u.init(cell_nx+1, cell_ny, cell_nz);
   v.init(cell_nx, cell_ny+1, cell_nz);
   w.init(cell_nx, cell_ny, cell_nz+1);
   pressure.init(cell_nx, cell_ny, cell_nz);
   marker.init(cell_nx, cell_ny, cell_nz);
   phi.init(cell_nx, cell_ny, cell_nz);
   du.init(cell_nx+1, cell_ny, cell_nz);
   dv.init(cell_nx, cell_ny+1, cell_nz);
   dw.init(cell_nx, cell_ny, cell_nz+1)
   poisson.init(cell_nx, cell_ny, cell_nz);
   preconditioner.init(cell_nx, cell_ny, cell_nz);
   m.init(cell_nx, cell_ny, cell_nz);
   r.init(cell_nx, cell_ny, cell_nz);
   z.init(cell_nx, cell_ny, cell_nz);
   s.init(cell_nx, cell_ny, cell_ny);
}

//TOTEST
float Grid::
CFL(void)
{
   float maxv2=max(h*gravity, sqr(u.infnorm())+sqr(v.infnorm())+ sqrt(w.infnorm()));
   if(maxv2<1e-16) maxv2=1e-16;
   return h/sqrt(maxv2);
}

void Grid::
save_velocities(void)
{
   u.copy_to(du);
   v.copy_to(dv);
   w.copy_to(dw);
}

//TODO: OPTIONAL spherical gravity
/* centered gravity is the spherical gravity I added. */
/* for the normal uniform gravity, see lines 107-108 */
void Grid::
add_gravity(float dt, bool centered, float cx, float cy)
{
   float dtg=dt*gravity;
   if( centered )
   {
      for(int i = 0; i < u.nx; ++i)
      {
         for( int j = 0; j < v.ny; ++j)
         {
            float x, y, dx, dy, dr2, dr;
            if ( j < u.ny)
            {
               x = ( i ) * h;
               y = ( 0.5 + j ) * h;
               dx = x - cx;
               dy = y - cy;
               dr2 = dx * dx + dy * dy;
               if (dr2 < 0.0001 * h * h){
                  printf("dr2 too small: %f \n", dr2);
                  dr2 = 0.0001 * h * h;
               }
               dr = sqrt(dr2);
               dx /= dr;
               dy /= dr;
               u(i,j) -= dtg * dx / dr2;
            }
            if( i < v.nx )
            {
               x = ( 0.5 + i ) * h;
               y = ( j ) * h;
               dx = x - cx;
               dy = y - cy;
               dr2 = dx * dx + dy * dy;
               if (dr2 < 0.0001 * h * h){
                  printf("dr2 too small: %f \n", dr2);
                  dr2 = 0.0001 * h * h;
               }
               dr = sqrt(dr2);
               dx /= dr;
               dy /= dr;
               v(i,j) -= dtg * dy / dr2;
            }
         }
      }
   }
   else
   {
      for(int i=0; i<v.size; ++i)
         v.data[i]-=dtg;
   }
}

//TODO: # of itr
void Grid::
compute_distance_to_fluid(void)
{
   init_phi();
   for(int i=0; i<2; ++i)
      sweep_phi();
}

//TODO: # of itr
void Grid::
extend_velocity(void)
{
   for(int i=0; i<4; ++i)
      sweep_velocity();
}

//TOTEST
void Grid::
apply_boundary_conditions(void)
{
    int i, j, k;
    // first mark where solid is
    for (k=0; k<marker.nz; ++k)
        marker(0, 0, k) = marker(marker.nx-1,marker.ny-1 ,k) = SOLIDCELL;
    for (j=0; j<marker.ny; ++j)
        marker(0,j, 0)=marker(marker.nx-1, j, marker.nz-1)=SOLIDCELL;
    for (i=0; i<marker.nx; ++i)
        marker(i,0, 0)=marker(i,marker.ny-1,marker.nz-1)=SOLIDCELL;

    // now makre sure nothing leaves the domain
    for(k=0; k<w.nz; ++k)
        w(0,0,k)=w(1,1,k)=w(w.nx-1, w.ny-1,k)=w(w.nx-2,w.ny-2,j)=0;
    for(j=0; j<u.ny; ++j)
        u(0,j,0)=u(1,j,1)=u(u.nx-1,j,u.nz-1)=u(u.nx-2,j,u.nz-2)=0;
    for(i=0; i<v.nx; ++i)
        v(i,0,0)=v(i,1,1)=v(i,v.ny-1,v.nz-1)=v(i,v.ny-2,v.nz-2)=0;
}

void Grid::
make_incompressible(void)
{
   find_divergence();
   form_poisson();
   form_preconditioner();
   solve_pressure(100, 1e-5);
   add_gradient();
}

//TOTEST
void Grid::
get_velocity_update(void)
{
   int i;
   for(i=0; i<u.size; ++i)
       du.data[i]=u.data[i]-du.data[i];
   for(i=0; i<v.size; ++i)
       dv.data[i]=v.data[i]-dv.data[i];
   for(i=0; i<w.size; ++i)
       dw.data[i]=w.data[i]-dw.data[i];
}

//====================================== private helper functions ============================

//TOTEST
void Grid::
init_phi(void)
{
   int i, j, k;
   // start off with indicator inside the fluid and overestimates of distance outside
   float large_distance=phi.nx+phi.ny+phi.nz+2;
   for(i=0; i<phi.size; ++i)
      phi.data[i]=large_distance;
   for(k = 1; k<phi.nz-1; ++k)for(j=1; j<phi.ny-1; ++j)for(i=1; i<phi.nx-1; ++i){
      if(marker(i,j,k)==FLUIDCELL){
         phi(i,j,k)=-0.5;
      }
   }
}

// 3D Eikonal equation
// Reference: Fluid Simulation for Computer Graphics 2nd Edition by Robert Bridson, Chapter 4 level set geometry, figure 4.4
//TOTEST
static inline void solve_distance(float phi0, float phi1, float phi2, float &r)
{
    //sort, ph0<phi1<phi2
    if (phi0 > phi1) std::swap(phi0,phi1);
    if (phi1 > phi2) std::swap(phi1,phi2);
    if (phi0 > phi1) std::swap(phi0,phi1);
    float d = phi0+1;
    if (d > phi1)
        d = 0.5*(phi0+phi1+sqrt(2-sqr(phi1-phi0)));
    if (d > phi2)
        d = (phi0+phi1+phi2 + sqrt(max(0, sqr(phi0+phi1+phi2)-3*(sqr(phi0)+ sqr(phi1)+sqr(phi2)-1))))/3.0;
    if(d<r) r=d;
}

inline void Grid::
sweep_phi(int i0, int i1, int j0, int j1, int k0, int k1)
{
    int di=(i0<i1) ? 1 : -1, dj=(j0<j1) ? 1 : -1, dk = (k0<k1)?1:-1;
    for (int k=k0; k!=k1; k+=dk) for( int j=j0; j!=j1; j+=dj) for(int i=i0; i!=i1; i+=di)
    {
        if (marker(i,j,k)!=FLUIDCELL)
            solve_distance(phi(i-di, j, k),phi(i, j-dj, k) , phi(i, j, k-dk),phi(i,j,k));
    }
}

//TOTEST
void Grid::
sweep_phi(void)
{
   // fast sweeping outside the fluid in all four sweep directions
    sweep_phi(1, phi.nx, 1, phi.ny, 1, phi.nz);
    sweep_phi(1, phi.nx, 1, phi.ny, phi.nz-2, -1);
    sweep_phi(1, phi.nx, phi.ny-2, -1, 1, phi.nz);
    sweep_phi(phi.nx-2, -1, 1, phi.ny, 1, phi.nz);
    sweep_phi(1, phi.nx, phi.ny-2, -1,  phi.nz-2, -1);
    sweep_phi(phi.nx-2, -1, 1, phi.ny, phi.nz-2, -1);
    sweep_phi(phi.nx-2, -1, phi.ny-2, -1, 1, phi.nz);
    sweep_phi(phi.nx-2, -1, phi.ny-2, -1, phi.nz-2, -1);
}

//TODO
void Grid::
sweep_u(int i0, int i1, int j0, int j1)
{
   int di=(i0<i1) ? 1 : -1, dj=(j0<j1) ? 1 : -1;
   float dp, dq, alpha;
   for(int j=j0; j!=j1; j+=dj) for(int i=i0; i!=i1; i+=di)
      if(marker(i-1,j)==AIRCELL && marker(i,j)==AIRCELL){
         dp=di*(phi(i,j)-phi(i-1,j));
         if(dp<0) continue; // not useful on this sweep direction
         dq=0.5*(phi(i-1,j)+phi(i,j)-phi(i-1,j-dj)-phi(i,j-dj));
         if(dq<0) continue; // not useful on this sweep direction
         if(dp+dq==0) alpha=0.5;
         else alpha=dp/(dp+dq);
         u(i,j)=alpha*u(i-di,j)+(1-alpha)*u(i,j-dj);
      }
}
//TODO
void Grid::
sweep_v(int i0, int i1, int j0, int j1)
{
   int di=(i0<i1) ? 1 : -1, dj=(j0<j1) ? 1 : -1;
   float dp, dq, alpha;
   for(int j=j0; j!=j1; j+=dj) for(int i=i0; i!=i1; i+=di)
      if(marker(i,j-1)==AIRCELL && marker(i,j)==AIRCELL){
         dq=dj*(phi(i,j)-phi(i,j-1));
         if(dq<0) continue; // not useful on this sweep direction
         dp=0.5*(phi(i,j-1)+phi(i,j)-phi(i-di,j-1)-phi(i-di,j));
         if(dp<0) continue; // not useful on this sweep direction
         if(dp+dq==0) alpha=0.5;
         else alpha=dp/(dp+dq);
         v(i,j)=alpha*v(i-di,j)+(1-alpha)*v(i,j-dj);
      }
}

//TODO, combine sweeping
void Grid::
sweep_w(int i0, int i1, int j0, int j1, int k0, int k1)
{
    int di=(i0<i1) ? 1 : -1, dj=(j0<j1) ? 1 : -1, dk = (k0<k1)?1:-1;
    float dp, dq, dr, alpha;
    for (int k = k0; k!=k1; k+=dk) for(int j=j0; j!=j1; j+=dj) for(int i=i0; i!=i1; i+=di)
        if(marker(i,j, k-1)==AIRCELL && marker(i,j,k)==AIRCELL)
        {
            dq=dj*(phi(i,j,k)-phi(i,j,k-1));
            if(dq<0) continue; // not useful on this sweep direction
            dp=0.5*(phi(i,j-1)+phi(i,j)-phi(i-di,j-1)-phi(i-di,j));
            if(dp<0) continue; // not useful on this sweep direction

            if(dr<0) continue;

            if(dp+dq==0) alpha=0.5;
            else alpha=dp/(dp+dq);
            w(i,j)=alpha*v(i-di,j)+(1-alpha)*v(i,j-dj);
        }
}


//TOTEST
static inline void sweep_velocity_boundary(Array3f& vfield)
{
    int i, j, k;

    for (i=0; i<vfield.nx; ++i){
        for (j = 0; j<vfield.ny; ++j){
            vfield(i,j, 0) = vfield(i,j,1);
            vfield(i,j, vfield.nz-1) = vfield(i,j,vfield.nz-2);
        }
    }

    for (i=0; i<vfield.nx; ++i){
        for (k=0; k<vfield.nz; ++k){
            vfield(i,0,k) = vfield(i,1,k);
            vfield(i,vfield.ny-1,k) = vfield(i,vfield.ny-2,k);
        }
    }

    for (j=0;j<vfield.ny;++j){
        for (k=0;k<vfield.nz;++k){
            vfield(0,j,k) = vfield(1,j,k);
            vfield(vfield.nx-1,j,k) = vfield(vfield.nx-2,j,k);
        }
    }
}

//TODO: add directions, add w, simplify?
void Grid::
sweep_velocity(void)
{
    // sweep u, only into the air
    sweep_u(1, u.nx-1, 1, u.ny-1);
    sweep_u(1, u.nx-1, u.ny-2, 0);
    sweep_u(u.nx-2, 0, 1, u.ny-1);
    sweep_u(u.nx-2, 0, u.ny-2, 0);
    sweep_velocity_boundary(u);

    // now the same for v
    sweep_v(1, v.nx-1, 1, v.ny-1);
    sweep_v(1, v.nx-1, v.ny-2, 0);
    sweep_v(v.nx-2, 0, 1, v.ny-1);
    sweep_v(v.nx-2, 0, v.ny-2, 0);
    sweep_velocity_boundary(v);

    // for w
    sweep_w(1, w.nx-1, 1, w.ny-1);
    sweep_w(1, w.nx-1, w.ny-2, 0);
    sweep_w(w.nx-2, 0, 1, w.ny-1);
    sweep_w(w.nx-2, 0, w.ny-2, 0);
    sweep_velocity_boundary(w);
}

//TOTEST
void Grid::
find_divergence(void)
{
   r.zero();
   for(int k = 0; k<r.nz; ++k)for(int j=0; j<r.ny; ++j) for(int i=0; i<r.nx; ++i){
      if(marker(i,j, k)==FLUIDCELL)
         r(i,j, k)=u(i+1,j,k)-u(i,j, k)+v(i,j+1,k)-v(i,j)+w(i,j,k+1)-w(i,j,k);
   }
}

//TOTEST check # of direction; change poisson to Array3x4
void Grid::
form_poisson(void)
{
   poisson.zero();
   for(int k = 1; k<poisson.nz-1; ++k)for(int j=1; j<poisson.ny-1; ++j) for(int i=1; i<poisson.nx-1; ++i){
      if(marker(i,j,k)==FLUIDCELL){

         if(marker(i-1,j,k)!=SOLIDCELL)
            poisson(i,j,k,0)+=1;
         if(marker(i+1,j, k)!=SOLIDCELL){
            poisson(i,j,k,0)+=1;
            if(marker(i+1,j,k)==FLUIDCELL)
               poisson(i,j,k,1)=-1;
         }

         if(marker(i,j-1,k)!=SOLIDCELL)
            poisson(i,j,k,0)+=1;
         if(marker(i,j+1,k)!=SOLIDCELL){
            poisson(i,j,k,0)+=1;
            if(marker(i,j+1,k)==FLUIDCELL)
               poisson(i,j,k,2)=-1;
         }

          if(marker(i,j,k-1)!=SOLIDCELL)
              poisson(i,j,k,0)+=1;
          if(marker(i,j,k+1)!=SOLIDCELL){
              poisson(i,j,k,0)+=1;
              if(marker(i,j,k+1)==FLUIDCELL)
                  poisson(i,j,k,3)=-1;
          }
      }
   }
}

//TOTEST
void Grid::
apply_poisson(const Array3d &x, Array3d &y)
{
    y.zero();
    for (int k=1; k<poisson.nz-1; ++k) for(int j=1; j<poisson.ny-1; ++j) for(int i=1; i<poisson.nx-1; ++i){
        if(marker(i,j,k)==FLUIDCELL){
            y(i,j,k) = poisson(i,j,k,0) * x(i,j,k)
                            + poisson(i-1,j,k,1) * x(i-1,j,k)
                            + poisson(i,j,k,1) * x(i+1,j,k)
                            + poisson(i,j-1,k,2) * x(i,j-1,k)
                            + poisson(i,j,k,2) * x(i,j+1,k)
                            + poisson(i,j,k-1,3) * x(i,j,k-1)
                            + poisson(i,j,k,3) * x(i,j,k+1)
        }
    }
}

//TOTEST blind
void Grid::
form_preconditioner()
{
    const double mic_parameter=0.99;
    double d;
    preconditioner.zero();
    for (int k=1; k<preconditioner.nz-1; ++k) for(int j=1; j<preconditioner.ny-1; ++j) for(int i=1; i<preconditioner.nx-1; ++i){
        if(marker(i,j,k)==FLUIDCELL){
            d = poisson(i,j,k,0)
                - sqr( poisson(i-1,j, k,1) * preconditioner(i-1,j, k))
                - sqr( poisson(i,j-1, k,2) * preconditioner(i,j-1, k))
                - sqr( poisson(i,j,k-1,3) * preconditioner(i,j,k-1))
                - mic_parameter * (poisson(i-1,j,k,1) * poisson(i-1,j,k,2) * poisson(i-1,j,k,3) * sqr(preconditioner(i-1,j, k))
                                    + poisson(i,j-1,k,2) * poisson(i,j-1,k,1) * poisson(i,j-1,k,3) * sqr(preconditioner(i,j-1, k))
                                    + poisson(i,j, k-1,3) * poisson(i,j,k-1,2)* poisson(i,j,k-1,1)  * sqr(preconditioner(i,j,k-1)));
            preconditioner(i,j,k)=1.0/sqrt(d+1e-6);
        }
    }
}å

//TOTEST blind
void Grid::
apply_preconditioner(const Array3d &x, Array3d &y, Array3d &m)
{
   int i, j, k;
   float d;
   m.zero();
    // solve L*m=x
    for(k=1; k<x.nz-1; ++k) for(j=1; j<x.ny-1; ++j) for(i=1; i<x.nx-1; ++i)
        if(marker(i,j,k)==FLUIDCELL)
        {
            d = x(i,j, k)
                - poisson(i-1,j,k, 1)*preconditioner(i-1,j,k)*m(i-1,j,k)
                - poisson(i,j-1,k,2)*preconditioner(i,j-1, k)*m(i,j-1, k)
                - poisson(i,j,k-1,3)*preconditioner(i,j, k-1)*m(i,j, k-1);
            m(i,j, k) = preconditioner(i,j, k)*d;
        }

    // solve L'*y=m
    y.zero();
    for(k=x.nz-2; k>0; --k) for(j=x.ny-2; j>0; --j) for(i=x.nx-2; i>0; --i)
        if(marker(i,j)==FLUIDCELL)
        {
            d = m(i,j,k)
                - poisson(i,j,k,1)*preconditioner(i,j,k)*y(i+1,j,k)
                - poisson(i,j,k,2)*preconditioner(i,j,k)*y(i,j+1, k)
                - poisson(i,j,k,3)*preconditioner(i,j,k)*y(i,j, k+1);
            y(i,j, k) = preconditioner(i,j, k)*d;
        }
}

//TOTEST
void Grid::
solve_pressure(int maxits, double tolerance)
{
   int its;
   double tol=tolerance*r.infnorm();
   pressure.zero();
   if(r.infnorm()==0)
      return;
   apply_preconditioner(r, z, m);
   z.copy_to(s);
   double rho=z.dot(r);
   if(rho==0)
      return;
   for(its=0; its<maxits; ++its){
      apply_poisson(s, z);
      double alpha=rho/s.dot(z);
      pressure.increment(alpha, s);
      r.increment(-alpha, z);
      if(r.infnorm()<=tol){
         printf("pressure converged to %g in %d iterations\n", r.infnorm(), its);
         return;
      }
      apply_preconditioner(r, z, m);
      double rhonew=z.dot(r);
      double beta=rhonew/rho;
      s.scale_and_increment(beta, z);
      rho=rhonew;
   }
   printf("Didn't converge in pressure solve (its=%d, tol=%g, |r|=%g)\n", its, tol, r.infnorm());
}

//TOTEST
void Grid::
add_gradient(void)
{
   int i, j, k;
   for(k=1; k<u.nz-1;++k) for(j=1; j<u.ny-1; ++j) for(i=2; i<u.nx-2; ++i){
      if(marker(i-1,j, k)|marker(i,j,k)==FLUIDCELL){ // if at least one is FLUID, neither is SOLID
         u(i,j,k)+=pressure(i,j, k)-pressure(i-1,j,k);
      }
   }

    for(k=1; k<v.nz-1;++k) for(j=2; j<v.ny-2; ++j) for(i=1; i<v.nx-1; ++i){
      if(marker(i,j-1, k)|marker(i,j, k)==FLUIDCELL){ // if at least one is FLUID, neither is SOLID
         v(i,j,k)+=pressure(i,j, k)-pressure(i,j-1, k);
      }
   }

    for(k=2; k<w.nz-2;++k) for(j=1; j<w.ny-1; ++j) for(i=1; i<w.nx-1; ++i){
        if(marker(i,j, k-1)|marker(i,j,k)==FLUIDCELL){ // if at least one is FLUID, neither is SOLID
            w(i,j, k)+=pressure(i,j, k)-pressure(i,j, k-1);
        }
    }
}

