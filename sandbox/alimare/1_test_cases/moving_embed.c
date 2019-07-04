/**
#Double embed boundary calculation

First test case with a double calculation inside and outside of an embedded boundary.
This is based on the example of the Bénard–von Kármán Vortex Street for flow 
around a cylinder at Re=160, I modified the geometry for fun and use a
Zalesak's notched disk. Both fluids have the same properties. We use cs for
the cell fraction in one phase and 1-cs in the other.

At the end of each iteration of the solver we swap the variables, change the fs,
the cs and do a solver iteration on the new set of variables.
We use the centered Navier-Stokes solver, with embedded boundaries and advect
the passive tracer *f*. 

I still have issues regarding the use of the same timestep for both phases. They
hare truly independant. Once this is done I will work on making both phases
fields dependant of one another and then make the boundary movement dependant on
both fields.


![Animation of cs*u.x + (1-cs)*u2.x.](update/movie.mp4)(loop)
*/

#define DOUBLE_EMBED  0
#define LevelSet      1

#include "embed.h"
#include "../centered_alex.h"
#include "../level_set.h"
#include "diffusion.h"
#include "tracer.h"
#include "view.h"

#define MIN_LEVEL 3
#define MAXLEVEL 8
#define T_eq         0.
#define TL_inf       1.


#define H0 -0.99*L0/4.
#define DT_MAX  L0/(1 << MAXLEVEL)*0.8

#define T_eq         0.
#define plane(x, y, H) (y - H)


scalar TL[], dist[];
scalar * tracers = {TL};
scalar * level_set = {dist};

face vector muv[];
mgstats mgT;

int     nb_cell_NB =  1 << 3 ;  // number of cells for the NB
double  NB_width ;              // length of the NB
  

/**
The domain is 4 units long, centered vertically. */

int main() {
  L0 = 4.;
  CFL = 0.5;
  origin (-L0/2., -L0/2.);
  
  N = 1 << MAXLEVEL;
  NB_width = L0*nb_cell_NB / (1<< MAXLEVEL);
  mu = muv;
  run(); 
}

/**
We set a constant viscosity corresponding to a Reynolds number of 160,
based on the cylinder diameter (0.125) and the inflow velocity (1). */

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*0.125/160.;
  boundary((scalar *) {muv});
}

TL[embed]  = dirichlet(T_eq);
TL[top]    = dirichlet(TL_inf); 

event init (t = 0)
{
  DT = 5.e-4;
  /**
  The domain is the intersection of a channel of width unity and a
  circle of diameter 0.125. */

#if TREE
    refine (level < MAXLEVEL && plane(x, y, (H0 - NB_width)) > 0.
            && plane(x, y, (H0 + NB_width)) < 0.);
#endif

  foreach_vertex(){
    dist[] = plane(x,y,H0);
  }
  boundary ({dist});
  restriction({dist});
  fractions (dist, cs, fs);

  boundary({cs,fs});
  restriction({cs,fs});
  
  foreach() {
    TL[] = cs[]*TL_inf;
  }

  boundary({TL});
  restriction({TL});
}

event stability(i++){
    double dtmax2 = DT_MAX;
    timestep (uf, dtmax2);
}

event tracer_diffusion(i++){

  foreach_face()
    muv.x[] = fm.x[]*0.125/160.;
  boundary((scalar *) {muv});
  mgT = diffusion(TL, dt, D = 1.);

}

/**
We produce an animation of the tracer field. */

event movies ( i +=10,last)
{
  boundary({TL});
  view (fov = 16.642, quat = {0,0,0,1}, tx = -0.0665815, 
    ty = -0.00665815, bg = {1,1,1}, width = 600, height = 600, samples = 1);
  // draw_vof("cs");
  squares("TL", min=0, max = 1);
  save ("TL.png");
}

/**
We check the number of iterations of the Poisson and viscous
problems. */

event logfile (i++;i<=20){

  fprintf (stderr, "%d %g %g %d %d\n", i, t, dt, mgp.i, mgu.i);
}

