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

#define DOUBLE_EMBED  1
#define LevelSet      1

#include "embed.h"
#include "../centered_alex.h"
#include "../level_set.h"
#include "diffusion.h"
#include "tracer.h"
#include "view.h"

#define MIN_LEVEL 3
#define MAXLEVEL 7
#define latent_heat 10.
#define T_eq         0.
#define TL_inf       1.
#define TS_inf       -1.


#define H0 -0.99*L0/4.
#define DT_MAX  L0/(1 << MAXLEVEL)*0.8

#define T_eq         0.
#define plane(x, y, H) (y - H)



scalar TL[], TS[], dist[];
scalar * tracers = {TL, TS};
scalar * level_set = {dist};
face vector v_pc[];
face vector muv[];
mgstats mgT;

int     nb_cell_NB =  1 << 3 ;  // number of cells for the NB
double  NB_width ;              // length of the NB
  

TL[embed]  = dirichlet(T_eq);
TL[top]    = dirichlet(TL_inf); 

TS[embed]  = dirichlet(T_eq);
TS[bottom] = dirichlet(TS_inf); 

void phase_change_velocity_LS_embed (scalar cs, face vector fs, scalar tr,
 scalar tr2, face vector v_pc, scalar dist, double L_H, 
 double NB_width) {
  
 /**
  The phase change velocity $\mathbf{v}_{pc}$ is

  $$
  \mathbf{v}_{pc} = \mathrm{Pe}\, D\, \nabla tr
  $$
  
  here we use the embed_gradient_face_x defined in embed that gives a proper
  definition of the gradients with embedded boundaries. */
  
  face vector gtr[], gtr2[];
  foreach_face()
    if(fabs(dist[])<NB_width)
    gtr.x[] = face_gradient_x(tr,0);
  boundary((scalar*){gtr});
  foreach(){
      cs[]      = 1.-cs[];
  }
  foreach_face()
    fs.x[]      = 1.-fs.x[];

  boundary({cs,fs});
  restriction({cs,fs});
  
  foreach_face()
    if(fabs(dist[])<NB_width)
    gtr2.x[] = face_gradient_x(tr2,0);
  boundary((scalar*){gtr2});
  
  foreach(){
      cs[]      = 1.-cs[];
  }
  foreach_face()
    fs.x[]      = 1.-fs.x[];
  boundary({cs,fs});
  restriction({cs,fs});

  /**
  With the the normal vector and the gradients of the tracers we can now compute
  the phase change velocity $\mathbf{v}_{pc}$, following the lines drawn in
  [meanflow.c](/sandbox/popinet/meanflow.c). We define it as the product between
  the density ratio and the diffusive flow. Note that $\mathbf{v}_{pc}$ is weighted
  by the face metric. */
  
  foreach_face() {

    v_pc.x[] = 0.;

    if (fabs(dist[])<NB_width){

      coord n = facet_normal (point, cs, fs);
      normalize(&n);

      v_pc.x[] = (gtr.x[]*n.x*fs.x[] + gtr2.x[]*n.x*(1.-fs.x[]))*1e-1;
    }
  }
  boundary((scalar *){v_pc});
  restriction((scalar *){v_pc});
}

void advection_LS (struct Advection p) //never takes into account the solid
// boundary, the velocity field is defined on both phases.
{

  /**
  If *src* is not provided we set all the source terms to zero. */
  
  scalar * lsrc = p.src;
  if (!lsrc) {
    const scalar zero[] = 0.;
    for (scalar s in p.tracers)
      lsrc = list_append (lsrc, zero);
  }

  assert (list_len(p.tracers) == list_len(lsrc));
  scalar f, src;
  for (f,src in p.tracers,lsrc) {
    face vector flux[];
    tracer_fluxes (f, p.u, flux, p.dt, src);
    foreach()
      foreach_dimension()
        f[] += p.dt*(flux.x[] - flux.x[1])/(Delta); // careful we have removed
        // cm[]
  }
  boundary (p.tracers);

  if (!p.src)
    free (lsrc);
}


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
    TS[] = (1. - cs[])*TS_inf;
  }

  boundary({TL,TS});
  restriction({TL,TS});
}

event stability(i++){
    double dtmax2 = DT_MAX;
    timestep (uf, dtmax2);
}


event tracer_advection(i++,last){
  if(i%2 == 1){
  double L_H       = latent_heat;
  foreach_face(){
    v_pc.x[] = 0.;                        
    fs.x[]   = 1.-fs.x[];
  }
  foreach()
    cs[] = 1.-cs[];

  
  boundary({cs,fs});
  restriction({cs,fs});
  /**
  cs and fs have just been inverted, we create temporary fields just for
  calculation purposes

  When this calculation is done (1-cs) => TL, cs => TS (same for fs,(1-fs) )*/
  
  phase_change_velocity_LS_embed (cs, fs ,TL, TS, v_pc, dist, L_H, NB_width);
  
  foreach_face(){
    fs.x[]  = 1.-fs.x[];
  }
  foreach()
    cs[] = 1.-cs[];

  boundary({cs,fs});
  restriction({cs,fs});

  advection_LS (level_set, v_pc, dt);
  boundary ({dist});
  restriction({dist});

  // fractions (dist, cs, fs);
  // boundary({cs,fs});
  // restriction({cs,fs});

  event ("properties");
  }
}

event tracer_diffusion(i++){

  foreach_face()
    muv.x[] = fm.x[]*0.125/160.;
  boundary((scalar *) {muv});

  if(i%2==0){
    mgT = diffusion(TL, dt, D = 1.);
  }
  else{
    mgT = diffusion(TS, dt, D = 1.);
  }
}

/**
We produce an animation of the tracer field. */

event movies ( i +=2,last)
{
  boundary({TL,TS});
  view (fov = 16.642, quat = {0,0,0,1}, tx = -0.0665815, 
    ty = -0.00665815, bg = {1,1,1}, width = 600, height = 600, samples = 1);
  // draw_vof("cs");
  squares("TL", min =0, max = 1);
  save ("TL.gif");
  // squares("u.y");
  // save ("u.gif");
  // squares("u2.y");
  // save ("u2.gif");
  squares("dist", min = -0.5*NB_width , max = 0.5*NB_width);
  save ("dist.gif");
  squares("v_pc.y");
  save ("v_pcy.gif");
}

/**
We check the number of iterations of the Poisson and viscous
problems. */

event logfile (i++;i<=40){

  fprintf (stderr, "%d %g %g %d %d\n", i, t, dt, mgp.i, mgu.i);
}

