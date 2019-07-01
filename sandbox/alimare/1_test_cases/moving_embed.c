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
#include "../elementary_body.h"
// #include "navier-stokes/perfs.h"
#include "tracer.h"
#include "basic_geom.h"
#include "view.h"

#define MIN_LEVEL 3
#define LEVEL     6
#define MAXLEVEL 8

#define H0 -L0/4.

#define T_eq         0.
#define D_L          0.12
#define D_S          1.
#define TL_inf       1.
#define TS_inf       -1.
#define peclet_L     0.253
#define peclet_S     0.124
#define rho_S        917.
#define rho_L        1000.
#define latent_heat  334
#define lambda_L     2.22e3
#define lambda_S     555.e3

#define plane(x, y, H) (y - H)


scalar TL[], TS[],dist[];
scalar * tracers = {TL, TS};
scalar * level_set = {dist};
face vector v_pc[];

face vector muv[];
mgstats mgT;

int     nb_cell_NB =  1 << 3 ;  // number of cells for the NB
double  NB_width ;              // length of the NB

TL.D             = D_L;
TL.peclet        = peclet_L;
TL.therm_conduct = lambda_L;
TL.rho           = rho_L;

TS.D             = D_S;
TS.peclet        = peclet_S;
TS.therm_conduct = lambda_S;
TS.rho           = rho_S;
  

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
}

/**
The fluid is injected on the left boundary with a unit velocity. The
tracer is injected in the lower-half of the left boundary. An outflow
condition is used on the right boundary. */

u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
The top and bottom walls are free-slip and the cylinder is no-slip. */

/** on the embed boundary : .n velocity is the velocity of the LS .t is 0.
*/

// u.n[embed] = fabs(y) > 1.5 ? neumann(0.) : 
//     dirichlet_embed_LS(point , cs , fs, v_pc);
u.n[embed] = fabs(y) == L0/2. ? neumann(0.) :\
  dirichlet_embed_LS(point , cs , fs, v_pc);
u.t[embed] = fabs(y) == L0/2. ? neumann(0.) : dirichlet(0.);

TL[embed]  = dirichlet(T_eq);
TL[top]    = dirichlet(TL_inf); 

TS[embed]    = dirichlet(T_eq);
TS[bottom]   = dirichlet(TS_inf); 

event init (t = 0)
{
  // DT = 5.e-4;
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

  /**
  We set the initial velocity field and set a rotating field inside the notched
  disk. */
  foreach_face()
    v_pc.x[] = 0.;

  boundary((scalar*){v_pc});
  
  foreach() {
    TL[] = cs[]*TL_inf;
    TS[] = (1. - cs[])*TS_inf;
  }

  boundary({TL, TS});
  // view (fov = 16.642, quat = {0,0,0,1}, tx = -0.0665815, 
    // ty = -0.00665815, bg = {1,1,1}, width = 600, height = 600, samples = 1);
  // draw_vof("cs");
  // squares("cs", min=0, max = 1);
  // save ("cs.png");
}

/**
We check the number of iterations of the Poisson and viscous
problems. */

event logfile (i++;t<=5.0){

  norm n  = normf (v_pc.x);

  fprintf (stderr, "%d %g %g %d %d %g %g\n", i, t, dt, mgp.i, mgu.i, n.max, n.avg);

}

event tracer_diffusion(i++){
  if(i%2 ==0){
    mgT =diffusion(TL, dt, D = TL.D);
    boundary({TL});
  }
  else{
    mgT = diffusion(TS, dt, D = TS.D);
    boundary({TS});
  }
  
  fprintf(ferr,"%t=%g\ti=%d\tDT=%g\n",t,i,dt);
}/**
Advection of the Level_set function. Based on the calculation of the velocity of
the phase_change 'v_pc' and of the gradients of both the tracers on the
interface using 'face_gradient_x'.
*/
event LS_advection(i++,last){
  if(i%2 == 1){
//     vertex scalar phi[];
//     foreach_vertex() {
//       phi[] = -geometry(x,y,0.13+0.025*t);
//     }
//     boundary ({phi});
//     fractions (phi, cs, fs);
//     /**
//     We set the velocity field in the covered cells. */
    
//     foreach()
//       u.x[] = cs[] ? cs[] : 0.; 
//     boundary ((scalar *){u});
//     trash ({uf});
//     foreach_face()
//       uf.x[] = fm.x[]*face_value (u.x, 0);
//     boundary ((scalar *){uf});

//     /**
//     We update fluid properties. */

//     event ("properties");

//     dtmax = DT;
//     event ("stability");

  /**
  In the evaporation examples
  ([static_drop.c](/sandbox/qmagdelaine/phase_change/1_elementary_body/static_drop.c) and
  [static_film.c](/sandbox/qmagdelaine/phase_change/1_elementary_body/static_film.c)),
  the evaporation velocity was due to the diffusive flux of the vapor in the
  gaz phase. Here, the solification is an equilibrium between the heat
  transfers of each side of the interface. */


  double L_H       = latent_heat;
  boundary({TL});
  restriction({TL});
  foreach_face(){
    v_pc.x[] = 0.;                        
    fs.x[]   = 1.-fs.x[];
  }
  foreach()
    cs[] = 1.-cs[];

  
  boundary({cs,fs});
  restriction({cs,fs});

  boundary({TS});
  restriction({TS});
  
  /**
  cs and fs have just been inverted, we create temporary fields just for
  calculation purposes

  When this calculation is done (1-cs) => TL, cs => TS (same for fs,(1-fs) )
  */
  phase_change_velocity_LS_embed (cs, fs ,TL, TS, v_pc, dist, L_H, NB_width);
  
  foreach_face(){
    fs.x[]  = 1.-fs.x[];
  }
  foreach()
    cs[] = 1.-cs[];

  boundary({cs,fs});
  restriction({cs,fs});

  // foreach_face()
  //   if(fabs(dist[])<NB_width)
  //     v_pc.x[] += fs.x[]*uf.x[] + (1.-fs.x[])*uf2.x[];


  advection (level_set, v_pc, dt);
  boundary ({dist});
  restriction({dist});

  fractions (dist, cs, fs);
  boundary({cs,fs});
  restriction({cs,fs});

  boundary ((scalar *){u});
  restriction ((scalar *){u});
  boundary ({pf,p,g});
  restriction ({pf,p,g});
  trash ({uf});
  foreach_face()
    uf.x[] = fm.x[]*face_value (u.x, 0);
  boundary ((scalar *){uf});
  restriction ((scalar *){uf});


  event ("properties");
  }
}

event LS_reinitialization(i++,last){
  if(i%2 == 1){
    LS_reinit2(dist,L0/(1 << MAXLEVEL), NB_width,
      1.4*(nb_cell_NB));
  }
}

/**
We produce an animation of the tracer field. */

event movies ( i +=10,last)
{
  scalar omega[];
  vector m[];
  foreach(){
    foreach_dimension()
      m.x[]  = cs[]*u.x[] + (1-cs[])*u2.x[];
  }
  boundary ((scalar *) {m});
  vorticity (m, omega);
  boundary ({omega});
  boundary({cs});
  boundary({TL});
  view (fov = 16.642, quat = {0,0,0,1}, tx = -0.0665815, 
    ty = -0.00665815, bg = {1,1,1}, width = 600, height = 600, samples = 1);
  draw_vof("cs");
  squares("m.y", min=-10, max = 10);
  save ("movie.mp4");
}

/**
We adapt according to the error on the embedded geometry, velocity and
tracer fields. */
#if DOUBLE_EMBED
// event adapt (i+=2) {
//   adapt_wavelet ({cs,u,u2}, (double[]){1e-2,3e-2,3e-2}, 8, 4);
// }
#else
event adapt (i++) {
  adapt_wavelet ({cs,u}, (double[]){1e-2,3e-2}, 8, 4);
}
#endif

/**
# Future work

* Synchronizing time steps in both phases
* Movement of the boundary must follow the movement of a physical quantity.
Use the velocity field near the boundary to advect a level set function and use
the 0-level set as the boundary for the embedded boundary.
*/
