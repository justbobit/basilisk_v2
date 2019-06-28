/**
# Phase change of an elementary body

This file proposes functions to handle the phase change of an elementary body
limited by the diffusion of a diffusive tracer (its vapor for an evaporation or
the temperature for solidification). This file has to be included after
[advection.h](/src/advection.h) or [centered.h](/src/navier-stokes/centered.h)
and [tracers.h](/src/tracers.h).

If not defined by the user we fixe a default value for F_ERR, the accepted error
over f to avoid division by zero. */

#ifndef F_ERR
  #define F_ERR 1e-10
#endif

/**
We use the VOF description of the interface, the Bell-Collela-Glaz advection
solver and the fully implicit reaction-diffusion solver with the following
header files:

* [vof.h](/src/vof.h): to advect the VOF tracer,
* [diffusion.h](/src/diffusion.h): to make diffuse the diffusive tracer,
* [curvature.h](/src/curvature.h): to use *interfacial()* function;
* [my_functions.h](/sandbox/qmagdelaine/my_functions.h): to use some general
functions like compute the normal in every cell. */

// #include "vof.h"
#include "diffusion.h"
#include "curvature.h"
#include "alex_functions.h"
#include "embed-tree.h"
/**
We add three attribute to the scalar field: a diffusion coefficient $D$, a 
equatilibrium value $t_eq$ and a Peclet number (for the evaporation this is the
ratio between the saturation vapour concentration and the liquid denstity, for
the solification is the ratio between the heat capacity times the a temperature
scale and the latent heat). */

attribute {
  double D;
  double rho;
  double therm_conduct;
  double tr_eq;
  double peclet;
}

/**
The modelisation of the evaporation of a pure liquid is devided here in two
parts:

* computation of the advection velocity field for the interface according to
Fick's law,
* diffusion of the diffusive tracer with an immersed dirichlet condition.

One function for each of these steps are defined.

## Phase change velocity 

This function evaluates the local gradient of the diffusive tracer in order to
compute the phase change velocity.

The inputs of the function are:

* $f$: VOF tracer,
* $tr$: diffusive tracer field,
* $\mathbf{v}_{pc}$: the phase change velocity. */

void phase_change_velocity_LS_embed (scalar cs, face vector fs, scalar tr,
 scalar tr2, face vector v_pc, scalar dist, double latent_heat) {
  
 /**
  The phase change velocity $\mathbf{v}_{pc}$ is

  $$
  \mathbf{v}_{pc} = \mathrm{Pe}\, D\, \nabla tr
  $$
  
  here we use the embed_gradient_face_x defined in embed that gives a proper
  definition of the gradients with embedded boundaries. */
  
  face vector gtr[], gtr2[];
  foreach_face()
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
  With the concentration gradient and the normal vector we can now compute the
  evaporation velocity $\mathbf{v}_{pc}$, following the lines drawn in
  [meanflow.c](/sandbox/popinet/meanflow.c). We define it as the product between
  the density ratio and the diffusive flow. Note that $\mathbf{v}_{pc}$ is weighted
  by the face metric. */
  
  foreach_face() {

    v_pc.x[] = 0.;
    if (cs[] > 0. && cs[] < 1.){

      /**
      To find the vapor neighbor cells and weight the averaging between them, we
      compute the normal (normalized w.r.t. the norm-2) to the interface. We
      have to compute it before the main loop, and not locally, to apply *boundary()*
      to it and get consistent values in the ghost cells. */

      coord n = facet_normal (point, cs, fs);
      normalize(&n);
      /**
      We compute the evaporation velocity. */
      v_pc.x[] = (gtr.x[]*n.x + gtr2.x[]*n.x)*fm.x[]*1e-4;
    }
  }
  boundary((scalar *){v_pc});
}


