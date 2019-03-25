/**
# Level-Set advection event

This event integrates advection equations of the form
$$
\partial_tLS+\mathbf{u_f}\cdot\nabla LS=0
$$
where $\mathbf{u_f}$ is the velocity field and $LS$ is the level set function.

The level set function is defined and initialized elsewhere (typically by the 
user), the face vector field `uf` and the timestep `dt` are defined by a
solver. */

extern scalar * LS;
extern face vector uf;
extern double dt;

/**
On adaptive meshes, the Level Set needs to use linear interpolation (rather
than the default bilinear interpolation) to ensure conservation when
refining cells. */

#if TREE
event defaults (i = 0) {
#if EMBED
    LS.refine = LS.prolongation = refine_embed_linear;
#else
    LS.refine  = refine_linear;
#endif
    LS.restriction = restriction_volume_average;
}
#endif

/**
The integration is performed using the Bell-Collela-Glaz scheme. */

#include "bcg.h"

event LS_advection (i++,last) {
  advection (LS, uf, dt);
}
/**
Reinitialization method can be chosen by overloading this function */

// event LS_reinitialization(i++,last) {


	
// }
