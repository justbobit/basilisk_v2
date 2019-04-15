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

extern scalar * level_set;
extern face vector uf;
extern double dt;


/**
The integration is performed using the Bell-Collela-Glaz scheme. */

#include "bcg.h"

event LS_advection (i++,last) {
  // printf("%d \n", i);
  advection (level_set, uf, dt);
}
/**
Reinitialization method can be chosen by overloading this function */

event LS_reinitialization(i++,last) {
}
