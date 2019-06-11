/**
# Time-reversed VOF and Level_set advection in a vortex

This classical test advects and stretches an initially circular
interface in a non-divergent vortical flow. The flow reverts in time
and the interface should come back to its original position. The
difference between the initial and final shapes is a measure of the
errors accumulated during advection.

We will need the advection solver combined with the VOF advection
scheme and the reinitialization function of the LS function. */

#include "utils.h"
#include "advection.h"
#include "vof.h"
#include "../level_set.h"
#include "basic_geom.h"

/**
The volume fraction is stored in scalar field `f` which is listed as
an *interface* for the VOF solver. The level set function is a *tracer* `dist`.


We do not advect any *level set* with
the default (diffusive) advection scheme of the advection solver.  
 */

scalar f[], dist[];
scalar * interfaces = {f}, * tracers = {dist};
scalar * level_set = NULL;
/**
Here are the parameters for the simulation. We use a narrow band (NB) approach 
meaning that the level set function has meaning only in the direct vicinity of 
the 0 value of the level set function. For this test case, the NB is made of 
only 4 cells. 
 */

int     MAXLEVEL;
int     nb_cell_NB =  1 << 2 ;  // number of cells for the NB
double  NB_width ;              // length of the NB

int     N_display = 6;          // maximum level of refinement for which we 
                                // will display the results
double mass_ls_init, mass_vof_init;

/**
We center the unit box on the origin and set a maximum timestep of 0.1 */

int main() {
  origin (-L0/2.,-L0/2.);
  
  /**
  We then run the simulation for different levels of refinement. */

  periodic(right);
  for (MAXLEVEL = 6; MAXLEVEL <= N_display; MAXLEVEL++) {
    NB_width = L0*nb_cell_NB / (1<<MAXLEVEL);
    init_grid (1 << MAXLEVEL);
    run();
  }
}

/**
The initial interface is a circle of radius 0.2 centered on
(-0.2,-0.236338) (for historical reasons). We use the levelset
function `circle()` to define this interface. 

The period of the stretching cycle is set to 15, which will lead to
strong stretching. Milder conditions can be obtained by decreasing it. */

// #define circle2(x,y) (sq(0.2) - (sq(x + 0.2) + sq(y + .236338)))
#define T 3000.


coord center_circle ={0.,0.};
double Radius   =  0.25;

/**
We define the auxiliary levelset function $\phi$ on each vertex of the grid and
compute the corresponding volume fraction field. 

The level set function `dist` is taken positive out of the circle and we 
clamp the distance due to our NB approach. We take a 2\% overshoot that prevents
 NB cells from appearing due to spurious oscillations.
*/

event init (i = 0) {
  fraction (f, circle(x,y,center_circle,Radius));


  foreach_vertex(){
    dist[] = -circle(x,y,center_circle,Radius);
  }

  if (N == (1<<N_display)) {
    scalar l[];
    // foreach()
    //   l[] = level;
    // output_ppm (l, file = "levels_init.png", n = 400, min = 0, max = 7);

    foreach()
      l[] = f[];
    output_ppm (l, file = "f_init.png", n = 400, min = 0, max = 1);

    foreach()
      l[] = dist[];
    output_ppm (l, file = "dist_init.png", n = 400);
  }

  foreach_face(x){
    u.x[] = 0.001;
  }

}

/**
At the start and end of the simulation we check the sum, min and max
values of the volume fraction field. The sum must be constant to
within machine precision and the volume fraction should be bounded by
zero and one. */

event logfile (t = {0,T}) {
  stats s = statsf (f);
  fprintf (stderr, "# %f %.12f %.9f %g\n", t, s.sum, s.min, s.max);
  if(N == (1<<N_display)) {
    mass_vof_init = s.sum;
  }
  scalar f2[];
  fractions (dist, f2);
  s = statsf (f2);
  fprintf (stderr, "# %f %.12f %.9f %g\n", t, s.sum, s.min, s.max);
  if(N == (1<<N_display)) {
    mass_ls_init = s.sum;
  }
}


/**
To compute the errors, we reinitialise field `e` and `e2` at the end of the
simulation with the initial shape and compute the differences with the
final shape. We output the norms as functions of the maximum
resolution `N`.  

Note that the error of the level set function is only studied in half of the NB width.
*/

event field (t = T) {
  scalar e[], e2[];
  fraction (e, circle(x,y,center_circle,Radius));
  foreach_vertex()
    e2[] = -circle(x,y,center_circle,Radius);
  foreach(){
    e[]  -= f[];
    e2[] -= fabs(e2[])< 0.5*NB_width ? dist[] : e2[];
  }
  norm n  = normf (e);
  norm n2 = normf (e2);
  fprintf (stderr, "%d %g %g %g %g %g %g\n", N, n.avg, n.rms, n.max, 
            n2.avg, n2.rms, n2.max);
  scalar l[];
    foreach()
      l[] = dist[];
    output_ppm (l, file = "dist_final.png", n = 400,min = -0.1*NB_width,
      max = 0.1*NB_width);
}

/**
We also output the shape of the reconstructed interface at regular
intervals (but only on the finest grid considered). */

event shape (t += T/4.) {
  if (N == (1<<N_display))
    output_facets (f);
}

/**
If we are using adaptivity, we also output the levels of refinement at
maximum stretching. */

#if TREE
event levels (t = T/2) {
  if (N == (1<<N_display)) {
    scalar l[];
    foreach()
      l[] = level;
    output_ppm (l, file = "levels.png", n = 400, min = 0, max = 7);
  }
}
#endif


event levels (t = T/2) {
  if (N == N_display) {
    scalar l[];
    // foreach()
    //   l[] = level;
    // output_ppm (l, file = "levels.png", n = 400, min = 0, max = 7);

    foreach()
      l[] = f[];
    output_ppm (l, file = "f_reversed.png", n = 400, min = 0, max = 1);

    foreach()
      l[] = dist[];
    output_ppm (l, file = "dist.png", n = 400, min = -0.05*NB_width,
     max=0.05*NB_width);
  }
}
/**
We make simple movies showing the field `f`, the cells used for the NB and the
level set function during the simulation.
*/

event levels2 (t += T/300){
  if (N == (1<<N_display)) {
  
    scalar l[];

    foreach()
      l[] = f[];
    output_ppm (l, file = "f_reversed2.gif", n = 400, 
      opt = "--delay 1",min = 0, max = 1);

    // foreach()
    //   l[] = fabs(dist[]) < NB_width ? 1. : 0.;
    // output_ppm (l, file = "NB.gif", n = 400, 
    //   opt = "--delay 1",min = 0, max = 1);

    foreach()
      l[] = dist[];
    output_ppm (l, file = "dist2.gif", n = 400, 
      opt = "--delay 1",min = -1.*NB_width, max=1.*NB_width);
  }
}

/**
Level set reinitialization event. The number of iteration of the LS_reinit2 
function is set to  $1.4 \times 2 \times nb_{cell NB}$ which is a bit more than 
the $NB_{width}$
*/

event LS_reinitialization(i++,last){

  // WIP
  // scalar norm_v[];
  // foreach()
  //   norm_v[] =  0.;

  // foreach_face()
  //   norm_v[] += sq(uf.x[]);

  // stats s         = statsf(norm_v);
  // d_max_advected += sqrt(s.max)*dt;
  // if(d_max_advected > NB_width){
  //    d_max_advected = 0.;
  // END OF WIP

  if(i>0){
    LS_reinit2(dist,L0/(1 << MAXLEVEL), NB_width,
      1.4*(nb_cell_NB));
  }
}


#if 0
event movie (i += 10)
{
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, 1 << MAXLEVEL, file = "level.mp4");
}
#endif

/**
## Results

We use gnuplot (see [reversed.plot]()) to compute the convergence rate
of the error norms with and without adaptation. The convergence rates
are comparable.

~~~gnuplot Convergence rates for constant- and adaptive grids

ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)

f(x)=a+b*x
fit f(x) 'log' u (log($1)):(log($4)) via a,b
f2(x)=a2+b2*x
fit f2(x) 'log' u (log($1)):(log($2)) via a2,b2

f3(x)=a3+b3*x
fit f3(x) 'log' u (log($1)):(log($7)) via a3,b3
f4(x)=a4+b4*x
fit f4(x) 'log' u (log($1)):(log($5)) via a4,b4

set xlabel 'Maximum resolution'
set ylabel 'Maximum error'
set key bottom left
set logscale
set xrange [16:256]
set xtics 16,2,256
set grid ytics
set cbrange [1:1]

plot 'log' u 1:4 t 'VOF max (cart)'  , exp(f(log(x)))       t ftitle(a ,b ), \
     'log' u 1:2 t 'VOF norm1 (cart)', exp(f2(log(x))) lw 3 t ftitle(a2,b2), \
     'log' u 1:7 t 'LS max (cart)'   , exp(f3(log(x)))      t ftitle(a3,b3), \
     'log' u 1:5 t 'LS norm1 (cart)' , exp(f4(log(x))) lw 3 t ftitle(a4,b4)

~~~

The shapes of the interface at $t=0$, $t=T/4$, $t=T/2$, $t=3T/4$ and
$t=T$ are displayed below for both sets of simulations (constant and
adaptive), for $N=N_{display}$. The shapes for $t=T/4$ should be identical to
those for $t=3T/4$ and similarly for $t=0$ and $t=T$ (for which we
measure the error). Note that the errors for $t=3T/4$ seem to be much
larger than those for $t=T$.

~~~gnuplot Shapes of the interface
reset
set size ratio -1
plot [-0.5:0.5][-0.5:0.5]'out' w l t "adaptive"
~~~


~~~gnuplot Mass loss during calculation
reset
set yrange [-0.6:0.6]
set ylabel 'Mass error (%)'
set xlabel 'Time'
plot 'mass' u 1:2 t 'LS' w l lw 3, \
     'mass' u 1:3 t 'VOF' w l lw 3


~~~

![Animation of the level set function](reversed/dist2.gif)(width="800" height="600") 


Initial Level set field    |  Final Level Set field
:-------------------------:|:-------------------------:
![Initial](reversed/dist_init.png)(width="400" height="300") |  ![Final](reversed/dist_final.png)(width="400" height="300")  
*/
