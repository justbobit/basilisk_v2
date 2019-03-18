
/**
# Solidification with a plane interface

We investigate the solidification of a block of water, between four plates, two
at -20°C and two at +20°C. It serves at minimum working example to show
how the functions defined in
[elementary_body.h](/sandbox/qmagdelaine/phase_change/elementary_body.h)
can be used to handle temperature fields and solidification.

![Solification of a block of water with the temperature
field](solidification_mwe_corner/video_solidification.mp4)

We define the geometrical, temporal and resolution parameters: */

#define H0 L0/5.
#define L 10. // size of the box

#define MIN_LEVEL 5
#define LEVEL     6
#define MAX_LEVEL 7
#define dH_refine (2.*L0/(1 << LEVEL))

#define F_ERR 1e-10

#define T_END   6.
#define DT_MAX  0.1
#define DELTA_T 0.1 // for videos and measurements
#define Pi 3.141592653589793

/**
We use an advection solver and the functions defined in
[elementary_body.h](/sandbox/qmagdelaine/phase_change/elementary_body.h): */
// #include "grid/cartesian.h"
//#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "alimare/elementary_body.h"
#include "tension.h"
// #include "view.h"
#define BG 0.7 // light gray for background
#define DG 0. // dark gray

/**
# Physical parameters

We non-dimensionalise the problem with the initial thickness of the film
$h_0$, the diffusion coefficient of temperature in the ice $D_S$ and the
command temperature $\theta_0$. A diffusive timescale appears: $\tau_{DS} =
\frac{h_0^2}{D_S}$.

Thus the diffusion equation for the temperature in the non-dimensional
space reads:
$$
\frac{\partial \theta}{\partial t} = \Delta \theta,
$$
appropriately completed with Dirichlet conditions at the surface of
the ice, at the top and bottom:
$$
\left\{
\begin{array}{ll}
\theta = 1 & \text{at the top}\\
\theta = -1 & \text{at the bottom}\\
\theta = 0 & \text{at the interface}
\end{array}
\right.
$$

The interface recedes at a (dimensional) velocity
$$
v_\text{e} = \frac{1}{L_H}\,\left(\frac{\lambda_S}{\rho_S}\, \nabla \theta\vert_S
+ \frac{\lambda_L}{\rho_L}\, \nabla \theta\vert_L\right)
\sim \frac{\lambda}{\rho\, L_H}\frac{\theta_0}{h_0} \equiv V
$$
a Peclet number Pe$= \frac{V\, h_0}{D_S} =
\frac{\lambda\, \theta_0}{L_H\, \rho_S\, D_S} = \frac{c_m^S\, \theta_0}{L_H}$
also enters in the problem description. This ratio is around 0.25 for
liquid water and 0.12 for ice.

We need time factor to set the Dirichlet condition, its role is specified in
[elementary_body.h](/sandbox/qmagdelaine/phase_change/elementary_body.h). */

#define T_eq                  0.
#define D_L                   0.12
#define D_S                   1.
#define TL_inf                -2.
#define TS_inf                -2.
#define peclet_L              0.253
#define peclet_S              0.124
#define dirichlet_time_factor 10.
#define SIGMA                 0.0007
#define VISC                  1.e-3
#define Ray_min               10.*L0
#define Precoeff              20.*(T_eq-TS_inf)/(SIGMA*Ray_min)

/**
We allocate several scalar fields to describe both the
interface and temperature fields. */

scalar f[], temperature_L[], temperature_S[], tr_eq[];
scalar * interfaces = {f}, * tracers = {temperature_L, temperature_S};
face vector v_pc[];

/**
We just specifie the dirichlet condition at the top, bottom, right and left: */

temperature_L[top]      = dirichlet(TL_inf );
temperature_L[bottom]   = dirichlet(TL_inf );
// temperature_L[left]     = periodic();
// temperature_L[right]    = periodic();
// temperature_L[bottom]= neumann(TL_inf );
// temperature_L[right] = neumann(TL_inf );

temperature_S[bottom]   = dirichlet(TS_inf );
temperature_S[top]      = dirichlet(TS_inf );
// temperature_S[left]     = periodic();
// temperature_S[right]    = periodic();
// temperature_S[left]  = neumann(TS_inf );
// temperature_S[right] = neumann(TS_inf );

/**
The main function of the program, where we set the domain geometry to
be ten times larger than the initial thickness of ice: */

int main()
{
  size (L);
  // origin (0., -L0/2.);
  origin (0., 0.);
  N  = 1 << LEVEL;
  init_grid (N);
  DT = DT_MAX;
  f.sigma = SIGMA;
  run();
}

/**
The initial position of the interface is defined with this function: */

// #define plane(x, y, H) (sqrt((x-L0/3.)*(x-L0/3.)/45.+(y)*(y)/15.) - H)
// #define plane(x, y, H) (sqrt((x-L0/2.)*(x-L0/2.)/350.+(y-L0/2.)*(y-L0/2.)/5.) - H )
// double hexagon (double x, double y, double h)
// {
//   double theta = atan2(x,-y);
//   double threshold1 = Pi/3.;
//   double threshold2 = 2.*Pi/3.;
//   return (theta < threshold1 ? -x+sqrt(3.*h)*(y+1.)  :
//         (theta < threshold2 ? -x+sqrt(3.*h)/2. : -x+sqrt(3.*h)*(1.-y)));
//   // return (-x+sqrt(3.)/2. );
//   // return ( x+sqrt(3.)*(y-1.) );
          
// }

double plane (double x, double y, double h)
{
  // double theta = atan2(x,-y);
  // double threshold1 = Pi/3.;
  // double threshold2 = 2.*Pi/3.;
  return (y+sin(6.*Pi*x/L0)-h);
  // return (-x+sqrt(3.)/2. );
  // return ( x+sqrt(3.)*(y-1.) );
          
}

/**
Before the first step, we initialize the temperature fields (after having
refined the grid around the future interface). */

event init (i = 0) {
  scalar curve[];
  #if TREE
    refine (level < MAX_LEVEL && plane(x, y, (H0 - dH_refine)) > 0.
            && plane(x, y, (H0 + dH_refine)) < 0.);
  // refine (level < MAX_LEVEL && plane(x, y, H0 - dH_refine ) < 0.
  //           && plane(x, y, H0 + dH_refine) > 0.);
    int iii = MAX_LEVEL;
    printf( " %i \n", iii);
  #endif
  fraction (f, plane(x, y, H0));

  boundary({f});
  curvature (f,curve);
  boundary({curve});
  foreach() {
    temperature_L[] = f[]*TL_inf;
    temperature_S[] = (1. - f[])*TS_inf;
    // tr_eq[] = T_eq;
     tr_eq[] = (f[] != 0. && f[] != 1. ? 
              -Precoeff*SIGMA*clamp(fabs(curve[]), 0., 1./Ray_min): 0.); 
  }
  foreach_face()
    uf.x[] = 0.;
  boundary({temperature_L, temperature_S, uf, tr_eq});
  CFL = 0.2;

  output_ppm (f, n=600, file = "init.png"\
  , min = 0., max = 1.); 
  output_ppm (tr_eq, n=600, file = "treq_init.png", min = -0.1, max = 0.); 
  stats s = statsf (tr_eq);
  printf ( " %f %f %f \n", s.min, s.max, s.stddev);

}

/**
## Solidification velocity

The velocity due to solidification is computed in the *stability()* event to
take into account this velocity in the CFL condition. */

event stability (i++) {

  /**
  The first step is used to make diffuse the temperature, because at the
  beginning, the concentration gradient diverges at the interface. We
  artificially make the simulation begin at a time $t_0 < 0$. In order to do
  this, we define a first timestep very short (we remain at $t \sim 0$), during
  which the interface does not move but the vapor diffuses over a finite time
  larger than the timestep. */

  if (i == 0)
    DT = 0.001; 
  else {
  
    /**
    In the evaporation examples
    ([static_drop.c](/sandbox/qmagdelaine/phase_change/1_elementary_body/static_drop.c) and
    [static_film.c](/sandbox/qmagdelaine/phase_change/1_elementary_body/static_film.c)),
    the evaporation velocity was due to the diffusive flux of the vapor in the
    gaz phase. Here, the solification is an equilibrium between the heat
    transfers of each side of the interface. */
  
    DT                    = DT_MAX;
    temperature_L.D       = D_L;
    temperature_L.inverse = false;
    temperature_L.peclet  = peclet_L;
    temperature_S.D       = D_S;
    temperature_S.inverse = true;
    temperature_S.peclet  = peclet_S;
    foreach_face()
      v_pc.x[] = 0.;

    for (scalar t in tracers) {    
      face vector tv[];
      phase_change_velocity (f, t, tv);
      foreach_face(){
        v_pc.x[] += (t.inverse ? tv.x[] : - tv.x[]);
        uf.x[]   += (t.inverse ? tv.x[] : - tv.x[]);
      }
    }
    boundary((scalar*){uf});
  }
}

/**
After the *vof()* event, the evaporation velocity has to be erased. */

event tracer_advection (i++) {
  foreach_face()
    uf.x[] = 0.;
  boundary((scalar*){uf});
}

/**
## Diffusion with immersed dirichlet condition

The temperature diffuses at each timestep. We need for that the maximal
level in the simulation. */

event tracer_diffusion(i++) {
  scalar curve[];
  #if TREE
    int max_level = MAX_LEVEL;
  #else
    int max_level = LEVEL;
  #endif

  temperature_L.D       = D_L;
  temperature_L.inverse = false;
  temperature_S.D       = D_S;
  temperature_S.inverse = true;
  // for (scalar t in tracers)  
    // t.tr_eq = T_eq;
  boundary({f});

  curvature (f,curve);

  boundary({curve});

  foreach(){
     // tr_eq[] = (f[] != 0. && f[] != 1. ? 
     //          -Precoeff*SIGMA*clamp(fabs(curve[]), 0., 1./Ray_min): 0.); 
    tr_eq[]  = T_eq;
    double nn = 0.;
    foreach_dimension()
      nn +=    (v_pc.x[]*v_pc.x[] + v_pc.x[1]*v_pc.x[1]);

#if dimension == 2
    tr_eq[] -= sqrt(nn)/4.*VISC;
#else 
    tr_eq[] -= sqrt(nn)/6.*VISC;
#endif
  }

  boundary({tr_eq}); 

  if (i == 0) {
    double dt_save = dt;
    dt = 0.05;
    for (scalar t in tracers)  
      dirichlet_diffusion (t, f, tr_eq, max_level, dt, dirichlet_time_factor);
    dt = dt_save;
  }
  else {
    for (scalar t in tracers)  
      dirichlet_diffusion (t, f, tr_eq, max_level, dt, dirichlet_time_factor);
  }
}

#if TREE
event adapt (i++) {
  adapt_wavelet ({f, temperature_L, temperature_S},\
                 (double[]){1e-3, 1e-3, 1e-3}, minlevel = MIN_LEVEL, \
                 maxlevel = MAX_LEVEL);
}
#endif

/**
## Video

We now juste have to save the video.
*/

event movie (t = 0.; t += max(DELTA_T, DT); t <= T_END)
{
  stats s = statsf (tr_eq);

  printf ( "%d %f %f %f %f \n", i, t, s.min, s.max, s.stddev);
  /**
  We mask out dead cells (i.e. cells for which `age` is zero). */

  // scalar m[];
  // scalar t_fonte[];

  // foreach()
    // m[] = age[] ? 1 : -1;

  // output_ppm (age, mask = m, n = 512, file = "age.gif", opt = "--delay 1",
  //   min = 0, max = 60);


  scalar temperature[];
  foreach() {
    f[] = clamp(f[], 0., 1.);
    temperature[] = (f[])*temperature_L[] + (1.-f[])*temperature_S[];
  }
  boundary ({f, temperature,tr_eq});

  output_ppm (f, n = 512, linear = true, file = "f.gif", opt = "--delay 1",\
              min = 0, max = 1);


  output_ppm (temperature, n = 512, linear = true, file = "T.gif", opt = "--delay 1",min = -2, max = 0);

  output_ppm (tr_eq, n = 512, linear = true, file = "T_solid.gif", \
   opt = "--delay 1", min = -0.4, max = 0);

  // output_ppm (tr_eq, n = 512, linear = true, file = "T_solid.gif", opt = "--delay 1",
  //   min = -0.2, max = 0);
}