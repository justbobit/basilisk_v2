/** We want to reproduce the geometry from this page:
https://en.wikipedia.org/wiki/Constructive_solid_geometry


 
#Utilisation of boolean operation for geometric construction

The goal is to define a geometry by using several boolean operation. This
geometry is describe
[here](https://en.wikipedia.org/wiki/Constructive_solid_geometry)
*/

#include "grid/octree.h"
#include "utils.h"
#include "fractions.h"
#include "view.h"

/**
##Definition of the boolean operation

Between 2 volumes (or surface in 2D), we can define a union, an intersection and a substraction of the volume.

Let use A and B, 2 volumes of the space, such that $A = f(x,y,z)$ and $B = g(x,y,z)$. Then:

$$A \cap B = max (f(x,y,z), g(x,y,z))$$
$$A \cup B = min (f(x,y,z), g(x,y,z))$$
$$A - B = max(f(x,y,z), -g(x,y,z))$$

Pay attention, in Basilisk, the surface expressions have to be continuous
function.

For example, a sphere is define with the following function:

$$f(x,y,z) = (x - x_c)^2 + (y-y_c)^2 + (z-z_c)^2 - R^2$$

In the sphere, $f$ is negative, allong the interface, $f$ is equal to zero and
outside of the sphere, $f$ is positive. The function is continuous. 

In basilisk, once you gave $f$ to the macros "fraction", the code will generate
an interface when $f$ is equal to zero.

##Definition of the geometric interface

We define first 2 basic geometric functions, a shpere and a cube.

The cube is centered on "center", and has a lenght of "size".
*/

double rectangle(double x, double y, coord center, coord size) {

  /**
  Our cube is defined as the intersection of 6 orthogonal planes. We
  define first 2 planes, $P1_{Plus}$ and $P1_{Minus}$.
  
  Then we define P1 has: 
  $$ P1 = P1_{Plus}\cap -P1_{Minus}$$*/

  double P1_Plus = x - size.x/2. + center.x;
  double P1_Minus = x + size.x/2. + center.x;

  double P1 = max (P1_Plus, -P1_Minus);

  /**
  We apply the same process to otain P2 and P3 */

  double P2_Plus = y - size.y/2. + center.y;
  double P2_Minus = y + size.y/2. + center.y;

  double P2 = max (P2_Plus, -P2_Minus);

  double c = max ( P1,P2 );

  return c;
}

/**
The sphere function will return a sphere, centered on "center", of
radius "radius"*/

double circle(double x, double y,  coord center, double radius) {
  return ( sq(x - center.x) + sq (y - center.y) - sq (radius));
}

/**
We will generate the CSG geometry examples. This geometry is define by
using a cube (C); a sphere (S) and 3 cylinders oriented allong $x$,
$y$ and $z$ (respectively $X$, $Y$ and $Z$).

$$ (C \cap S) - (X \cup Y \cup Z)$$*/

double geometry(double x, double y) {

  coord center_circle, center_rectangle, size_rectangle;
  center_circle.x = center_rectangle.x = 0.5;
  center_circle.y = 0.75;

  center_rectangle.y = 0.725;

  size_rectangle.x = 0.05;
  size_rectangle.y = 0.25; 

  /**
  We define the sphere and the cube */
  double s = circle (x, y, center_circle, 0.15);
  double r = rectangle (x, y, center_rectangle, size_rectangle);

  /**
  sIc is the intersection between the sphere and the cube. */
  double zalesak = max (s, -r);

  return zalesak;
}


/**
We have define all the important function for the geometry
generation. We can now implemet that into a basilisk mesh. */

int main() {

  /**
  We shift the origin so that the simulation domain will be centered
  on $(0,0,0)$ */

  /**
  We initialise the grid with 7 levels. */
  init_grid(1<<7);

  /**
  We create the scalar field that will get the geometry. */
  scalar f[];

  /**
  To have a refine mesh on the geometry, we initialise an iterative compt */
  int iteration = 0;
  do {
    iteration++;
    fraction(f, geometry (x, y));
  }while (adapt_wavelet({f}, (double []){0.2},
    maxlevel = 9, 2).nf != 0 && iteration <= 10);

  /**
  The geometry obtain with basilisk can be observe by using bview. The
  generated picture is:
  
  ![Reconstructed VOF surface.](csgBool/vof.png)
 */
  
  scalar l[];
    // foreach()
    //   l[] = level;
    // output_ppm (l, file = "levels.png", n = 400, min = 0, max = 7);

    foreach()
      l[] = f[];
    output_ppm (l, file = "vof.png", n = 400, min = 0, max = 1);
}
