#include <math.h>
#include <stdio.h>
#include "prototypes.h"

#define gamma 6.670e-11L
#define twopi 6.2831853L
#define si2mg 1.0e5L
#define G_TEMP ((gamma) * (si2mg) )
#define SMALL 1e-10L
#define G_TEMP_x_DENSITY(d) ((G_TEMP) * (d))	


/* double gbox(double *x0,double *y0,double *z0, double *x1,double *y1,double *z1,double *x2,double *y2,double *z2,double *rho) { */
double gbox(POINT *pt, PRISM *pr, PARAMETER *pa) { 
  /*
    Function gbox computes the vertical attraction of a 
    rectangular prism.  Sides of prism are parallel to x,y,z axes,
    and z axis is vertical down.  
    
    Input parameters:
    Observation point is (x0,y0,z0) (pt->easting, pt->northing, pt->elev).  The prism extends from x1 (west)
    to x2(east), from y1(south) to y2(north), and from z1(top) to z2(bottom) in the x, y, and z 
    directions, respectively.  Density of prism is rho.  All 
    distance parameters in units of m; rho in units of kg/(m**3). 
    
     pt->easting(x) and pt->northing(y) 
   *  are the geographic coordinates of the point where the calculation is made.
   *  These points can be on a grid or random 
   *  (mimicing an actual survey for example).  
   
   *  The following parameters are properties of a prism. A prism
   *  must have all of the following parameters defined:
   
   *  south_edge, north_edge  represent the length of the prism in meters 
   *  west_edge, east_edge  represent the width of the prism in meters
   *  surf_to_top and surf_to_bot represent the depth of the prism in meters
   *  density is the rock density of the prism
   *
   
    Output parameters:
    Vertical attraction of gravity, g, in mGal, summed over all prisms in the model. 
    
  */
  int i;
  int x,y,z;
  double sum;
  double  rijk, ijk;
  double arg1, arg2, arg3;
  double  xs[2],ys[2],zs[2],isign[2];
  double g = 0;
  /*(void) fprintf (stderr, "In gbox\n"); */
  
  for (i=0; i < pa->N_units; i++) {
  
    /* (void) fprintf (stderr, "%f %f %f %f %f %f %f\n", *x0, *y0, *z0, *x1, *y1, *z1, *rho);*/
    /* xs[0] = *x0 - *x1; */
    xs[0] = pt->easting - (pr+i)->west;
   
    /* ys[0] = *y0 - *y1; */
    ys[0] = pt->northing - (pr+i)->south;
  
    /* zs[0] = *z0 - *z1; */
    zs[0] = pt->elev - pa->depth_to_top;
  
    /* xs[1] = *x0 - *x2; */
    xs[1] = pt->easting - (pr+i)->east;
  
    /* ys[1] = *y0 - *y2; */
    ys[1] = pt->northing - (pr+i)->north;
  
    /* zs[1] = *z0 - *z2; */
    zs[1] = pt->elev - (pr+i)->depth_to_bottom;
  
    isign[0] = -1.0;
    isign[1] = 1.0;
    /*(void) fprintf (stderr, "%lf %lf %lf %lf %lf %lf %lf\n", xs[0], xs[1], ys[0], ys[1], zs[0], zs[1], *rho  );*/
  
  
    sum=0.0;
    for (x=0; x<2; x++) {
      for (y=0; y<2; y++) {
        for (z=0; z<2; z++) {
	       rijk = sqrt(xs[x]*xs[x] + ys[y]*ys[y] + zs[z]*zs[z]);
	       ijk = isign[x]*isign[y]*isign[z];
	       arg1 = atan2((xs[x]*ys[y]),(zs[z]*rijk));
	       if (arg1 < 0.0) arg1 = arg1 + twopi;
	       arg2 = rijk+ys[y];
	       arg3 = rijk+xs[x];
	       if(arg2 <= 0.0) arg2 = SMALL;
	       if (arg3 <= 0.0) arg3 = SMALL;
	       arg2 = log(arg2);
	       arg3 = log(arg3);
	       sum += ijk*(zs[z]*arg1-xs[x]*arg2-ys[y]*arg3);
        }
      }
    }
    /* g = *rho * sum * G_TEMP; */
    /*g += (pa->density * sum * G_TEMP); */
    g += sum;
  }
  g *=  G_TEMP_x_DENSITY(pa->density);
  return g;
}

