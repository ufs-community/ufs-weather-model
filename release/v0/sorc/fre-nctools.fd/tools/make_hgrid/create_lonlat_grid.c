#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include "mpp.h"
#include "mosaic_util.h"
#include "interp.h"
#include "tool_util.h"
#include "constant.h"
#include "create_hgrid.h"
#define  D2R (M_PI/180.)
#define  R2D (180./M_PI)

/*********************************************************************************
   some private routines used in this file
*********************************************************************************/
void set_regular_lonlat_grid( int nxp, int nyp, int isc, int iec, int jsc, int jec, double *xb, double *yb,
			      double *x, double *y, double *dx, double *dy, double *area, double *angle,
			      int use_great_circle_algorithm);
void set_f_plane_grid( int nxp, int nyp, int isc, int iec, int jsc, int jec, double *xb, double *yb, double f_plane_latitude, 
		       double *x, double *y, double *dx, double *dy, double *area, double *angle);
double cartesian_dist(double x1, double y1, double x2, double y2, double f_plane_latitude);
double cartesian_box_area(double x1, double y1, double x2, double y2, double f_plane_latitude);

/************************************************************************************
   void create_regular_lonlat_grid( int *nxbnds, int *nybnds, double *xbnds, double *ybnds,
                                    int *nlon, int *nlat, int *isc, int *iec,
                                    int *jsc, int *jec, double *x, double *y, double *dx,
                                    double *dy, double *area, double *angle_dx, int use_great_circle_algorithm )
   calculate grid location, length, area and rotation angle.
   The routine takes the following arguments

   INPUT:

   OUTPUT:

************************************************************************************/
void create_regular_lonlat_grid( int *nxbnds, int *nybnds, double *xbnds, double *ybnds,
		         	 int *nlon, int *nlat, double *dlon, double *dlat,
				 int use_legacy, int *isc, int *iec,
				 int *jsc, int *jec, double *x, double *y, double *dx,
				 double *dy, double *area, double *angle_dx, const char *center,
				 int use_great_circle_algorithm)
{
  int nx, ny, nxp, nyp, nxb, nyb;
  double *xb=NULL, *yb=NULL;
  int refine;
  double stretch = 1;
  
  /* use cubic-spline interpolation algorithm to calculate nominal zonal grid location. */
  nxb = *nxbnds;
  nyb = *nybnds;

  if( use_legacy ) {
    xb = compute_grid_bound_legacy(nxb, xbnds, dlon, stretch, &nx, center);
    yb = compute_grid_bound_legacy(nyb, ybnds, dlat, stretch, &ny, center);
  }
  else {
    xb = compute_grid_bound(nxb, xbnds, nlon, &nx, center);
    yb = compute_grid_bound(nyb, ybnds, nlat, &ny, center);
  }
  nxp = nx + 1;
  nyp = ny + 1;
  /*
    Test if left and right hand resolutions differ. Inform user of potential problem.
  */
  if ( fabs(xb[1] - xb[0] - xb[nx] + xb[nx-1])/(xb[1] - xb[0]) > 1.e-6 ) {
    printf("%d\n",nxp);
    printf("%s\n" ,"Note: End point resolutions differ for x-axis. Not suitable for periodic axes");
    printf("%s\n" ,"      See documentation for generating periodic axes when  center = 'c_cell'");
  }
  
  if ( fabs(yb[1] - yb[0] - yb[ny] + yb[ny-1])/(yb[1] - yb[0]) > 1.e-6 ) {
    printf("%d\n",nyp);
    printf("%s\n" ,"Note: End point resolutions differ for y-axis. Not suitable for periodic axes");
    printf("%s\n" ,"      See documentation for generating periodic axes when  center = 'c_cell'");
  }   
     
  set_regular_lonlat_grid( nxp, nyp, *isc, *iec, *jsc, *jec, xb, yb, x, y, dx, dy, area, angle_dx, use_great_circle_algorithm);
  free(xb);
  free(yb);
    
}; /* create_regular_lonlat_grid */


/************************************************************************************
   void create_simple_cartesian_grid( double *xbnds, double *ybnds, int *nlon, int *nlat,
                                      double *simple_dx, *simple_dy, int *isc, int *iec,
				      int *jsc, int *jec, double *x, double *y,
				      double *dx, double *dy, double *area, double *angle_dx )
   calculate grid location, length, area and rotation angle.
   The routine takes the following arguments

   INPUT:

   OUTPUT:

************************************************************************************/
void create_simple_cartesian_grid( double *xbnds, double *ybnds, int *nlon, int *nlat,
				   double *simple_dx, double *simple_dy, int *isc, int *iec,
				   int *jsc, int *jec, double *x, double *y,
				   double *dx, double *dy, double *area, double *angle_dx )
{
  int nx, ny, nxp, nyp, i, j, n, nxc, nyc, nxb, nyb;
  double *grid1=NULL, *grid2=NULL, *xb=NULL, *yb=NULL;

  nxb = 2;
  nyb = 2;
  nx = *nlon;
  ny = *nlat;
  nxp = nx + 1;
  nyp = ny + 1;  
  /* use cubic-spline interpolation algorithm to calculate nominal zonal grid location. */
  xb    = (double *)malloc(nxp*sizeof(double));
  grid1 = (double *)malloc(nxb*sizeof(double));
  grid1[0] = 1;
  grid1[1] = nxp;
  grid2 = (double *)malloc(nxp*sizeof(double));
  for(i=0;i<nxp;i++) grid2[i] = i + 1.0;
  cubic_spline( nxb, nxp, grid1, grid2, xbnds, xb, 1e30,1e30);
  free(grid1);
  free(grid2);

  /* use cubic-spline interpolation algorithm to calculate nominal meridinal grid location. */
  yb    = (double *)malloc(nyp*sizeof(double));
  grid1 = (double *)malloc(nyb*sizeof(double));
  grid1[0] = 1;
  grid1[1] = nyp;
  grid2 = (double *)malloc(nyp*sizeof(double));
  for(j=0;j<nyp;j++) grid2[j] = j + 1.0;
  cubic_spline( nyb, nyp, grid1, grid2, ybnds, yb, 1e30,1e30);
  free(grid1);
  free(grid2);

  n = 0;
  for(j=0; j<nyp; j++) {
    for(i=0; i<nxp; i++) {
      x[n] = xb[i];
      y[n++] = yb[j];
    }
  }

  nxc = *iec - *isc + 1;
  nyc = *jec - *jsc + 1;
  
  for(n = 0; n< nxc*(nyc+1);     n++) dx[n] = *simple_dx;
  for(n = 0; n< (nxc+1)*nyc;     n++) dy[n] = *simple_dy;
  for(n = 0; n< nxc*nyc;         n++) area[n] = (*simple_dx)*(*simple_dy);
  for(n = 0; n< (nxc+1)*(nyc+1); n++) angle_dx[n] = 0;
  
  free(xb);
  free(yb);
    
}; /* create_simple_cartesian_grid */


/*******************************************************************************
   void create_spectral_grid( nt *nlon, int *nlat, int *isc, int *iec,
                              int *jsc, int *jec, double *x, double *y, double *dx,
                              double *dy, double *area, double *angle_dx )
   generate spectral horizontal grid
*******************************************************************************/
void create_spectral_grid( int *nlon, int *nlat, int *isc, int *iec,
			   int *jsc, int *jec, double *x, double *y, double *dx,
			   double *dy, double *area, double *angle_dx, int use_great_circle_algorithm )
{
  const int itermax = 10;
  const double epsln = 1e-15;
  int ni, nj, nx, ny, nxp, nyp, i, j, converge, iter;
  double dlon, z, p1, p2, p3, z1, pp, a, b, c, d, sum_wts; 
  double *xb, *yb, *lon, *lonb, *lat, *latb;
  double *sin_hem, *wts_hem, *sin_lat, *wts_lat;
  
  nx = *nlon;
  ny = *nlat;
  nxp = nx + 1;
  nyp = ny + 1;
  ni = nx/2;
  nj = ny/2;

  xb = (double *) malloc(nxp*sizeof(double));
  yb = (double *) malloc(nyp*sizeof(double));
  lon     = (double *) malloc(ni    *sizeof(double));
  lonb    = (double *) malloc((ni+1)*sizeof(double));
  lat     = (double *) malloc(nj    *sizeof(double));
  latb    = (double *) malloc((nj+1)*sizeof(double));
  sin_hem = (double *) malloc(nj/2*sizeof(double));
  wts_hem = (double *) malloc(nj/2*sizeof(double));
  sin_lat = (double *) malloc(nj*sizeof(double));
  wts_lat = (double *) malloc(nj*sizeof(double));
  dlon = 360./ni;
  for(i=0;i< ni;i++) lon[i] = i*dlon;
  for(i=0;i<=ni;i++) lonb[i] = (i-0.5)*dlon;

  for(j=0;j<nj/2;j++) {
    converge = 0;
    z = cos(M_PI*(j +0.75)/(nj + 0.5));
    for(iter=1; iter<=itermax; i++){
      p1 = 1.0;
      p2 = 0.0;

      for(i=1;i<=nj;i++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0*i - 1.0)*z*p2 - (i - 1.0)*p3)/i;
      }

      pp = nj*(z*p1 - p2)/(z*z - 1.0E+00);
      z1 = z;
      z  = z1 - p1/pp;
      if(fabs(z - z1) < epsln ) {
        converge = 1;
	break;
      }
    }
    if(! converge) mpp_error("create_spectral_grid: abscissas failed to converge "
			     "in itermax iterations");
    sin_hem [j]     = z;
    wts_hem [j]     = 2.0/((1.0 - z*z)*pp*pp);
  }

  for(j=0;j<nj/2;j++) {
    sin_lat[j]      = - sin_hem[j];
    sin_lat[nj-1-j] =   sin_hem[j];
    wts_lat[j]      =   wts_hem[j];
    wts_lat[nj-1-j] =   wts_hem[j];
  }

  for(j=0;j<nj;j++){
    lat[j]    = asin(sin_lat[j])*R2D;
  }

  latb[0] = -90.;
  latb[nj] = 90.;
  for(j=1;j<nj;j++) {
    sum_wts = sum_wts + wts_lat[j-1];
    latb[j] = asin(sum_wts-1.)*R2D;
  }

  for(i=0;i<=ni;i++) xb[i*2]   = lonb[i];
  for(i=0;i<ni; i++) xb[i*2+1] = lon[i];
  for(j=0;j<=nj;j++) yb[j*2]   = latb[j];
  for(j=0;j<nj; j++) yb[j*2+1] = lat[j];

  set_regular_lonlat_grid( nxp, nyp, *isc, *iec, *jsc, *jec, xb, yb, x, y, dx, dy, area, angle_dx, use_great_circle_algorithm );
  free(xb);
  free(yb);
  free(lon);
  free(lonb);
  free(lat);
  free(latb);
  free(sin_hem);
  free(wts_hem);
  free(sin_lat);
  free(wts_lat); 
  
}; /* create_spectral_grid */

/*******************************************************************************
   void set_regular_lonlat_grid( int nxp, int nyp, int isc, int iec, int jsc, int jec,
                                 double *xb, double *yb, double *x, double *y,
                                 double *dx, double *dy, double *area, double *angle )
   set geographic grid location, calculate grid length, area and rotation angle
   x and y are on global domain, the other fields are on compute domain 
*******************************************************************************/
void set_regular_lonlat_grid( int nxp, int nyp, int isc, int iec, int jsc, int jec, double *xb, double *yb,
			      double *x, double *y, double *dx, double *dy, double *area, double *angle,
			      int use_great_circle_algorithm)
{
  int n, i, j;
  double lon[4], lat[4];
  
  n = 0;
  for(j=0; j<nyp; j++) {
    for(i=0; i<nxp; i++) {
      x[n]   = xb[i];
      y[n++] = yb[j];
    }
  }
  /* zonal length */
  n = 0;
  for(j=jsc; j<=jec+1; j++) {
    for(i=isc; i<=iec; i++ ) {
      dx[n++] = spherical_dist(x[j*nxp+i], y[j*nxp+i], x[j*nxp+i+1], y[j*nxp+i+1] );
    }
  }

  /* meridinal length */
  n = 0;
  for(j=jsc; j<=jec; j++) {
    for(i=isc; i<=iec+1; i++ ) {
      dy[n++] = spherical_dist(x[j*nxp+i], y[j*nxp+i], x[(j+1)*nxp+i], y[(j+1)*nxp+i] );
    }
  }

  /* cell area */
  if(use_great_circle_algorithm) {
    double *x_rad=NULL, *y_rad=NULL;
    int nx, ny;

    nx = nxp-1;
    ny = nyp-1;
    
    /* since make_hgrid is limited to be run on 1 processor, we could assume nx = iec-isc+1 and ny=jec-jsc+1 */
    x_rad  = (double *)malloc(nxp*nyp*sizeof(double));
    y_rad  = (double *)malloc(nxp*nyp*sizeof(double));
    for(i=0; i<nxp*nyp; i++){
      x_rad[i] = x[i]*D2R;
      y_rad[i] = y[i]*D2R;
    }
    get_grid_great_circle_area(&nx, &ny, x_rad, y_rad, area);
    free(x_rad);
    free(y_rad);
  }
  else {
    n = 0;
    for(j=jsc; j<=jec; j++) {
      for(i=isc; i<=iec; i++ ) {
	area[n++] = box_area(x[j*nxp+i]*D2R, y[j*nxp+i]*D2R, x[(j+1)*nxp+i+1]*D2R, y[(j+1)*nxp+i+1]*D2R );
      }
    }
  }
  
  /* rotation angle */
  n = 0;
  for(j=jsc; j<=jec+1; j++) {
    for(i=isc; i<=iec+1; i++ ) angle[n++] = 0;   
  }

};  /* set_regular_lonlat_grid */


/*******************************************************************************
   void create_tripolar_grid( int int nxbnds, int nybnds, double *xbnds, double *ybnds,
                              int *nlon, int *nlat, double *lat_join_in, int *isc, int *iec,
                              int *jsc, int *jec, double *x, double *y, double *dx, double *dy,
			      double *area, double *angle_dx )
   create tripolar grid infomrmation
********************************************************************************/
void create_tripolar_grid( int *nxbnds, int *nybnds, double *xbnds, double *ybnds,
			   int *nlon, int *nlat, double *dlon, double *dlat,
			   int use_legacy, double *lat_join_in, int *isc, int *iec,
			   int *jsc, int *jec, double *x, double *y, double *dx, double *dy,
			   double *area, double *angle_dx, const char *center, unsigned int verbose,
			   int use_great_circle_algorithm)
{
  int nxb, nyb, i, j, nx, ny, nxp, nyp, j_join, n, ip1, ii, n_count;
  double lat_join, lon_start, lon_end, lon_bpeq, lon_bpnp, lon_bpsp;
  double lam0, rp, lon_last, lon_scale;
  double *xb, *yb;
  double x_poly[20], y_poly[20];
  double stretch = 1;
  
  nxb = *nxbnds;
  nyb = *nybnds;

  if( use_legacy ) {
    xb = compute_grid_bound_legacy(nxb, xbnds, dlon, stretch, &nx, center);
    yb = compute_grid_bound_legacy(nyb, ybnds, dlat, stretch, &ny, center);
  }
  else {
    xb = compute_grid_bound(nxb, xbnds, nlon, &nx, center);
    yb = compute_grid_bound(nyb, ybnds, nlat, &ny, center);
  }
  
  nxp = nx + 1;  
  nyp = ny + 1;

  if ( fabs(xb[1] - xb[0] - xb[nx] + xb[nx-1])/(xb[1] - xb[0]) > 1.e-6 ) {
    printf("%d\n",nxp);
    printf("%s\n" ,"Note: End point resolutions differ for x-axis. Not suitable for periodic axes");
    printf("%s\n" ,"      See documentation for generating periodic axes when  center = 'c_cell'");
  }
  
  n = 0;
  for(j=0; j<nyp; j++) {
    for(i=0; i<nxp; i++) {
      x[n] = xb[i];
      y[n++] = yb[j];
    }
  }

  lat_join = *lat_join_in;
  j_join = nearest_index(lat_join, yb, nyp);
  lat_join  = yb[j_join];
  lon_start = xbnds[0];
  lon_end   = xbnds[1];  
  lon_bpeq = lon_start + 90. ;
  lon_bpnp = lon_start;
  lon_bpsp = lon_start+180.;
  
  if(*lat_join_in != lat_join) {
    if(mpp_pe() == mpp_root_pe() && verbose )printf("NOTE: Change join latitude from %f to %f\n",*lat_join_in, lat_join);
  }
  /**********************************************************************
     Transform from bipolar grid coordinates (bp_lon, bp_lat) to
     geographic coordinates (geolon_t, geolat_t) following R. Murray,
     "Explicit generation of orthogonal grids for ocean models",
     1996, J.Comp.Phys., v. 126, p. 251-273.  All equation
     numbers refer to the Murray paper.
  **********************************************************************/

  lam0 = fmod(lon_bpeq*D2R + 2*M_PI, 2*M_PI);
  rp = tan((0.5*M_PI-lat_join*D2R)/2.);        /* eqn. 2 */

  /*   calculate zonal length */
  n = 0;
  for(j=*jsc;j<=*jec+1;j++){
    for(i=*isc;i<=*iec;i++){
      if(j < j_join )  /* regular grid region */
        dx[n++] = spherical_dist(x[j*nxp+i], y[j*nxp+i], x[j*nxp+i+1], y[j*nxp+i+1]);
      else /*bipolar region */
	dx[n++] = bipolar_dist(x[j*nxp+i], y[j*nxp+i], x[j*nxp+i+1], y[j*nxp+i+1], lon_bpeq, lon_bpsp, lon_bpnp, rp );
    }
  }
  
  /*   calculate meridinal length */
  n = 0;
  for(j=*jsc;j<=*jec;j++){
    for(i=*isc;i<=*iec+1;i++){
      if(j < j_join )  /* regular grid region */
        dy[n++] = spherical_dist(x[j*nxp+i], y[j*nxp+i], x[(j+1)*nxp+i], y[(j+1)*nxp+i]);
      else /*bipolar region */
	dy[n++] = bipolar_dist(x[j*nxp+i], y[j*nxp+i], x[(j+1)*nxp+i], y[(j+1)*nxp+i], lon_bpeq, lon_bpsp, lon_bpnp, rp );
    }
  }
  
  /*define tripolar grid location */
  for(j=j_join;j<nyp;j++) {
    for(i=0;i<nxp;i++) {
      lon_last = x[j*nxp+max(i-1,0)];
      tp_trans(&(x[j*nxp+i]),&(y[j*nxp+i]),lon_last,lon_start, lam0, lon_bpeq, lon_bpsp, lon_bpnp, rp);
    }
  }

 /*calculte cell area */
  n = 0;
  if(use_great_circle_algorithm) {
    double *cart_x=NULL, *cart_y=NULL, *cart_z=NULL;
    double *x_rad=NULL, *y_rad=NULL;
    
    /* since make_hgrid is limited to be run on 1 processor, we could assume nx = iec-isc+1 and ny=jec-jsc+1 */
    cart_x = (double *)malloc(nxp*nyp*sizeof(double));
    cart_y = (double *)malloc(nxp*nyp*sizeof(double));
    cart_z = (double *)malloc(nxp*nyp*sizeof(double));
    x_rad  = (double *)malloc(nxp*nyp*sizeof(double));
    y_rad  = (double *)malloc(nxp*nyp*sizeof(double));
    for(i=0; i<nxp*nyp; i++){
      x_rad[i] = x[i]*D2R;
      y_rad[i] = y[i]*D2R;
    }
    latlon2xyz(nxp*nyp, x_rad, y_rad, cart_x, cart_y, cart_z);
    get_grid_great_circle_area(&nx, &ny, cart_x, cart_y, cart_z, area);
    free(x_rad);
    free(y_rad);
    free(cart_x);
    free(cart_y);
    free(cart_z);
  }
  else { 
    for(j=*jsc;j<=*jec;j++){
      for(i=*isc;i<=*iec;i++){
	if(j < j_join) 
	  area[n++] = box_area(x[j*nxp+i]*D2R, y[j*nxp+i]*D2R, x[(j+1)*nxp+i+1]*D2R, y[(j+1)*nxp+i+1]*D2R);
	else {
	  x_poly[0] = x[j*nxp+i]*D2R;       y_poly[0] = y[j*nxp+i]*D2R;
	  x_poly[1] = x[j*nxp+i+1]*D2R;     y_poly[1] = y[j*nxp+i+1]*D2R;
	  x_poly[2] = x[(j+1)*nxp+i+1]*D2R; y_poly[2] = y[(j+1)*nxp+i+1]*D2R;
	  x_poly[3] = x[(j+1)*nxp+i]*D2R;   y_poly[3] = y[(j+1)*nxp+i]*D2R;
	  n_count = fix_lon(x_poly, y_poly, 4, M_PI);
	  area[n++] = poly_area(x_poly, y_poly, n_count);
	}
      }
    }  
  }
  /*calculte rotation angle at cell vertex */
  n = 0;
  for(j=*jsc;j<=(*jec)+1;j++){
    for(i=*isc;i<=(*iec)+1;i++){
      if(j < j_join) 
	angle_dx[n++] = 0.0;
      else {
	lon_scale = cos(y[j*nxp+i]*D2R);
	if(i==0 || i == nx ) /* always at the boundary */
	  angle_dx[n++] = 0;
	else
	  angle_dx[n++] = atan2(y[j*nxp+i+1]-y[j*nxp+i-1],(x[j*nxp+i+1]-x[j*nxp+i-1])*lon_scale)*R2D;
      }
    }
  }
  
}; /* create_tripolar_grid */

void set_f_plane_grid( int nxp, int nyp, int isc, int iec, int jsc, int jec, double *xb, double *yb, double f_plane_latitude, 
		       double *x, double *y, double *dx, double *dy, double *area, double *angle)
{
  int n, i, j;
  
  n = 0;
  for(j=0; j<nyp; j++) {
    for(i=0; i<nxp; i++) {
      x[n]   = xb[i];
      y[n++] = yb[j];
    }
  }
  /* zonal length */
  
  n = 0;
  for(j=jsc; j<=jec+1; j++) {
    for(i=isc; i<=iec; i++ ) {
      dx[n++] = cartesian_dist(x[j*nxp+i], y[j*nxp+i], x[j*nxp+i+1], y[j*nxp+i+1], f_plane_latitude );
    }
  }

  /* meridinal length */
  n = 0;
  for(j=jsc; j<=jec; j++) {
    for(i=isc; i<=iec+1; i++ ) {
      dy[n++] = cartesian_dist(x[j*nxp+i], y[j*nxp+i], x[(j+1)*nxp+i], y[(j+1)*nxp+i], f_plane_latitude );
    }
  }
  
  /* cell area */
  n = 0;
  for(j=jsc; j<=jec; j++) {
    for(i=isc; i<=iec; i++ ) {
      area[n++] = cartesian_box_area(x[j*nxp+i], y[j*nxp+i], x[(j+1)*nxp+i+1], y[(j+1)*nxp+i+1], f_plane_latitude );
    }
  }

  /* rotation angle */
  n = 0;
  for(j=jsc; j<=jec+1; j++) {
    for(i=isc; i<=iec+1; i++ ) angle[n++] = 0;   
  }
  
}

  
/************************************************************************************/
void create_f_plane_grid( int *nxbnds, int *nybnds, double *xbnds, double *ybnds,
			  int *nlon, int *nlat, double *dlon, double *dlat,
			  int use_legacy, double f_plane_latitude, int *isc, int *iec,
			  int *jsc, int *jec, double *x, double *y, double *dx,
			  double *dy, double *area, double *angle_dx, const char *center )
{
  int nx, ny, nxp, nyp, nxb, nyb;
  double *xb=NULL, *yb=NULL;
  int refine;
  double stretch = 1;
  
  /* use cubic-spline interpolation algorithm to calculate nominal zonal grid location. */
  nxb = *nxbnds;
  nyb = *nybnds;

  if( use_legacy ) {
    xb = compute_grid_bound_legacy(nxb, xbnds, dlon, stretch, &nx, center);
    yb = compute_grid_bound_legacy(nyb, ybnds, dlat, stretch, &ny, center);
  }
  else {
    xb = compute_grid_bound(nxb, xbnds, nlon, &nx, center);
    yb = compute_grid_bound(nyb, ybnds, nlat, &ny, center);
  }
  nxp = nx + 1;
  nyp = ny + 1;
  /*
    Test if left and right hand resolutions differ. Inform user of potential problem.
  */
  if ( fabs(xb[1] - xb[0] - xb[nx] + xb[nx-1])/(xb[1] - xb[0]) > 1.e-6 ) {
    printf("%d\n",nxp);
    printf("%s\n" ,"Note: End point resolutions differ for x-axis. Not suitable for periodic axes");
    printf("%s\n" ,"      See documentation for generating periodic axes when  center = 'c_cell'");
  }
  
  if ( fabs(yb[1] - yb[0] - yb[ny] + yb[ny-1])/(yb[1] - yb[0]) > 1.e-6 ) {
    printf("%d\n",nyp);
    printf("%s\n" ,"Note: End point resolutions differ for y-axis. Not suitable for periodic axes");
    printf("%s\n" ,"      See documentation for generating periodic axes when  center = 'c_cell'");
  }   
     
  set_f_plane_grid( nxp, nyp, *isc, *iec, *jsc, *jec, xb, yb, f_plane_latitude, x, y, dx, dy, area, angle_dx);
  free(xb);
  free(yb);
    
}; /* create_f_plane_grid */


/* distance between cartesian grid on the earth */
double cartesian_dist(double x1, double y1, double x2, double y2, double f_plane_latitude)
{
  double dist = 0.0;

  if(x1 == x2)
    dist = fabs(y2-y1)*D2R*RADIUS;
  else if(y1 == y2) 
    dist = fabs(x2-x1)*D2R*RADIUS*cos(f_plane_latitude*D2R);
  else
    mpp_error("create_lonlat_grid: This is not rectangular grid");

  return dist;
}

 
/* rectangular grid box area for cartesian grid */
double cartesian_box_area(double x1, double y1, double x2, double y2, double f_plane_latitude)
{
  double area=0;
  double dx, dy;
  dx = x2-x1;
  dx = dx * D2R;

  dy = y2-y1;
  dy = dy * D2R;
  area = dx*dy*cos(f_plane_latitude*D2R)*RADIUS*RADIUS;

  return area;
}

