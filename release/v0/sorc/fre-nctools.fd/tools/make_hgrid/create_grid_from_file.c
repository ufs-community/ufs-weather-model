/*
  Copyright 2011 NOAA Geophysical Fluid Dynamics Lab, Princeton, NJ
  This program is distributed under the terms of the GNU General Public
  License. See the file COPYING contained in this directory
*/
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include "mpp.h"
#include "mosaic_util.h"
#include "tool_util.h"
#include "constant.h"
#include "mpp_io.h"
#include "create_hgrid.h"
#define  D2R (M_PI/180.)
#define  R2D (180./M_PI)

void create_grid_from_text_file( char *file, int *nlon, int *nlat, double *x, double *y, double *dx, double *dy,
				 double *area, double *angle_dx );
void get_grid_v1(int fid, double *xt, double *yt, double *xc, double *yc);
void get_grid_v2(int fid, int nlon, int nlat, double *xt, double *yt, double *xc, double *yc);
void get_grid_v3(int fid, int nlon, int nlat, double *xt, double *yt, double *xc, double *yc);
void get_grid_v4(int fid, int nlon, int nlat, double *xt, double *yt, double *xc, double *yc);

/***********************************************************************
  void create_grid_from_file( char *file, int *nlon, int *nlat, double *x, double *y, double *dx,
                              double *dy, double *area, double *angle_dx )
   the grid location is defined through ascii file. calculate cell length,
   cell area and rotation angle  
************************************************************************/
void create_grid_from_file( char *file, int *nlon, int *nlat, double *x, double *y, double *dx, double *dy,
           		    double *area, double *angle_dx, int use_great_circle_algorithm )
{
  double *xt, *yt, *xc, *yc;
  double p1[4], p2[4], cart_x[4], cart_y[4], cart_z[4];
  int nx, ny, nxp, nyp, ni, nj, i, j, n;
  int is_uniform;
  char mesg[256], txt[128];
  
  /************************************************************
     identify the grid_file is ascii file or netcdf file,
     if the file name contains ".nc", it is a netcdf file,
     otherwise it is ascii file.
  *********************************************************/
  nx  = *nlon;
  ny  = *nlat;
  nxp = nx + 1;
  nyp = ny + 1;
  if(strstr(file, ".nc") ) {
    int fid, vid;
    
    ni  = *nlon/2;
    nj  = *nlat/2;
    xc = (double *)malloc((ni+1)*(nj+1)*sizeof(double));
    yc = (double *)malloc((ni+1)*(nj+1)*sizeof(double));
    xt = (double *)malloc( ni   * nj   *sizeof(double));
    yt = (double *)malloc( ni   * nj   *sizeof(double));

    fid = mpp_open(file, MPP_READ);
    if(mpp_dim_exist(fid, "grid_xt") ) 
      get_grid_v1(fid, xt, yt, xc, yc);
    else if(mpp_dim_exist(fid, "rlon") || mpp_dim_exist(fid, "i") )
      get_grid_v2(fid, ni, nj, xt, yt, xc, yc);
    else if(mpp_dim_exist(fid, "lon") )
      get_grid_v3(fid, ni, nj, xt, yt, xc, yc);
    else if(mpp_dim_exist(fid, "x") )
      get_grid_v4(fid, ni, nj, xt, yt, xc, yc);
    
    mpp_close(fid);
    for(j=0; j<nj+1; j++) for(i=0; i<ni+1; i++) {
      x[j*2*nxp+i*2] = xc[j*(ni+1)+i];
      y[j*2*nxp+i*2] = yc[j*(ni+1)+i];
    }
    for(j=0; j<nj; j++) for(i=0; i<ni; i++) {
      x[(j*2+1)*nxp+i*2+1] = xt[j*ni+i];
      y[(j*2+1)*nxp+i*2+1] = yt[j*ni+i];
    }
    /* check if the grid is uniform of not */
    is_uniform = 1;
    for(j=0; j<nj; j++) {
      for(i=1; i<ni; i++) {
	if(yt[j*ni+i] != yt[j*ni]) is_uniform = 0;
      }
    }
    
    for(i=0; i<ni; i++) {
      for(j=1; j<nj; j++) {
	if(xt[j*ni+i] != xt[i]) is_uniform = 0;
      }
    }

    if(is_uniform) {
       for(j=0; j<nj+1; j++) for(i=0; i<ni; i++) {
         x[j*2*nxp+i*2+1] = xt[i];   
         y[j*2*nxp+i*2+1] = yc[j*(ni+1)+i];
       }
       for(j=0; j<nj; j++) for(i=0; i<ni+1; i++) {
	 x[(j*2+1)*nxp+i*2] = xc[j*(ni+1)+i];
	 y[(j*2+1)*nxp+i*2] = yt[j*ni];
       }
    }
    else {
      for(j=0; j<nj+1; j++) for(i=0; i<ni; i++) {
	x[j*2*nxp+i*2+1] = (xc[j*(ni+1)+i]+xc[j*(ni+1)+i+1])*0.5;
	y[j*2*nxp+i*2+1] = (yc[j*(ni+1)+i]+yc[j*(ni+1)+i+1])*0.5;
      }    
      for(j=0; j<nj; j++) for(i=0; i<ni+1; i++) {
	x[(j*2+1)*nxp+i*2] = (xc[j*(ni+1)+i]+xc[(j+1)*(ni+1)+i])*0.5;
	y[(j*2+1)*nxp+i*2] = (yc[j*(ni+1)+i]+yc[(j+1)*(ni+1)+i])*0.5;
      }
    }
    for(j=0; j<nyp; j++) for(i=0; i<nx; i++) {
      p1[0] = x[j*nxp+i]*D2R; p2[0] = x[j*nxp+i+1]*D2R;
      p1[1] = y[j*nxp+i]*D2R; p2[1] = y[j*nxp+i+1]*D2R;
      dx[j*nx+i] = great_circle_distance(p1, p2);
    }
    for(j=0; j<ny; j++) for(i=0; i<nxp; i++) {
      p1[0] = x[j*nxp+i]*D2R; p2[0] = x[(j+1)*nxp+i]*D2R;
      p1[1] = y[j*nxp+i]*D2R; p2[1] = y[(j+1)*nxp+i]*D2R;
      dy[j*nxp+i] = great_circle_distance(p1, p2);
    }
    if(use_great_circle_algorithm) {
      for(j=0; j<ny; j++) for(i=0; i<nx; i++) {
        p1[0] = x[j*nxp+i]*D2R; p1[1] = x[(j+1)*nxp+i]*D2R; p1[2] = x[(j+1)*nxp+i+1]*D2R; p1[3] = x[j*nxp+i+1]*D2R;
	p2[0] = y[j*nxp+i]*D2R; p2[1] = y[(j+1)*nxp+i]*D2R; p2[2] = y[(j+1)*nxp+i+1]*D2R; p2[3] = y[j*nxp+i+1]*D2R;
        latlon2xyz(4, p1, p2, cart_x, cart_y, cart_z);
        area[n++] = great_circle_area(4, cart_x, cart_y, cart_z);
      }
    }
    else { 
      for(j=0; j<ny; j++) for(i=0; i<nx; i++) {
	p1[0] = x[j*nxp+i]*D2R; p1[1] = x[j*nxp+i+1]*D2R; p1[2] = x[(j+1)*nxp+i+1]*D2R; p1[3] = x[(j+1)*nxp+i]*D2R;
	p2[0] = y[j*nxp+i]*D2R; p2[1] = y[j*nxp+i+1]*D2R; p2[2] = y[(j+1)*nxp+i+1]*D2R; p2[3] = y[(j+1)*nxp+i]*D2R;
	area[j*nx+i] = poly_area(p1, p2, 4);
      }
    }
    /* currently set angle to 0 */
    for(j=0; j<nyp; j++) for(i=0; i<nxp; i++) angle_dx[j*nxp+i] = 0;    
    free(xt);
    free(yt);
    free(xc);
    free(yc);
  }
  else {
    create_grid_from_text_file(file, nlon, nlat, x, y, dx, dy, area, angle_dx );
  }
}; /* create_grid_from_file */



void create_grid_from_text_file( char *file, int *nlon, int *nlat, double *x, double *y, double *dx, double *dy,
           		    double *area, double *angle_dx )
{
  double *xb, *yb;
  int nx, ny, nxp, nyp, i, j, n;
  FILE *pFile;
  char mesg[256], txt[128];
  
  nx  = *nlon;
  ny  = *nlat;
  nxp = nx + 1;
  nyp = ny + 1;

  pFile = fopen (file,"r");
  if(!pFile) {
    strcpy(mesg, "create_grid_from_file: Can not open ascii file ");
    strcat(mesg,file);
    mpp_error(mesg);
  }

  fscanf(pFile, "%s%*[^\n]",txt); /* header line (ignored) */
  xb = (double *) malloc(nxp*sizeof(double));
  yb = (double *) malloc(nyp*sizeof(double));
  for(i=0;i<nxp;i++) 
    fscanf(pFile, "%lg", xb+i); /* longitude */ 
  fscanf(pFile, "%s%*[^\n]",txt); /* header line (ignored) */   
  for(j=0;j<nyp;j++) fscanf(pFile, "%lg", yb+j); /* latitude */
    
  fclose(pFile);
  n=0;
  for(j=0; j<nyp; j++) {
    for(i=0; i<nxp; i++) {
      x[n]       = xb[i];
      y[n]       = yb[j];
      angle_dx[n++] = 0;   
    }
  }
  /* zonal length */
  n = 0;
  for(j=0; j<nyp; j++) {
    for(i=0; i<nx; i++ ) {
      dx[n++] = spherical_dist(x[j*nxp+i], y[j*nxp+i], x[j*nxp+i+1], y[j*nxp+i+1] );
    }
  }

  /* meridinal length */
  n = 0;
  for(j=0; j<ny; j++) {
    for(i=0; i<nxp; i++ ) {
      dy[n++] = spherical_dist(x[j*nxp+i], y[j*nxp+i], x[(j+1)*nxp+i], y[(j+1)*nxp+i] );
    }
  }

  /* cell area */
  n = 0;
  for(j=0; j<ny; j++) {
    for(i=0; i<nx; i++ ) {
      area[n++] = box_area(x[j*nxp+i]*D2R, y[j*nxp+i]*D2R, x[(j+1)*nxp+i+1]*D2R, y[(j+1)*nxp+i+1]*D2R );
    }
  }

  free(xb);
  free(yb);

}; /* create_grid_from_text_file */

void get_grid_v1(int fid, double *xt, double *yt, double *xc, double *yc)
{
  int vid;
  
  vid = mpp_get_varid(fid, "grid_lon");
  mpp_get_var_value(fid, vid, xc);
  vid = mpp_get_varid(fid, "grid_lat");
  mpp_get_var_value(fid, vid, yc);
  vid = mpp_get_varid(fid, "grid_lont");
  mpp_get_var_value(fid, vid, xt);
  vid = mpp_get_varid(fid, "grid_latt");
  mpp_get_var_value(fid, vid, yt); 

}

void get_grid_v2(int fid, int nlon, int nlat, double *xt, double *yt, double *xc, double *yc)
{
  int vid, i, j;
  double *x_vert=NULL, *y_vert=NULL;
  
  vid = mpp_get_varid(fid, "lon");
  mpp_get_var_value(fid, vid, xt);
  vid = mpp_get_varid(fid, "lat");
  mpp_get_var_value(fid, vid, yt);

  x_vert = (double *)malloc(nlon*nlat*4*sizeof(double));
  y_vert = (double *)malloc(nlon*nlat*4*sizeof(double));
  vid = mpp_get_varid(fid, "lon_vertices");
  mpp_get_var_value(fid, vid, x_vert);
  vid = mpp_get_varid(fid, "lat_vertices");
  mpp_get_var_value(fid, vid, y_vert);  

  for(j=0; j<nlat; j++) for(i=0; i<nlon; i++) {
    xc[j*(nlon+1)+i] = x_vert[(j*nlon+i)*4];
    yc[j*(nlon+1)+i] = y_vert[(j*nlon+i)*4];
  }

  for(j=0; j<nlat; j++) {
    xc[j*(nlon+1)+nlon] = x_vert[(j*nlon+nlon-1)*4+1];
    yc[j*(nlon+1)+nlon] = y_vert[(j*nlon+nlon-1)*4+1];
  }
      
  for(i=0; i<nlon; i++) {
    xc[nlat*(nlon+1)+i] = x_vert[((nlat-1)*nlon+i)*4+3];
    yc[nlat*(nlon+1)+i] = y_vert[((nlat-1)*nlon+i)*4+3];
  }

  xc[nlat*(nlon+1)+nlon] = x_vert[(nlat*nlon-1)*4+2];
  yc[nlat*(nlon+1)+nlon] = y_vert[(nlat*nlon-1)*4+2];

  free(x_vert);
  free(y_vert);
}

void get_grid_v3(int fid, int nlon, int nlat, double *xt, double *yt, double *xc, double *yc)
{
  int vid, i, j;
  double *lon=NULL, *lat=NULL;
  double *lon_bnds=NULL, *lat_bnds=NULL;

  lon = (double *)malloc(nlon*sizeof(double));
  lat = (double *)malloc(nlat*sizeof(double));
  lon_bnds = (double *)malloc(2*nlon*sizeof(double));
  lat_bnds = (double *)malloc(2*nlat*sizeof(double));	 
  
  vid = mpp_get_varid(fid, "lon");
  mpp_get_var_value(fid, vid, lon);
  vid = mpp_get_varid(fid, "lat");
  mpp_get_var_value(fid, vid, lat);
  vid = mpp_get_varid(fid, "lon_bnds");
  mpp_get_var_value(fid, vid, lon_bnds);
  vid = mpp_get_varid(fid, "lat_bnds");
  mpp_get_var_value(fid, vid, lat_bnds);  

  for(j=0; j<nlat; j++) for(i=0; i<nlon; i++) {
    xt[j*nlon+i] = lon[i];
    yt[j*nlon+i] = lat[j];
  }
    
  for(j=0; j<nlat+1; j++) {
    for(i=0; i<nlon; i++) {
      xc[j*(nlon+1)+i] = lon_bnds[2*i];
    }
    xc[j*(nlon+1)+nlon] = lon_bnds[2*nlon-1];
  }

  for(i=0; i<nlon+1; i++) {
    for(j=0; j<nlat; j++) {
      yc[j*(nlon+1)+i] = lat_bnds[2*j];
    }
    yc[nlat*(nlon+1)+i] = lat_bnds[2*nlat-1];
  }

  free(lon);
  free(lat);
  free(lon_bnds);
  free(lat_bnds);
  
}


void get_grid_v4(int fid, int nlon, int nlat, double *xt, double *yt, double *xc, double *yc)
{
  int vid, i, j;
  double *x_vert=NULL, *y_vert=NULL;

#define PHI1 50
#define PHI2 50.000001
#define CENTRAL_LON -107
#define CENTRAL_LAT 50
  
  double phi1, phi2;
  double central_lon, central_lat;
  double n, F, rho0, rho;
  double *lon=NULL, *lat=NULL, *x1=NULL, *y1=NULL, *x2=NULL, *y2=NULL;
  double theta;
  
  phi1 = PHI1;
  phi2 = PHI2;
  central_lon = CENTRAL_LON;
  central_lat = CENTRAL_LAT;
  phi1 *= D2R;
  phi2 *= D2R;
  central_lon *=D2R;
  central_lat *= D2R;

  printf("******************************************************************************\n");
  printf("NOTE from create_grid_from_file: this routine is specially for NARR data sets \n");
  printf("     The tools assume the four corner of grid is \n");
  printf("     (1.000N, 145.500W), (46.635N, 148.639E), (46.352N, 2.566W), (0.897N, 68.318W) \n");
  printf("     The standard parallel is 50 and 50.000001 \n");
  printf("     The reference (or central) longitude is -107 and the reference latitude is 50 \n");
  printf("******************************************************************************\n");
  
  /* read the lon and lat */
  lon = (double *) malloc(nlon*nlat*sizeof(double));
  lat = (double *) malloc(nlon*nlat*sizeof(double));
  x1   = (double *) malloc(nlon*nlat*sizeof(double));
  y1   = (double *) malloc(nlon*nlat*sizeof(double));
  vid = mpp_get_varid(fid, "lon");
  mpp_get_var_value(fid, vid, lon);
  vid = mpp_get_varid(fid, "lat");
  mpp_get_var_value(fid, vid, lat);

  /* Always use negative value for longitude */
  for(i=0; i<nlon*nlat; i++) {
    if(lon[i] > 0) lon[i] -= 360;
    xt[i]  = lon[i];
    yt[i]  = lat[i];
    lon[i] *= D2R;
    lat[i] *= D2R;
  }

  
  n = log(cos(phi1)/cos(phi2))/log(tan(0.25*M_PI+0.5*phi2)/tan(0.25*M_PI+0.5*phi1));
  F = cos(phi1)*pow(tan(0.25*M_PI+0.5*phi1),n);
  rho0 = F*pow(1/tan(0.25*M_PI+0.5*central_lat),n);
  for(i=0; i<nlon*nlat; i++) {
    rho  = F*pow(1/tan(0.25*M_PI+0.5*lat[i]),n);
    x1[i] = rho*sin(n*(lon[i]-central_lon));
    y1[i] = rho0 - rho*cos(n*(lon[i]-central_lon));
  }

  /* check if it is uniform */
  for(i=0; i<nlon; i++) {
    for(j=1; j<nlat; j++) {
      if( fabs(x1[j*nlon+i] - x1[i]) > 1.e-6 ) {
	printf(" i = %d, j = %d, x1 = %g, x2 = %g \n", i, j, x1[j*nlon+i], x1[i]);
	mpp_error(" x is not uniform \n");
      }
    }
  }

  printf(" x is uniform \n \n \n");

  for(j=0; i<nlat; j++) {
    for(i=1; i<nlon; i++) {
      if( fabs(y1[j*nlon+i] - y1[j*nlon]) > 1.e-6 ) {
	printf(" i = %d, j = %d, x1 = %g, x2 = %g \n", i, j, y1[j*nlon+i], y1[j*nlon]);
	mpp_error(" y is not uniform \n");
      }
    }
  }  

  printf(" y is uniform \n \n \n");

  /* get the x and y value at corner */

  x2 = (double *)malloc((nlon+1)*(nlat+1)*sizeof(double));
  y2 = (double *)malloc((nlon+1)*(nlat+1)*sizeof(double));

  for(i=1; i<nlon; i++) x2[i] = 0.5*(x1[i-1]+x1[i]);
  x2[0] = 2*x1[0] - x2[1];
  x2[nlon] = 2*x1[nlon-1] - x2[nlon-1];
  for(j=1; j<=nlat; j++) for(i=0; i<=nlon; i++) x2[j*(nlon+1)+i] = x2[i];

  for(j=1; j<nlat; j++) y2[j*(nlon+1)] = 0.5*(y1[(j-1)*nlon]+y1[j*nlon]);
  y2[0] = 2*y1[0] - y2[nlon+1];
  y2[nlat*(nlon+1)] = 2*y1[(nlat-1)*nlon] - y2[(nlat-1)*(nlon+1)];
  for(i=1; i<=nlon; i++) for(j=0; j<=nlat; j++) y2[j*(nlon+1)+i] = y2[j*(nlon+1)];
  
  /* convert x2 and y2 back to lon, lat */
  for(i=0; i<(nlon+1)*(nlat+1); i++) {
    rho = sqrt(x2[i]*x2[i]+pow(rho0-y2[i],2));
    if( n<0 ) rho = -rho;
    theta = atan(x2[i]/(rho0-y2[i]));
    xc[i] = central_lon + theta/n;
    yc[i] =  2*atan( pow(F/rho, 1/n)) - 0.5*M_PI;
  }

  /* convert back to degree */
  for(i=0; i<(nlon+1)*(nlat+1); i++) {
    xc[i] *= R2D;
    yc[i] *= R2D;
  }

  free(x1);
  free(x2);
  free(y1);
  free(y2);
  free(lon);
  free(lat);

}
