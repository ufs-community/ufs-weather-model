#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mosaic_util.h"
#include "create_xgrid.h"
#include "constant.h"

#define AREA_RATIO_THRESH (1.e-6)  
#define MASK_THRESH       (0.5)
#define EPSLN8            (1.e-8)
#define EPSLN30           (1.0e-30)
#define EPSLN10           (1.0e-10)
#define R2D (180/M_PI)
#define TPI (2.0*M_PI)
double grid_box_radius(const double *x, const double *y, const double *z, int n);
double dist_between_boxes(const double *x1, const double *y1, const double *z1, int n1,
			  const double *x2, const double *y2, const double *z2, int n2);
int inside_edge(double x0, double y0, double x1, double y1, double x, double y);
int line_intersect_2D_3D(double *a1, double *a2, double *q1, double *q2, double *q3,
		         double *intersect, double *u_a, double *u_q, int *inbound);


/*******************************************************************************
  int get_maxxgrid
  return constants MAXXGRID.
*******************************************************************************/
int get_maxxgrid(void)  
{
  return MAXXGRID;
}

int get_maxxgrid_(void)
{
  return get_maxxgrid();
}

/*******************************************************************************
void get_grid_area(const int *nlon, const int *nlat, const double *lon, const double *lat, const double *area)
  return the grid area.
*******************************************************************************/
#ifndef __AIX
void get_grid_area_(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area)
{
  get_grid_area(nlon, nlat, lon, lat, area);
}
#endif

void get_grid_area(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area)
{
  int nx, ny, nxp, i, j, n_in;
  double x_in[20], y_in[20];
  
  nx = *nlon;
  ny = *nlat;
  nxp = nx + 1;

  for(j=0; j<ny; j++) for(i=0; i < nx; i++) {
    x_in[0] = lon[j*nxp+i];
    x_in[1] = lon[j*nxp+i+1];
    x_in[2] = lon[(j+1)*nxp+i+1];
    x_in[3] = lon[(j+1)*nxp+i];
    y_in[0] = lat[j*nxp+i];
    y_in[1] = lat[j*nxp+i+1];
    y_in[2] = lat[(j+1)*nxp+i+1];
    y_in[3] = lat[(j+1)*nxp+i];
    n_in = fix_lon(x_in, y_in, 4, M_PI);    
    area[j*nx+i] = poly_area(x_in, y_in, n_in);
  }

};  /* get_grid_area */

#ifndef __AIX
void get_grid_great_circle_area_(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area)
{
  get_grid_great_circle_area(nlon, nlat, lon, lat, area);

}
#endif

void get_grid_great_circle_area(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area)
{
  int nx, ny, nxp, nyp, i, j, n_in;
  int n0, n1, n2, n3;
  double x_in[20], y_in[20], z_in[20];
  struct Node *grid=NULL;
  double *x=NULL, *y=NULL, *z=NULL;


  nx = *nlon;
  ny = *nlat;
  nxp = nx + 1;
  nyp = ny + 1;

  x = (double *)malloc(nxp*nyp*sizeof(double));
  y = (double *)malloc(nxp*nyp*sizeof(double));
  z = (double *)malloc(nxp*nyp*sizeof(double));

  latlon2xyz(nxp*nyp, lon, lat, x, y, z);
  
  for(j=0; j<ny; j++) for(i=0; i < nx; i++) {
    /* clockwise */
    n0 = j*nxp+i;
    n1 = (j+1)*nxp+i;
    n2 = (j+1)*nxp+i+1;
    n3 = j*nxp+i+1;
    rewindList();
    grid = getNext();
    addEnd(grid, x[n0], y[n0], z[n0], 0, 0, 0, -1);
    addEnd(grid, x[n1], y[n1], z[n1], 0, 0, 0, -1);
    addEnd(grid, x[n2], y[n2], z[n2], 0, 0, 0, -1);
    addEnd(grid, x[n3], y[n3], z[n3], 0, 0, 0, -1);
    area[j*nx+i] = gridArea(grid);
  }

  free(x);
  free(y);
  free(z);
  
};  /* get_grid_great_circle_area */


void get_grid_area_dimensionless(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area)
{
  int nx, ny, nxp, i, j, n_in;
  double x_in[20], y_in[20];
  
  nx = *nlon;
  ny = *nlat;
  nxp = nx + 1;

  for(j=0; j<ny; j++) for(i=0; i < nx; i++) {
    x_in[0] = lon[j*nxp+i];
    x_in[1] = lon[j*nxp+i+1];
    x_in[2] = lon[(j+1)*nxp+i+1];
    x_in[3] = lon[(j+1)*nxp+i];
    y_in[0] = lat[j*nxp+i];
    y_in[1] = lat[j*nxp+i+1];
    y_in[2] = lat[(j+1)*nxp+i+1];
    y_in[3] = lat[(j+1)*nxp+i];
    n_in = fix_lon(x_in, y_in, 4, M_PI);    
    area[j*nx+i] = poly_area_dimensionless(x_in, y_in, n_in);
  }

};  /* get_grid_area */



void get_grid_area_no_adjust(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area)
{
  int nx, ny, nxp, i, j, n_in;
  double x_in[20], y_in[20];
  
  nx = *nlon;
  ny = *nlat;
  nxp = nx + 1;

  for(j=0; j<ny; j++) for(i=0; i < nx; i++) {
    x_in[0] = lon[j*nxp+i];
    x_in[1] = lon[j*nxp+i+1];
    x_in[2] = lon[(j+1)*nxp+i+1];
    x_in[3] = lon[(j+1)*nxp+i];
    y_in[0] = lat[j*nxp+i];
    y_in[1] = lat[j*nxp+i+1];
    y_in[2] = lat[(j+1)*nxp+i+1];
    y_in[3] = lat[(j+1)*nxp+i];
    n_in = 4;
    area[j*nx+i] = poly_area_no_adjust(x_in, y_in, n_in);
  }

};  /* get_grid_area_no_adjust */

/*******************************************************************************
  void create_xgrid_1dx2d_order1
  This routine generate exchange grids between two grids for the first order
  conservative interpolation. nlon_in,nlat_in,nlon_out,nlat_out are the size of the grid cell
  and lon_in,lat_in are 1-D grid bounds, lon_out,lat_out are geographic grid location of grid cell bounds. 
*******************************************************************************/
int create_xgrid_1dx2d_order1_(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
			       const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
			       const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out, double *xgrid_area)
{
  int nxgrid;
  
  nxgrid = create_xgrid_1dx2d_order1(nlon_in, nlat_in, nlon_out, nlat_out, lon_in, lat_in, lon_out, lat_out, mask_in,
			       i_in, j_in, i_out, j_out, xgrid_area);
  return nxgrid;
    
};  

int create_xgrid_1dx2d_order1(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out, const double *lon_in,
			      const double *lat_in, const double *lon_out, const double *lat_out,
			      const double *mask_in, int *i_in, int *j_in, int *i_out,
			      int *j_out, double *xgrid_area)
{

  int nx1, ny1, nx2, ny2, nx1p, nx2p;
  int i1, j1, i2, j2, nxgrid;
  double ll_lon, ll_lat, ur_lon, ur_lat, x_in[MV], y_in[MV], x_out[MV], y_out[MV];
  double *area_in, *area_out, min_area;
  double *tmpx, *tmpy;
  
  nx1 = *nlon_in;
  ny1 = *nlat_in;
  nx2 = *nlon_out;
  ny2 = *nlat_out;

  nxgrid = 0;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;

  area_in = (double *)malloc(nx1*ny1*sizeof(double));
  area_out = (double *)malloc(nx2*ny2*sizeof(double));
  tmpx = (double *)malloc((nx1+1)*(ny1+1)*sizeof(double));
  tmpy = (double *)malloc((nx1+1)*(ny1+1)*sizeof(double));
  for(j1=0; j1<=ny1; j1++) for(i1=0; i1<=nx1; i1++) {
    tmpx[j1*nx1p+i1] = lon_in[i1];
    tmpy[j1*nx1p+i1] = lat_in[j1];
  }
  /* This is just a temporary fix to solve the issue that there is one point in zonal direction */
  if(nx1 > 1)
     get_grid_area(nlon_in, nlat_in, tmpx, tmpy, area_in);
  else
    get_grid_area_no_adjust(nlon_in, nlat_in, tmpx, tmpy, area_in);
  
  get_grid_area(nlon_out, nlat_out, lon_out, lat_out, area_out);  
  free(tmpx);
  free(tmpy);  

  for(j1=0; j1<ny1; j1++) for(i1=0; i1<nx1; i1++) if( mask_in[j1*nx1+i1] > MASK_THRESH ) {

    ll_lon = lon_in[i1];   ll_lat = lat_in[j1];
    ur_lon = lon_in[i1+1]; ur_lat = lat_in[j1+1];
    for(j2=0; j2<ny2; j2++) for(i2=0; i2<nx2; i2++) {
      int n_in, n_out;
      double Xarea;
      
      y_in[0] = lat_out[j2*nx2p+i2];
      y_in[1] = lat_out[j2*nx2p+i2+1];
      y_in[2] = lat_out[(j2+1)*nx2p+i2+1];
      y_in[3] = lat_out[(j2+1)*nx2p+i2];
      if (  (y_in[0]<=ll_lat) && (y_in[1]<=ll_lat)
	    && (y_in[2]<=ll_lat) && (y_in[3]<=ll_lat) ) continue;
      if (  (y_in[0]>=ur_lat) && (y_in[1]>=ur_lat)
	    && (y_in[2]>=ur_lat) && (y_in[3]>=ur_lat) ) continue;

      x_in[0] = lon_out[j2*nx2p+i2];
      x_in[1] = lon_out[j2*nx2p+i2+1];
      x_in[2] = lon_out[(j2+1)*nx2p+i2+1];
      x_in[3] = lon_out[(j2+1)*nx2p+i2];
      n_in = fix_lon(x_in, y_in, 4, (ll_lon+ur_lon)/2);
      
      if ( (n_out = clip ( x_in, y_in, n_in, ll_lon, ll_lat, ur_lon, ur_lat, x_out, y_out )) > 0 ) {
	Xarea = poly_area (x_out, y_out, n_out ) * mask_in[j1*nx1+i1];
	min_area = min(area_in[j1*nx1+i1], area_out[j2*nx2+i2]);
	if( Xarea/min_area > AREA_RATIO_THRESH ) {
      	  xgrid_area[nxgrid] = Xarea;
	  i_in[nxgrid]    = i1;
	  j_in[nxgrid]    = j1;
	  i_out[nxgrid]   = i2;
	  j_out[nxgrid]   = j2;
	  ++nxgrid;
	  if(nxgrid > MAXXGRID) error_handler("nxgrid is greater than MAXXGRID, increase MAXXGRID");
	}
      }
    }
  }

  free(area_in);
  free(area_out);
  
  return nxgrid;
  
}; /* create_xgrid_1dx2d_order1 */


/********************************************************************************
  void create_xgrid_1dx2d_order2
  This routine generate exchange grids between two grids for the second order
  conservative interpolation. nlon_in,nlat_in,nlon_out,nlat_out are the size of the grid cell
  and lon_in,lat_in are 1-D grid bounds, lon_out,lat_out are geographic grid location of grid cell bounds.
********************************************************************************/
int create_xgrid_1dx2d_order2_(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
			       const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
			       const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
			       double *xgrid_area, double *xgrid_clon, double *xgrid_clat)
{
  int nxgrid;
  nxgrid = create_xgrid_1dx2d_order2(nlon_in, nlat_in, nlon_out, nlat_out, lon_in, lat_in, lon_out, lat_out, mask_in, i_in,
                                     j_in, i_out, j_out, xgrid_area, xgrid_clon, xgrid_clat);
  return nxgrid;

};
int create_xgrid_1dx2d_order2(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
			      const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
			      const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
			      double *xgrid_area, double *xgrid_clon, double *xgrid_clat)
{

  int nx1, ny1, nx2, ny2, nx1p, nx2p;
  int i1, j1, i2, j2, nxgrid, n;
  double ll_lon, ll_lat, ur_lon, ur_lat, x_in[MV], y_in[MV], x_out[MV], y_out[MV];
  double *area_in, *area_out, min_area;
  double *tmpx, *tmpy;
  
  nx1 = *nlon_in;
  ny1 = *nlat_in;
  nx2 = *nlon_out;
  ny2 = *nlat_out;

  nxgrid = 0;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;

  area_in      = (double *)malloc(nx1*ny1*sizeof(double));
  area_out     = (double *)malloc(nx2*ny2*sizeof(double));
  tmpx = (double *)malloc((nx1+1)*(ny1+1)*sizeof(double));
  tmpy = (double *)malloc((nx1+1)*(ny1+1)*sizeof(double));
  for(j1=0; j1<=ny1; j1++) for(i1=0; i1<=nx1; i1++) {
    tmpx[j1*nx1p+i1] = lon_in[i1];
    tmpy[j1*nx1p+i1] = lat_in[j1];
  }
  get_grid_area(nlon_in, nlat_in, tmpx, tmpy, area_in);     
  get_grid_area(nlon_out, nlat_out, lon_out, lat_out, area_out);  
  free(tmpx);
  free(tmpy);    
  
  for(j1=0; j1<ny1; j1++) for(i1=0; i1<nx1; i1++) if( mask_in[j1*nx1+i1] > MASK_THRESH ) {

    ll_lon = lon_in[i1];   ll_lat = lat_in[j1];
    ur_lon = lon_in[i1+1]; ur_lat = lat_in[j1+1];
    for(j2=0; j2<ny2; j2++) for(i2=0; i2<nx2; i2++) {
      int n_in, n_out;
      double xarea, lon_in_avg;
      
      y_in[0] = lat_out[j2*nx2p+i2];
      y_in[1] = lat_out[j2*nx2p+i2+1];
      y_in[2] = lat_out[(j2+1)*nx2p+i2+1];
      y_in[3] = lat_out[(j2+1)*nx2p+i2];
      if (  (y_in[0]<=ll_lat) && (y_in[1]<=ll_lat)
	    && (y_in[2]<=ll_lat) && (y_in[3]<=ll_lat) ) continue;
      if (  (y_in[0]>=ur_lat) && (y_in[1]>=ur_lat)
	    && (y_in[2]>=ur_lat) && (y_in[3]>=ur_lat) ) continue;

      x_in[0] = lon_out[j2*nx2p+i2];
      x_in[1] = lon_out[j2*nx2p+i2+1];
      x_in[2] = lon_out[(j2+1)*nx2p+i2+1];
      x_in[3] = lon_out[(j2+1)*nx2p+i2];
      n_in = fix_lon(x_in, y_in, 4, (ll_lon+ur_lon)/2);
      lon_in_avg = avgval_double(n_in, x_in);
      
      if (  (n_out = clip ( x_in, y_in, n_in, ll_lon, ll_lat, ur_lon, ur_lat, x_out, y_out )) > 0 ) {
	xarea = poly_area (x_out, y_out, n_out ) * mask_in[j1*nx1+i1];	
        min_area = min(area_in[j1*nx1+i1], area_out[j2*nx2+i2]);
	if(xarea/min_area > AREA_RATIO_THRESH ) {	  
	  xgrid_area[nxgrid] = xarea;
	  xgrid_clon[nxgrid] = poly_ctrlon(x_out, y_out, n_out, lon_in_avg);
	  xgrid_clat[nxgrid] = poly_ctrlat (x_out, y_out, n_out );
	  i_in[nxgrid]    = i1;
	  j_in[nxgrid]    = j1;
	  i_out[nxgrid]   = i2;
	  j_out[nxgrid]   = j2;
	  ++nxgrid;
	  if(nxgrid > MAXXGRID) error_handler("nxgrid is greater than MAXXGRID, increase MAXXGRID");
	}
      }
    }
  }
  free(area_in);
  free(area_out);
  
  return nxgrid;
  
}; /* create_xgrid_1dx2d_order2 */

/*******************************************************************************
  void create_xgrid_2dx1d_order1
  This routine generate exchange grids between two grids for the first order
  conservative interpolation. nlon_in,nlat_in,nlon_out,nlat_out are the size of the grid cell
  and lon_out,lat_out are 1-D grid bounds, lon_in,lat_in are geographic grid location of grid cell bounds.
  mask is on grid lon_in/lat_in. 
*******************************************************************************/
int create_xgrid_2dx1d_order1_(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
			       const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
			       const double *mask_in, int *i_in, int *j_in, int *i_out,
			       int *j_out, double *xgrid_area)
{
  int nxgrid;
  
  nxgrid = create_xgrid_2dx1d_order1(nlon_in, nlat_in, nlon_out, nlat_out, lon_in, lat_in, lon_out, lat_out, mask_in,
			       i_in, j_in, i_out, j_out, xgrid_area);
  return nxgrid;
    
};  
int create_xgrid_2dx1d_order1(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out, const double *lon_in,
			      const double *lat_in, const double *lon_out, const double *lat_out,
			      const double *mask_in, int *i_in, int *j_in, int *i_out,
			      int *j_out, double *xgrid_area)
{

  int nx1, ny1, nx2, ny2, nx1p, nx2p;
  int i1, j1, i2, j2, nxgrid;
  double ll_lon, ll_lat, ur_lon, ur_lat, x_in[MV], y_in[MV], x_out[MV], y_out[MV];
  double *area_in, *area_out, min_area;
  double *tmpx, *tmpy;
  
  nx1 = *nlon_in;
  ny1 = *nlat_in;
  nx2 = *nlon_out;
  ny2 = *nlat_out;

  nxgrid = 0;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;
  area_in = (double *)malloc(nx1*ny1*sizeof(double));
  area_out = (double *)malloc(nx2*ny2*sizeof(double));
  tmpx = (double *)malloc((nx2+1)*(ny2+1)*sizeof(double));
  tmpy = (double *)malloc((nx2+1)*(ny2+1)*sizeof(double));
  for(j2=0; j2<=ny2; j2++) for(i2=0; i2<=nx2; i2++) {
    tmpx[j2*nx2p+i2] = lon_out[i2];
    tmpy[j2*nx2p+i2] = lat_out[j2];
  }
  get_grid_area(nlon_in, nlat_in, lon_in, lat_in, area_in);     
  get_grid_area(nlon_out, nlat_out, tmpx, tmpy, area_out);  

  free(tmpx);
  free(tmpy);
  
  for(j2=0; j2<ny2; j2++) for(i2=0; i2<nx2; i2++) {

    ll_lon = lon_out[i2];   ll_lat = lat_out[j2];
    ur_lon = lon_out[i2+1]; ur_lat = lat_out[j2+1];
    for(j1=0; j1<ny1; j1++) for(i1=0; i1<nx1; i1++) if( mask_in[j1*nx1+i1] > MASK_THRESH ) {
      int n_in, n_out;
      double Xarea;
      
      y_in[0] = lat_in[j1*nx1p+i1];
      y_in[1] = lat_in[j1*nx1p+i1+1];
      y_in[2] = lat_in[(j1+1)*nx1p+i1+1];
      y_in[3] = lat_in[(j1+1)*nx1p+i1];
      if (  (y_in[0]<=ll_lat) && (y_in[1]<=ll_lat)
	    && (y_in[2]<=ll_lat) && (y_in[3]<=ll_lat) ) continue;
      if (  (y_in[0]>=ur_lat) && (y_in[1]>=ur_lat)
	    && (y_in[2]>=ur_lat) && (y_in[3]>=ur_lat) ) continue;

      x_in[0] = lon_in[j1*nx1p+i1];
      x_in[1] = lon_in[j1*nx1p+i1+1];
      x_in[2] = lon_in[(j1+1)*nx1p+i1+1];
      x_in[3] = lon_in[(j1+1)*nx1p+i1];

      n_in = fix_lon(x_in, y_in, 4, (ll_lon+ur_lon)/2);
      
      if ( (n_out = clip ( x_in, y_in, n_in, ll_lon, ll_lat, ur_lon, ur_lat, x_out, y_out )) > 0 ) {
	Xarea = poly_area ( x_out, y_out, n_out ) * mask_in[j1*nx1+i1];
	min_area = min(area_in[j1*nx1+i1], area_out[j2*nx2+i2]);
	if( Xarea/min_area > AREA_RATIO_THRESH ) {
      	  xgrid_area[nxgrid] = Xarea;
	  i_in[nxgrid]    = i1;
	  j_in[nxgrid]    = j1;
	  i_out[nxgrid]   = i2;
	  j_out[nxgrid]   = j2;
	  ++nxgrid;
	  if(nxgrid > MAXXGRID) error_handler("nxgrid is greater than MAXXGRID, increase MAXXGRID");
	}
      }
    }
  }

  free(area_in);
  free(area_out);
  
  return nxgrid;
  
}; /* create_xgrid_2dx1d_order1 */


/********************************************************************************
  void create_xgrid_2dx1d_order2
  This routine generate exchange grids between two grids for the second order
  conservative interpolation. nlon_in,nlat_in,nlon_out,nlat_out are the size of the grid cell
  and lon_out,lat_out are 1-D grid bounds, lon_in,lat_in are geographic grid location of grid cell bounds.
  mask is on grid lon_in/lat_in. 
********************************************************************************/
int create_xgrid_2dx1d_order2_(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
			       const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
			       const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
			       double *xgrid_area, double *xgrid_clon, double *xgrid_clat)
{
  int nxgrid;
  nxgrid = create_xgrid_2dx1d_order2(nlon_in, nlat_in, nlon_out, nlat_out, lon_in, lat_in, lon_out, lat_out, mask_in, i_in,
                                     j_in, i_out, j_out, xgrid_area, xgrid_clon, xgrid_clat);
  return nxgrid;

};

int create_xgrid_2dx1d_order2(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
			      const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
			      const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
			      double *xgrid_area, double *xgrid_clon, double *xgrid_clat)
{

  int nx1, ny1, nx2, ny2, nx1p, nx2p;
  int i1, j1, i2, j2, nxgrid, n;
  double ll_lon, ll_lat, ur_lon, ur_lat, x_in[MV], y_in[MV], x_out[MV], y_out[MV];
  double *tmpx, *tmpy;
  double *area_in, *area_out, min_area;
  double  lon_in_avg;
  
  nx1 = *nlon_in;
  ny1 = *nlat_in;
  nx2 = *nlon_out;
  ny2 = *nlat_out;

  nxgrid = 0;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;

  area_in      = (double *)malloc(nx1*ny1*sizeof(double));
  area_out     = (double *)malloc(nx2*ny2*sizeof(double));
  tmpx = (double *)malloc((nx2+1)*(ny2+1)*sizeof(double));
  tmpy = (double *)malloc((nx2+1)*(ny2+1)*sizeof(double));
  for(j2=0; j2<=ny2; j2++) for(i2=0; i2<=nx2; i2++) {
    tmpx[j2*nx2p+i2] = lon_out[i2];
    tmpy[j2*nx2p+i2] = lat_out[j2];
  }
  get_grid_area(nlon_in, nlat_in, lon_in, lat_in, area_in);     
  get_grid_area(nlon_out, nlat_out, tmpx, tmpy, area_out);  

  free(tmpx);
  free(tmpy);  
  
  for(j2=0; j2<ny2; j2++) for(i2=0; i2<nx2; i2++) {

    ll_lon = lon_out[i2];   ll_lat = lat_out[j2];
    ur_lon = lon_out[i2+1]; ur_lat = lat_out[j2+1];
    for(j1=0; j1<ny1; j1++) for(i1=0; i1<nx1; i1++) if( mask_in[j1*nx1+i1] > MASK_THRESH ) {
      int n_in, n_out;
      double xarea;
      
      y_in[0] = lat_in[j1*nx1p+i1];
      y_in[1] = lat_in[j1*nx1p+i1+1];
      y_in[2] = lat_in[(j1+1)*nx1p+i1+1];
      y_in[3] = lat_in[(j1+1)*nx1p+i1];
      if (  (y_in[0]<=ll_lat) && (y_in[1]<=ll_lat)
	    && (y_in[2]<=ll_lat) && (y_in[3]<=ll_lat) ) continue;
      if (  (y_in[0]>=ur_lat) && (y_in[1]>=ur_lat)
	    && (y_in[2]>=ur_lat) && (y_in[3]>=ur_lat) ) continue;

      x_in[0] = lon_in[j1*nx1p+i1];
      x_in[1] = lon_in[j1*nx1p+i1+1];
      x_in[2] = lon_in[(j1+1)*nx1p+i1+1];
      x_in[3] = lon_in[(j1+1)*nx1p+i1];

      n_in = fix_lon(x_in, y_in, 4, (ll_lon+ur_lon)/2);
      lon_in_avg = avgval_double(n_in, x_in);
      
      if (  (n_out = clip ( x_in, y_in, n_in, ll_lon, ll_lat, ur_lon, ur_lat, x_out, y_out )) > 0 ) {
	xarea = poly_area (x_out, y_out, n_out ) * mask_in[j1*nx1+i1];	
	min_area = min(area_in[j1*nx1+i1], area_out[j2*nx2+i2]);
	if(xarea/min_area > AREA_RATIO_THRESH ) {	  
	  xgrid_area[nxgrid] = xarea;
	  xgrid_clon[nxgrid] = poly_ctrlon(x_out, y_out, n_out, lon_in_avg);
	  xgrid_clat[nxgrid] = poly_ctrlat (x_out, y_out, n_out );
	  i_in[nxgrid]  = i1;
	  j_in[nxgrid]  = j1;
	  i_out[nxgrid] = i2;
	  j_out[nxgrid] = j2;
	  ++nxgrid;
	  if(nxgrid > MAXXGRID) error_handler("nxgrid is greater than MAXXGRID, increase MAXXGRID");
	}
      }
    }
  }

  free(area_in);
  free(area_out);  
  
  return nxgrid;
  
}; /* create_xgrid_2dx1d_order2 */

/*******************************************************************************
  void create_xgrid_2DX2D_order1
  This routine generate exchange grids between two grids for the first order
  conservative interpolation. nlon_in,nlat_in,nlon_out,nlat_out are the size of the grid cell
  and lon_in,lat_in, lon_out,lat_out are geographic grid location of grid cell bounds.
  mask is on grid lon_in/lat_in.
*******************************************************************************/
#ifndef __AIX
int create_xgrid_2dx2d_order1_(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
			       const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
			       const double *mask_in, int *i_in, int *j_in, int *i_out,
			       int *j_out, double *xgrid_area)
{
  int nxgrid;
  
  nxgrid = create_xgrid_2dx2d_order1(nlon_in, nlat_in, nlon_out, nlat_out, lon_in, lat_in, lon_out, lat_out, mask_in,
			       i_in, j_in, i_out, j_out, xgrid_area);
  return nxgrid;
    
};  
#endif
int create_xgrid_2dx2d_order1(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
			      const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
			      const double *mask_in, int *i_in, int *j_in, int *i_out,
			      int *j_out, double *xgrid_area)
{

#define MAX_V 8  
  int nx1, nx2, ny1, ny2, nx1p, nx2p, nxgrid;
  double *area_in, *area_out;
  int nblocks =1;
  int *istart2=NULL, *iend2=NULL;
  int npts_left, nblks_left, pos, m, npts_my, ij;
  double *lon_out_min_list,*lon_out_max_list,*lon_out_avg,*lat_out_min_list,*lat_out_max_list;  
  double *lon_out_list, *lat_out_list;
  int *pnxgrid=NULL, *pstart;
  int *pi_in=NULL, *pj_in=NULL, *pi_out=NULL, *pj_out=NULL;
  double *pxgrid_area=NULL;
  int    *n2_list;
  int nthreads, nxgrid_block_max;
  
  nx1 = *nlon_in;
  ny1 = *nlat_in;
  nx2 = *nlon_out;
  ny2 = *nlat_out;  
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;

  area_in  = (double *)malloc(nx1*ny1*sizeof(double));
  area_out = (double *)malloc(nx2*ny2*sizeof(double));
  get_grid_area(nlon_in, nlat_in, lon_in, lat_in, area_in);     
  get_grid_area(nlon_out, nlat_out, lon_out, lat_out, area_out);

  nthreads = 1;
#if defined(_OPENMP)
#pragma omp parallel
  nthreads = omp_get_num_threads();
#endif  

  nblocks = nthreads;

  istart2 = (int *)malloc(nblocks*sizeof(int));
  iend2 = (int *)malloc(nblocks*sizeof(int));

  pstart = (int *)malloc(nblocks*sizeof(int));
  pnxgrid = (int *)malloc(nblocks*sizeof(int));

  nxgrid_block_max = MAXXGRID/nblocks;
  
  for(m=0; m<nblocks; m++) {
    pnxgrid[m] = 0;
    pstart[m] = m*nxgrid_block_max;
  }

  if(nblocks == 1) {
    pi_in = i_in;
    pj_in = j_in;
    pi_out = i_out;
    pj_out = j_out;
    pxgrid_area = xgrid_area;
  }
  else {
    pi_in = (int *)malloc(MAXXGRID*sizeof(int));
    pj_in = (int *)malloc(MAXXGRID*sizeof(int));
    pi_out = (int *)malloc(MAXXGRID*sizeof(int));
    pj_out = (int *)malloc(MAXXGRID*sizeof(int));
    pxgrid_area = (double *)malloc(MAXXGRID*sizeof(double));
  }
  
  npts_left = nx2*ny2;
  nblks_left = nblocks;
  pos = 0;
  for(m=0; m<nblocks; m++) {
    istart2[m] = pos;
    npts_my = npts_left/nblks_left;
    iend2[m] = istart2[m] + npts_my - 1;
    pos = iend2[m] + 1;
    npts_left -= npts_my;
    nblks_left--;
  }

  lon_out_min_list = (double *)malloc(nx2*ny2*sizeof(double));
  lon_out_max_list = (double *)malloc(nx2*ny2*sizeof(double));
  lat_out_min_list = (double *)malloc(nx2*ny2*sizeof(double));
  lat_out_max_list = (double *)malloc(nx2*ny2*sizeof(double));
  lon_out_avg = (double *)malloc(nx2*ny2*sizeof(double));
  n2_list     = (int *)malloc(nx2*ny2*sizeof(int));
  lon_out_list = (double *)malloc(MAX_V*nx2*ny2*sizeof(double));
  lat_out_list = (double *)malloc(MAX_V*nx2*ny2*sizeof(double));
#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(nx2,ny2,nx2p,lon_out,lat_out,lat_out_min_list, \
                                              lat_out_max_list,lon_out_min_list,lon_out_max_list, \
                                              lon_out_avg,n2_list,lon_out_list,lat_out_list) 
#endif                        
  for(ij=0; ij<nx2*ny2; ij++){
    int i2, j2, n, n0, n1, n2, n3, n2_in, l;
    double x2_in[MV], y2_in[MV];
    i2 = ij%nx2;
    j2 = ij/nx2;
    n = j2*nx2+i2;
    n0 = j2*nx2p+i2; n1 = j2*nx2p+i2+1;
    n2 = (j2+1)*nx2p+i2+1; n3 = (j2+1)*nx2p+i2;
    x2_in[0] = lon_out[n0]; y2_in[0] = lat_out[n0];
    x2_in[1] = lon_out[n1]; y2_in[1] = lat_out[n1];
    x2_in[2] = lon_out[n2]; y2_in[2] = lat_out[n2];
    x2_in[3] = lon_out[n3]; y2_in[3] = lat_out[n3];
    
    lat_out_min_list[n] = minval_double(4, y2_in);
    lat_out_max_list[n] = maxval_double(4, y2_in);
    n2_in = fix_lon(x2_in, y2_in, 4, M_PI);
    if(n2_in > MAX_V) error_handler("create_xgrid.c: n2_in is greater than MAX_V");
    lon_out_min_list[n] = minval_double(n2_in, x2_in);
    lon_out_max_list[n] = maxval_double(n2_in, x2_in);
    lon_out_avg[n] = avgval_double(n2_in, x2_in);
    n2_list[n] = n2_in;
    for(l=0; l<n2_in; l++) {
      lon_out_list[n*MAX_V+l] = x2_in[l];
      lat_out_list[n*MAX_V+l] = y2_in[l];
    }    
  }

nxgrid = 0;

#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(nblocks,nx1,ny1,nx1p,mask_in,lon_in,lat_in, \
                                              istart2,iend2,nx2,lat_out_min_list,lat_out_max_list, \
                                              n2_list,lon_out_list,lat_out_list,lon_out_min_list, \
                                              lon_out_max_list,lon_out_avg,area_in,area_out, \
                                              pxgrid_area,pnxgrid,pi_in,pj_in,pi_out,pj_out,pstart,nthreads)
#endif  
  for(m=0; m<nblocks; m++) {
    int i1, j1, ij;
    for(j1=0; j1<ny1; j1++) for(i1=0; i1<nx1; i1++) if( mask_in[j1*nx1+i1] > MASK_THRESH ) {
      int n0, n1, n2, n3, l,n1_in;
      double lat_in_min,lat_in_max,lon_in_min,lon_in_max,lon_in_avg;
      double x1_in[MV], y1_in[MV], x_out[MV], y_out[MV];
 
      n0 = j1*nx1p+i1;       n1 = j1*nx1p+i1+1;
      n2 = (j1+1)*nx1p+i1+1; n3 = (j1+1)*nx1p+i1;      
      x1_in[0] = lon_in[n0]; y1_in[0] = lat_in[n0];
      x1_in[1] = lon_in[n1]; y1_in[1] = lat_in[n1];
      x1_in[2] = lon_in[n2]; y1_in[2] = lat_in[n2];
      x1_in[3] = lon_in[n3]; y1_in[3] = lat_in[n3];
      lat_in_min = minval_double(4, y1_in);
      lat_in_max = maxval_double(4, y1_in);
      n1_in = fix_lon(x1_in, y1_in, 4, M_PI);
      lon_in_min = minval_double(n1_in, x1_in);
      lon_in_max = maxval_double(n1_in, x1_in);
      lon_in_avg = avgval_double(n1_in, x1_in);
      for(ij=istart2[m]; ij<=iend2[m]; ij++) {
	int n_in, n_out, i2, j2, n2_in;
	double xarea, dx, lon_out_min, lon_out_max;
	double x2_in[MAX_V], y2_in[MAX_V];
	
	i2 = ij%nx2;
	j2 = ij/nx2;
	
	if(lat_out_min_list[ij] >= lat_in_max || lat_out_max_list[ij] <= lat_in_min ) continue;
	/* adjust x2_in according to lon_in_avg*/
	n2_in = n2_list[ij];
	for(l=0; l<n2_in; l++) {
	  x2_in[l] = lon_out_list[ij*MAX_V+l];
	  y2_in[l] = lat_out_list[ij*MAX_V+l];
	}
	lon_out_min = lon_out_min_list[ij];
	lon_out_max = lon_out_max_list[ij];  
        dx = lon_out_avg[ij] - lon_in_avg;
	if(dx < -M_PI ) {
	  lon_out_min += TPI;
	  lon_out_max += TPI;
	  for (l=0; l<n2_in; l++) x2_in[l] += TPI;
	}
        else if (dx >  M_PI) {
	  lon_out_min -= TPI;
	  lon_out_max -= TPI;
	  for (l=0; l<n2_in; l++) x2_in[l] -= TPI;
	}

	/* x2_in should in the same range as x1_in after lon_fix, so no need to
	   consider cyclic condition
	*/
	if(lon_out_min >= lon_in_max || lon_out_max <= lon_in_min ) continue;
	if (  (n_out = clip_2dx2d( x1_in, y1_in, n1_in, x2_in, y2_in, n2_in, x_out, y_out )) > 0) {
          double min_area;
	  int    nn;
	  xarea = poly_area (x_out, y_out, n_out ) * mask_in[j1*nx1+i1];	
	  min_area = min(area_in[j1*nx1+i1], area_out[j2*nx2+i2]);
	  if( xarea/min_area > AREA_RATIO_THRESH ) {
	    pnxgrid[m]++;
            if(pnxgrid[m]>= MAXXGRID/nthreads)
	      error_handler("nxgrid is greater than MAXXGRID/nthreads, increase MAXXGRID, decrease nthreads, or increase number of MPI ranks");
	    nn = pstart[m] + pnxgrid[m]-1;	    

	    pxgrid_area[nn] = xarea;
	    pi_in[nn]       = i1;
	    pj_in[nn]       = j1;
	    pi_out[nn]      = i2;
	    pj_out[nn]      = j2;
	  }
	  
	}
	
      }
    }
  }

  /*copy data if nblocks > 1 */
  if(nblocks == 1) {
    nxgrid = pnxgrid[0];
    pi_in = NULL;
    pj_in = NULL;
    pi_out = NULL;
    pj_out = NULL;
    pxgrid_area = NULL;
  }
  else {
    int nn, i;
    nxgrid = 0;
    for(m=0; m<nblocks; m++) {
      for(i=0; i<pnxgrid[m]; i++) {
	nn = pstart[m] + i;
	i_in[nxgrid] = pi_in[nn];
	j_in[nxgrid] = pj_in[nn];
	i_out[nxgrid] = pi_out[nn];
	j_out[nxgrid] = pj_out[nn];
	xgrid_area[nxgrid] = pxgrid_area[nn];
	nxgrid++;
      }
    }
    free(pi_in);
    free(pj_in);
    free(pi_out);
    free(pj_out);
    free(pxgrid_area);
  }
  
  free(area_in);
  free(area_out);  
  free(lon_out_min_list);
  free(lon_out_max_list);
  free(lat_out_min_list);
  free(lat_out_max_list);
  free(lon_out_avg);
  free(n2_list);  
  free(lon_out_list);
  free(lat_out_list);

  return nxgrid;
  
};/* get_xgrid_2Dx2D_order1 */

/********************************************************************************
  void create_xgrid_2dx1d_order2
  This routine generate exchange grids between two grids for the second order
  conservative interpolation. nlon_in,nlat_in,nlon_out,nlat_out are the size of the grid cell
  and lon_in,lat_in, lon_out,lat_out are geographic grid location of grid cell bounds.
  mask is on grid lon_in/lat_in. 
********************************************************************************/
#ifndef __AIX
int create_xgrid_2dx2d_order2_(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
			       const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
			       const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
			       double *xgrid_area, double *xgrid_clon, double *xgrid_clat)
{
  int nxgrid;
  nxgrid = create_xgrid_2dx2d_order2(nlon_in, nlat_in, nlon_out, nlat_out, lon_in, lat_in, lon_out, lat_out, mask_in, i_in,
                                     j_in, i_out, j_out, xgrid_area, xgrid_clon, xgrid_clat);
  return nxgrid;

};
#endif
int create_xgrid_2dx2d_order2(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
			      const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
			      const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
			      double *xgrid_area, double *xgrid_clon, double *xgrid_clat)
{

#define MAX_V 8  
  int nx1, nx2, ny1, ny2, nx1p, nx2p, nxgrid;
  double xctrlon, xctrlat;
  double *area_in, *area_out;
  int nblocks =1;
  int *istart2=NULL, *iend2=NULL;
  int npts_left, nblks_left, pos, m, npts_my, ij;
  double *lon_out_min_list,*lon_out_max_list,*lon_out_avg,*lat_out_min_list,*lat_out_max_list;  
  double *lon_out_list, *lat_out_list;
  int *pnxgrid=NULL, *pstart;
  int *pi_in=NULL, *pj_in=NULL, *pi_out=NULL, *pj_out=NULL;
  double *pxgrid_area=NULL, *pxgrid_clon=NULL, *pxgrid_clat=NULL;
  int    *n2_list;
  int nthreads, nxgrid_block_max;
  
  nx1 = *nlon_in;
  ny1 = *nlat_in;
  nx2 = *nlon_out;
  ny2 = *nlat_out;  
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;

  area_in  = (double *)malloc(nx1*ny1*sizeof(double));
  area_out = (double *)malloc(nx2*ny2*sizeof(double));
  get_grid_area(nlon_in, nlat_in, lon_in, lat_in, area_in);     
  get_grid_area(nlon_out, nlat_out, lon_out, lat_out, area_out);

  nthreads = 1;
#if defined(_OPENMP)
#pragma omp parallel
  nthreads = omp_get_num_threads();
#endif  

  nblocks = nthreads;

  istart2 = (int *)malloc(nblocks*sizeof(int));
  iend2 = (int *)malloc(nblocks*sizeof(int));

  pstart = (int *)malloc(nblocks*sizeof(int));
  pnxgrid = (int *)malloc(nblocks*sizeof(int));

  nxgrid_block_max = MAXXGRID/nblocks;
  
  for(m=0; m<nblocks; m++) {
    pnxgrid[m] = 0;
    pstart[m] = m*nxgrid_block_max;
  }

  if(nblocks == 1) {
    pi_in = i_in;
    pj_in = j_in;
    pi_out = i_out;
    pj_out = j_out;
    pxgrid_area = xgrid_area;
    pxgrid_clon = xgrid_clon;
    pxgrid_clat = xgrid_clat;
  }
  else {
    pi_in = (int *)malloc(MAXXGRID*sizeof(int));
    pj_in = (int *)malloc(MAXXGRID*sizeof(int));
    pi_out = (int *)malloc(MAXXGRID*sizeof(int));
    pj_out = (int *)malloc(MAXXGRID*sizeof(int));
    pxgrid_area = (double *)malloc(MAXXGRID*sizeof(double));
    pxgrid_clon = (double *)malloc(MAXXGRID*sizeof(double));
    pxgrid_clat = (double *)malloc(MAXXGRID*sizeof(double));    
  }
    
  npts_left = nx2*ny2;
  nblks_left = nblocks;
  pos = 0;
  for(m=0; m<nblocks; m++) {
    istart2[m] = pos;
    npts_my = npts_left/nblks_left;
    iend2[m] = istart2[m] + npts_my - 1;
    pos = iend2[m] + 1;
    npts_left -= npts_my;
    nblks_left--;
  }

  lon_out_min_list = (double *)malloc(nx2*ny2*sizeof(double));
  lon_out_max_list = (double *)malloc(nx2*ny2*sizeof(double));
  lat_out_min_list = (double *)malloc(nx2*ny2*sizeof(double));
  lat_out_max_list = (double *)malloc(nx2*ny2*sizeof(double));
  lon_out_avg = (double *)malloc(nx2*ny2*sizeof(double));
  n2_list     = (int *)malloc(nx2*ny2*sizeof(int));
  lon_out_list = (double *)malloc(MAX_V*nx2*ny2*sizeof(double));
  lat_out_list = (double *)malloc(MAX_V*nx2*ny2*sizeof(double));
#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(nx2,ny2,nx2p,lon_out,lat_out,lat_out_min_list, \
                                              lat_out_max_list,lon_out_min_list,lon_out_max_list, \
                                              lon_out_avg,n2_list,lon_out_list,lat_out_list) 
#endif                        
  for(ij=0; ij<nx2*ny2; ij++){
    int i2, j2, n, n0, n1, n2, n3, n2_in, l;
    double x2_in[MV], y2_in[MV];
    i2 = ij%nx2;
    j2 = ij/nx2;
    n = j2*nx2+i2;
    n0 = j2*nx2p+i2; n1 = j2*nx2p+i2+1;
    n2 = (j2+1)*nx2p+i2+1; n3 = (j2+1)*nx2p+i2;
    x2_in[0] = lon_out[n0]; y2_in[0] = lat_out[n0];
    x2_in[1] = lon_out[n1]; y2_in[1] = lat_out[n1];
    x2_in[2] = lon_out[n2]; y2_in[2] = lat_out[n2];
    x2_in[3] = lon_out[n3]; y2_in[3] = lat_out[n3];
    
    lat_out_min_list[n] = minval_double(4, y2_in);
    lat_out_max_list[n] = maxval_double(4, y2_in);
    n2_in = fix_lon(x2_in, y2_in, 4, M_PI);
    if(n2_in > MAX_V) error_handler("create_xgrid.c: n2_in is greater than MAX_V");
    lon_out_min_list[n] = minval_double(n2_in, x2_in);
    lon_out_max_list[n] = maxval_double(n2_in, x2_in);
    lon_out_avg[n] = avgval_double(n2_in, x2_in);
    n2_list[n] = n2_in;
    for(l=0; l<n2_in; l++) {
      lon_out_list[n*MAX_V+l] = x2_in[l];
      lat_out_list[n*MAX_V+l] = y2_in[l];
    }    
  }

nxgrid = 0;

#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(nblocks,nx1,ny1,nx1p,mask_in,lon_in,lat_in, \
                                              istart2,iend2,nx2,lat_out_min_list,lat_out_max_list, \
                                              n2_list,lon_out_list,lat_out_list,lon_out_min_list, \
                                              lon_out_max_list,lon_out_avg,area_in,area_out, \
                                              pxgrid_area,pnxgrid,pxgrid_clon,pxgrid_clat,pi_in, \
                                              pj_in,pi_out,pj_out,pstart,nthreads)
#endif  
  for(m=0; m<nblocks; m++) {
    int i1, j1, ij;
    for(j1=0; j1<ny1; j1++) for(i1=0; i1<nx1; i1++) if( mask_in[j1*nx1+i1] > MASK_THRESH ) {
      int n0, n1, n2, n3, l,n1_in;
      double lat_in_min,lat_in_max,lon_in_min,lon_in_max,lon_in_avg;
      double x1_in[MV], y1_in[MV], x_out[MV], y_out[MV];
 
      n0 = j1*nx1p+i1;       n1 = j1*nx1p+i1+1;
      n2 = (j1+1)*nx1p+i1+1; n3 = (j1+1)*nx1p+i1;      
      x1_in[0] = lon_in[n0]; y1_in[0] = lat_in[n0];
      x1_in[1] = lon_in[n1]; y1_in[1] = lat_in[n1];
      x1_in[2] = lon_in[n2]; y1_in[2] = lat_in[n2];
      x1_in[3] = lon_in[n3]; y1_in[3] = lat_in[n3];
      lat_in_min = minval_double(4, y1_in);
      lat_in_max = maxval_double(4, y1_in);
      n1_in = fix_lon(x1_in, y1_in, 4, M_PI);
      lon_in_min = minval_double(n1_in, x1_in);
      lon_in_max = maxval_double(n1_in, x1_in);
      lon_in_avg = avgval_double(n1_in, x1_in);
      for(ij=istart2[m]; ij<=iend2[m]; ij++) {
	int n_in, n_out, i2, j2, n2_in;
	double xarea, dx, lon_out_min, lon_out_max;
	double x2_in[MAX_V], y2_in[MAX_V];
	
	i2 = ij%nx2;
	j2 = ij/nx2;
	
	if(lat_out_min_list[ij] >= lat_in_max || lat_out_max_list[ij] <= lat_in_min ) continue;
	/* adjust x2_in according to lon_in_avg*/
	n2_in = n2_list[ij];
	for(l=0; l<n2_in; l++) {
	  x2_in[l] = lon_out_list[ij*MAX_V+l];
	  y2_in[l] = lat_out_list[ij*MAX_V+l];
	}
	lon_out_min = lon_out_min_list[ij];
	lon_out_max = lon_out_max_list[ij];  
        dx = lon_out_avg[ij] - lon_in_avg;
	if(dx < -M_PI ) {
	  lon_out_min += TPI;
	  lon_out_max += TPI;
	  for (l=0; l<n2_in; l++) x2_in[l] += TPI;
	}
        else if (dx >  M_PI) {
	  lon_out_min -= TPI;
	  lon_out_max -= TPI;
	  for (l=0; l<n2_in; l++) x2_in[l] -= TPI;
	}

	/* x2_in should in the same range as x1_in after lon_fix, so no need to
	   consider cyclic condition
	*/
	if(lon_out_min >= lon_in_max || lon_out_max <= lon_in_min ) continue;
	if (  (n_out = clip_2dx2d( x1_in, y1_in, n1_in, x2_in, y2_in, n2_in, x_out, y_out )) > 0) {
          double min_area;
	  int nn;
	  xarea = poly_area (x_out, y_out, n_out ) * mask_in[j1*nx1+i1];	
	  min_area = min(area_in[j1*nx1+i1], area_out[j2*nx2+i2]);
	  if( xarea/min_area > AREA_RATIO_THRESH ) {
	    pnxgrid[m]++;
            if(pnxgrid[m]>= MAXXGRID/nthreads)
	      error_handler("nxgrid is greater than MAXXGRID/nthreads, increase MAXXGRID, decrease nthreads, or increase number of MPI ranks");
	    nn = pstart[m] + pnxgrid[m]-1;
	    pxgrid_area[nn] = xarea;
	    pxgrid_clon[nn] = poly_ctrlon(x_out, y_out, n_out, lon_in_avg);
	    pxgrid_clat[nn] = poly_ctrlat (x_out, y_out, n_out );
	    pi_in[nn]       = i1;
	    pj_in[nn]       = j1;
	    pi_out[nn]      = i2;
	    pj_out[nn]      = j2;
	  }	  
	}
      }
    }
  }

  /*copy data if nblocks > 1 */
  if(nblocks == 1) {
    nxgrid = pnxgrid[0];
    pi_in = NULL;
    pj_in = NULL;
    pi_out = NULL;
    pj_out = NULL;
    pxgrid_area = NULL;
    pxgrid_clon = NULL;
    pxgrid_clat = NULL;
  }
  else {
    int nn, i;
    nxgrid = 0;
    for(m=0; m<nblocks; m++) {
      for(i=0; i<pnxgrid[m]; i++) {
	nn = pstart[m] + i;
	i_in[nxgrid] = pi_in[nn];
	j_in[nxgrid] = pj_in[nn];
	i_out[nxgrid] = pi_out[nn];
	j_out[nxgrid] = pj_out[nn];
	xgrid_area[nxgrid] = pxgrid_area[nn];
	xgrid_clon[nxgrid] = pxgrid_clon[nn];
	xgrid_clat[nxgrid] = pxgrid_clat[nn];
	nxgrid++;
      }
    }
    free(pi_in);
    free(pj_in);
    free(pi_out);
    free(pj_out);
    free(pxgrid_area);
    free(pxgrid_clon);
    free(pxgrid_clat);
  }

  free(area_in);
  free(area_out);  
  free(lon_out_min_list);
  free(lon_out_max_list);
  free(lat_out_min_list);
  free(lat_out_max_list);
  free(lon_out_avg);
  free(n2_list);  
  free(lon_out_list);
  free(lat_out_list);

  return nxgrid;
  
};/* get_xgrid_2Dx2D_order2 */


/*******************************************************************************
   Sutherland-Hodgeman algorithm sequentially clips parts outside 4 boundaries
*******************************************************************************/

int clip(const double lon_in[], const double lat_in[], int n_in, double ll_lon, double ll_lat,
	 double ur_lon, double ur_lat, double lon_out[], double lat_out[])
{
  double x_tmp[MV], y_tmp[MV], x_last, y_last;
  int i_in, i_out, n_out, inside_last, inside;

  /* clip polygon with LEFT boundary - clip V_IN to V_TMP */
  x_last = lon_in[n_in-1];
  y_last = lat_in[n_in-1];
  inside_last = (x_last >= ll_lon);
  for (i_in=0,i_out=0;i_in<n_in;i_in++) {
 
    /* if crossing LEFT boundary - output intersection */
    if ((inside=(lon_in[i_in] >= ll_lon))!=inside_last) {
      x_tmp[i_out] = ll_lon;
      y_tmp[i_out++] = y_last + (ll_lon - x_last) * (lat_in[i_in] - y_last) / (lon_in[i_in] - x_last);
    }

    /* if "to" point is right of LEFT boundary, output it */
    if (inside) {
      x_tmp[i_out]   = lon_in[i_in];
      y_tmp[i_out++] = lat_in[i_in];
    }
    x_last = lon_in[i_in];
    y_last = lat_in[i_in];
    inside_last = inside;
  }
  if (!(n_out=i_out)) return(0);

  /* clip polygon with RIGHT boundary - clip V_TMP to V_OUT */
  x_last = x_tmp[n_out-1];
  y_last = y_tmp[n_out-1];
  inside_last = (x_last <= ur_lon);
  for (i_in=0,i_out=0;i_in<n_out;i_in++) {
 
    /* if crossing RIGHT boundary - output intersection */
    if ((inside=(x_tmp[i_in] <= ur_lon))!=inside_last) {
      lon_out[i_out]   = ur_lon;
      lat_out[i_out++] = y_last + (ur_lon - x_last) * (y_tmp[i_in] - y_last)
                                                 / (x_tmp[i_in] - x_last);
    }

    /* if "to" point is left of RIGHT boundary, output it */
    if (inside) {
      lon_out[i_out]   = x_tmp[i_in];
      lat_out[i_out++] = y_tmp[i_in];
    }
    
    x_last = x_tmp[i_in];
    y_last = y_tmp[i_in];
    inside_last = inside;
  }
  if (!(n_out=i_out)) return(0);

  /* clip polygon with BOTTOM boundary - clip V_OUT to V_TMP */
  x_last = lon_out[n_out-1];
  y_last = lat_out[n_out-1];
  inside_last = (y_last >= ll_lat);
  for (i_in=0,i_out=0;i_in<n_out;i_in++) {
 
    /* if crossing BOTTOM boundary - output intersection */
    if ((inside=(lat_out[i_in] >= ll_lat))!=inside_last) {
      y_tmp[i_out]   = ll_lat;
      x_tmp[i_out++] = x_last + (ll_lat - y_last) * (lon_out[i_in] - x_last) / (lat_out[i_in] - y_last);
    }

    /* if "to" point is above BOTTOM boundary, output it */
    if (inside) {
      x_tmp[i_out]   = lon_out[i_in];
      y_tmp[i_out++] = lat_out[i_in];
    }
    x_last = lon_out[i_in];
    y_last = lat_out[i_in];
    inside_last = inside;
  }
  if (!(n_out=i_out)) return(0);

  /* clip polygon with TOP boundary - clip V_TMP to V_OUT */
  x_last = x_tmp[n_out-1];
  y_last = y_tmp[n_out-1];
  inside_last = (y_last <= ur_lat);
  for (i_in=0,i_out=0;i_in<n_out;i_in++) {
 
    /* if crossing TOP boundary - output intersection */
    if ((inside=(y_tmp[i_in] <= ur_lat))!=inside_last) {
      lat_out[i_out]   = ur_lat;
      lon_out[i_out++] = x_last + (ur_lat - y_last) * (x_tmp[i_in] - x_last) / (y_tmp[i_in] - y_last);
    }

    /* if "to" point is below TOP boundary, output it */
    if (inside) {
      lon_out[i_out]   = x_tmp[i_in];
      lat_out[i_out++] = y_tmp[i_in];
    }
    x_last = x_tmp[i_in];
    y_last = y_tmp[i_in];
    inside_last = inside;
  }
  return(i_out);
}; /* clip */  


/*******************************************************************************
   Revise Sutherland-Hodgeman algorithm to find the vertices of the overlapping
   between any two grid boxes. It return the number of vertices for the exchange grid.
*******************************************************************************/

int clip_2dx2d(const double lon1_in[], const double lat1_in[], int n1_in, 
	 const double lon2_in[], const double lat2_in[], int n2_in, 
	 double lon_out[], double lat_out[])
{
  double lon_tmp[MV], lat_tmp[MV];
  double x1_0, y1_0, x1_1, y1_1, x2_0, y2_0, x2_1, y2_1;
  double dx1, dy1, dx2, dy2, determ, ds1, ds2;
  int i_out, n_out, inside_last, inside, i1, i2;

  /* clip polygon with each boundary of the polygon */
  /* We treat lon1_in/lat1_in as clip polygon and lon2_in/lat2_in as subject polygon */
  n_out = n1_in;
  for(i1=0; i1<n1_in; i1++) {
    lon_tmp[i1] = lon1_in[i1];
    lat_tmp[i1] = lat1_in[i1];
  }
  x2_0 = lon2_in[n2_in-1];
  y2_0 = lat2_in[n2_in-1];
  for(i2=0; i2<n2_in; i2++) {
    x2_1 = lon2_in[i2];
    y2_1 = lat2_in[i2];
    x1_0 = lon_tmp[n_out-1];
    y1_0 = lat_tmp[n_out-1];
    inside_last = inside_edge( x2_0, y2_0, x2_1, y2_1, x1_0, y1_0);
    for(i1=0, i_out=0; i1<n_out; i1++) {
      x1_1 = lon_tmp[i1];
      y1_1 = lat_tmp[i1];
      if((inside = inside_edge(x2_0, y2_0, x2_1, y2_1, x1_1, y1_1)) != inside_last ) {
        /* there is intersection, the line between <x1_0,y1_0> and  <x1_1,y1_1>
           should not parallel to the line between <x2_0,y2_0> and  <x2_1,y2_1>
           may need to consider truncation error */
	dy1 = y1_1-y1_0;
	dy2 = y2_1-y2_0;
	dx1 = x1_1-x1_0;
	dx2 = x2_1-x2_0;
	ds1 = y1_0*x1_1 - y1_1*x1_0;
	ds2 = y2_0*x2_1 - y2_1*x2_0;
	determ = dy2*dx1 - dy1*dx2;
        if(fabs(determ) < EPSLN30) {
	  error_handler("the line between <x1_0,y1_0> and  <x1_1,y1_1> should not parallel to "
				     "the line between <x2_0,y2_0> and  <x2_1,y2_1>");
	}
	lon_out[i_out]   = (dx2*ds1 - dx1*ds2)/determ;
	lat_out[i_out++] = (dy2*ds1 - dy1*ds2)/determ;
	

      }
      if(inside) {
	lon_out[i_out]   = x1_1;
	lat_out[i_out++] = y1_1;	
      }
      x1_0 = x1_1;
      y1_0 = y1_1;
      inside_last = inside;
    }
    if(!(n_out=i_out)) return 0;
    for(i1=0; i1<n_out; i1++) {
      lon_tmp[i1] = lon_out[i1];
      lat_tmp[i1] = lat_out[i1];
    }
    /* shift the starting point */
    x2_0 = x2_1;
    y2_0 = y2_1;
  }
  return(n_out);
}; /* clip */
    
/*#define debug_test_create_xgrid*/  

#ifndef __AIX
int create_xgrid_great_circle_(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
			      const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
			      const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
			      double *xgrid_area, double *xgrid_clon, double *xgrid_clat)
{
  int nxgrid;
  nxgrid = create_xgrid_great_circle(nlon_in, nlat_in, nlon_out, nlat_out, lon_in, lat_in, lon_out, lat_out,
			      mask_in, i_in, j_in, i_out, j_out, xgrid_area, xgrid_clon, xgrid_clat);

  return nxgrid;
};
#endif
  
int create_xgrid_great_circle(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
			      const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
			      const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
			      double *xgrid_area, double *xgrid_clon, double *xgrid_clat)
{

  int nx1, nx2, ny1, ny2, nx1p, nx2p, ny1p, ny2p, nxgrid, n1_in, n2_in;
  int n0, n1, n2, n3, i1, j1, i2, j2, l, n;
  double x1_in[MV], y1_in[MV], z1_in[MV];
  double x2_in[MV], y2_in[MV], z2_in[MV];
  double x_out[MV], y_out[MV], z_out[MV];
  double *x1=NULL, *y1=NULL, *z1=NULL;
  double *x2=NULL, *y2=NULL, *z2=NULL;

  double xctrlon, xctrlat;
  double *area1, *area2, min_area;
  
  nx1 = *nlon_in;
  ny1 = *nlat_in;
  nx2 = *nlon_out;
  ny2 = *nlat_out;  
  nxgrid = 0;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;
  ny1p = ny1 + 1;
  ny2p = ny2 + 1;

  /* first convert lon-lat to cartesian coordinates */
  x1 = (double *)malloc(nx1p*ny1p*sizeof(double));
  y1 = (double *)malloc(nx1p*ny1p*sizeof(double));
  z1 = (double *)malloc(nx1p*ny1p*sizeof(double));
  x2 = (double *)malloc(nx2p*ny2p*sizeof(double));
  y2 = (double *)malloc(nx2p*ny2p*sizeof(double));
  z2 = (double *)malloc(nx2p*ny2p*sizeof(double));
  
  latlon2xyz(nx1p*ny1p, lon_in, lat_in, x1, y1, z1);
  latlon2xyz(nx2p*ny2p, lon_out, lat_out, x2, y2, z2);
  
  area1  = (double *)malloc(nx1*ny1*sizeof(double));
  area2 = (double *)malloc(nx2*ny2*sizeof(double));
  get_grid_great_circle_area(nlon_in, nlat_in, lon_in, lat_in, area1);     
  get_grid_great_circle_area(nlon_out, nlat_out, lon_out, lat_out, area2); 
  n1_in = 4;
  n2_in = 4;
  
  for(j1=0; j1<ny1; j1++) for(i1=0; i1<nx1; i1++) if( mask_in[j1*nx1+i1] > MASK_THRESH ) {
    /* clockwise */
    n0 = j1*nx1p+i1;       n1 = (j1+1)*nx1p+i1;
    n2 = (j1+1)*nx1p+i1+1; n3 = j1*nx1p+i1+1;      
    x1_in[0] = x1[n0]; y1_in[0] = y1[n0]; z1_in[0] = z1[n0];
    x1_in[1] = x1[n1]; y1_in[1] = y1[n1]; z1_in[1] = z1[n1];
    x1_in[2] = x1[n2]; y1_in[2] = y1[n2]; z1_in[2] = z1[n2];
    x1_in[3] = x1[n3]; y1_in[3] = y1[n3]; z1_in[3] = z1[n3];
      
    for(j2=0; j2<ny2; j2++) for(i2=0; i2<nx2; i2++) {
      int n_in, n_out;
      double xarea;

      n0 = j2*nx2p+i2;       n1 = (j2+1)*nx2p+i2;
      n2 = (j2+1)*nx2p+i2+1; n3 = j2*nx2p+i2+1;
      x2_in[0] = x2[n0]; y2_in[0] = y2[n0]; z2_in[0] = z2[n0];
      x2_in[1] = x2[n1]; y2_in[1] = y2[n1]; z2_in[1] = z2[n1];
      x2_in[2] = x2[n2]; y2_in[2] = y2[n2]; z2_in[2] = z2[n2];
      x2_in[3] = x2[n3]; y2_in[3] = y2[n3]; z2_in[3] = z2[n3];

      if (  (n_out = clip_2dx2d_great_circle( x1_in, y1_in, z1_in, n1_in, x2_in, y2_in, z2_in, n2_in,
					      x_out, y_out, z_out)) > 0) {
	xarea = great_circle_area ( n_out, x_out, y_out, z_out ) * mask_in[j1*nx1+i1];
	min_area = min(area1[j1*nx1+i1], area2[j2*nx2+i2]);
	if( xarea/min_area > AREA_RATIO_THRESH ) {
#ifdef debug_test_create_xgrid	  
	  printf("(i2,j2)=(%d,%d), (i1,j1)=(%d,%d), xarea=%g\n", i2, j2, i1, j1, xarea);
#endif
	  xgrid_area[nxgrid] = xarea;
	  xgrid_clon[nxgrid] = 0; /*z1l: will be developed very soon */
	  xgrid_clat[nxgrid] = 0; 
	  i_in[nxgrid]       = i1;
	  j_in[nxgrid]       = j1;
	  i_out[nxgrid]      = i2;
	  j_out[nxgrid]      = j2;
	  ++nxgrid;
	  if(nxgrid > MAXXGRID) error_handler("nxgrid is greater than MAXXGRID, increase MAXXGRID");
	}
      }
    }
  }

  
  free(area1);
  free(area2);  

  free(x1);
  free(y1);
  free(z1);
  free(x2);
  free(y2);
  free(z2);
  
  return nxgrid;
  
};/* create_xgrid_great_circle */

/*******************************************************************************
   Revise Sutherland-Hodgeman algorithm to find the vertices of the overlapping
   between any two grid boxes. It return the number of vertices for the exchange grid.
   Each edge of grid box is a part of great circle. All the points are cartesian
   coordinates. Here we are assuming each polygon is convex.
   RANGE_CHECK_CRITERIA is used to determine if the two grid boxes are possible to be
   overlap. The size should be between 0 and 0.5. The larger the range_check_criteria,
   the more expensive of the computatioin. When the value is close to 0,
   some small exchange grid might be lost. Suggest to use value 0.05 for C48. 
*******************************************************************************/

int clip_2dx2d_great_circle(const double x1_in[], const double y1_in[], const double z1_in[], int n1_in, 
			    const double x2_in[], const double y2_in[], const double z2_in [], int n2_in, 
			    double x_out[], double y_out[], double z_out[])
{
  struct Node *subjList=NULL;
  struct Node *clipList=NULL;
  struct Node *grid1List=NULL;
  struct Node *grid2List=NULL;
  struct Node *intersectList=NULL;
  struct Node *polyList=NULL;
  struct Node *curList=NULL;
  struct Node *firstIntersect=NULL, *curIntersect=NULL;
  struct Node *temp1=NULL, *temp2=NULL, *temp=NULL;
  
  int    i1, i2, i1p, i2p, i2p2, npts1, npts2;
  int    nintersect, n_out;
  int    maxiter1, maxiter2, iter1, iter2;
  int    found1, found2, curListNum;
  int    has_inbound, inbound;
  double pt1[MV][3], pt2[MV][3];
  double *p1_0=NULL, *p1_1=NULL;
  double *p2_0=NULL, *p2_1=NULL, *p2_2=NULL;
  double intersect[3];
  double u1, u2;
  double min_x1, max_x1, min_y1, max_y1, min_z1, max_z1;
  double min_x2, max_x2, min_y2, max_y2, min_z2, max_z2;
  static int first_call=1;

  
  /* first check the min and max of (x1_in, y1_in, z1_in) with (x2_in, y2_in, z2_in) */
  min_x1 = minval_double(n1_in, x1_in);
  max_x2 = maxval_double(n2_in, x2_in);  
  if(min_x1 >= max_x2+RANGE_CHECK_CRITERIA) return 0;
  max_x1 = maxval_double(n1_in, x1_in); 
  min_x2 = minval_double(n2_in, x2_in);
  if(min_x2 >= max_x1+RANGE_CHECK_CRITERIA) return 0;

  min_y1 = minval_double(n1_in, y1_in);
  max_y2 = maxval_double(n2_in, y2_in);  
  if(min_y1 >= max_y2+RANGE_CHECK_CRITERIA) return 0;
  max_y1 = maxval_double(n1_in, y1_in); 
  min_y2 = minval_double(n2_in, y2_in);
  if(min_y2 >= max_y1+RANGE_CHECK_CRITERIA) return 0;  

  min_z1 = minval_double(n1_in, z1_in);
  max_z2 = maxval_double(n2_in, z2_in);  
  if(min_z1 >= max_z2+RANGE_CHECK_CRITERIA) return 0;
  max_z1 = maxval_double(n1_in, z1_in); 
  min_z2 = minval_double(n2_in, z2_in);
  if(min_z2 >= max_z1+RANGE_CHECK_CRITERIA) return 0;    

  rewindList();

  grid1List = getNext();
  grid2List = getNext();  
  intersectList = getNext();
  polyList = getNext();
    
  /* insert points into SubjList and ClipList */
  for(i1=0; i1<n1_in; i1++) addEnd(grid1List, x1_in[i1], y1_in[i1], z1_in[i1], 0, 0, 0, -1);
  for(i2=0; i2<n2_in; i2++) addEnd(grid2List, x2_in[i2], y2_in[i2], z2_in[i2], 0, 0, 0, -1);
  npts1 = length(grid1List);
  npts2 = length(grid2List);

  n_out = 0;
  /* set the inside value */
#ifdef debug_test_create_xgrid
  printf("\nNOTE from clip_2dx2d_great_circle: begin to set inside value grid1List\n"); 
#endif  
  /* first check number of points in grid1 is inside grid2 */

  temp = grid1List;
  while(temp) {    
    if(insidePolygon(temp, grid2List))
      temp->isInside = 1;
    else
      temp->isInside = 0;
    temp = getNextNode(temp);
  }

#ifdef debug_test_create_xgrid
  printf("\nNOTE from clip_2dx2d_great_circle: begin to set inside value of grid2List\n"); 
#endif      
  /* check if grid2List is inside grid1List */
  temp = grid2List;

  while(temp) {
    if(insidePolygon(temp, grid1List))
      temp->isInside = 1;
    else
      temp->isInside = 0;
    temp = getNextNode(temp);
  }
      
  /* make sure the grid box is clockwise */
  
  /*make sure each polygon is convex, which is equivalent that the great_circle_area is positive */
  if( gridArea(grid1List) <= 0 )
    error_handler("create_xgrid.c(clip_2dx2d_great_circle): grid box 1 is not convex");
  if( gridArea(grid2List) <= 0 )
    error_handler("create_xgrid.c(clip_2dx2d_great_circle): grid box 2 is not convex");
  
#ifdef debug_test_create_xgrid
  printNode(grid1List, "grid1List");
  printNode(grid2List, "grid2List");
#endif  

  /* get the coordinates from grid1List and grid2List.
     Please not npts1 might not equal n1_in, npts2 might not equal n2_in because of pole
  */

  temp = grid1List;
  for(i1=0; i1<npts1; i1++) {
    getCoordinates(temp, pt1[i1]);
    temp = temp->Next;
  }
  temp = grid2List;
  for(i2=0; i2<npts2; i2++) {
    getCoordinates(temp, pt2[i2]);
    temp = temp->Next;
  }  
  
  firstIntersect=getNext();
  curIntersect = getNext();

#ifdef debug_test_create_xgrid  
  printf("\n\n************************ Start line_intersect_2D_3D ******************************\n");
#endif
  /* first find all the intersection points */
  nintersect = 0;
  for(i1=0; i1<npts1; i1++) {
    i1p = (i1+1)%npts1;
    p1_0 = pt1[i1];
    p1_1 = pt1[i1p];
    for(i2=0; i2<npts2; i2++) {
      i2p = (i2+1)%npts2;
      i2p2 = (i2+2)%npts2;
      p2_0 = pt2[i2];
      p2_1 = pt2[i2p];
      p2_2 = pt2[i2p2];
#ifdef debug_test_create_xgrid
      printf("\n******************************************************************************\n");
      printf(" i1 = %d, i2 = %d \n", i1, i2);
      printf("********************************************************************************\n");
#endif
      if( line_intersect_2D_3D(p1_0, p1_1, p2_0, p2_1, p2_2, intersect, &u1, &u2, &inbound) ) {
	int n_prev, n_cur;
	int is_in_subj, is_in_clip;

	/* from the value of u1, u2 and inbound, we can partially decide if a point is inside or outside of polygon */        
	
	/* add the intersection into intersetList, The intersection might already be in
	   intersectList and will be taken care addIntersect
	*/
	if(addIntersect(intersectList, intersect[0], intersect[1], intersect[2], 1, u1, u2, inbound, i1, i1p, i2, i2p)) {
	  /* add the intersection into the grid1List */
    
	  if(u1 == 1) {
	    insertIntersect(grid1List, intersect[0], intersect[1], intersect[2], 0.0, u2, inbound, p1_1[0], p1_1[1], p1_1[2]);
	  }
	  else
	    insertIntersect(grid1List, intersect[0], intersect[1], intersect[2], u1, u2, inbound, p1_0[0], p1_0[1], p1_0[2]);
	  /* when u1 == 0 or 1, need to adjust the vertice to intersect value for roundoff error */
	  if(u1==1) {
	    p1_1[0] = intersect[0];
	    p1_1[1] = intersect[1];
	    p1_1[2] = intersect[2];
	  }
	  else if(u1 == 0) {
	    p1_0[0] = intersect[0];
	    p1_0[1] = intersect[1];
	    p1_0[2] = intersect[2];
	  }
	  /* add the intersection into the grid2List */
	  if(u2==1)
	    insertIntersect(grid2List, intersect[0], intersect[1], intersect[2], 0.0, u1, 0, p2_1[0], p2_1[1], p2_1[2]);
	  else
	    insertIntersect(grid2List, intersect[0], intersect[1], intersect[2], u2, u1, 0, p2_0[0], p2_0[1], p2_0[2]);
	  /* when u2 == 0 or 1, need to adjust the vertice to intersect value for roundoff error */
	  if(u2==1) {
	    p2_1[0] = intersect[0];
	    p2_1[1] = intersect[1];
	    p2_1[2] = intersect[2];
	  } 
	  else if(u2 == 0) {
	    p2_0[0] = intersect[0];
	    p2_0[1] = intersect[1];
	    p2_0[2] = intersect[2];
	  }
	}
      }
    }
  }
  
  /* set inbound value for the points in intersectList that has inbound == 0,
     this will also set some inbound value of the points in grid1List
  */

  /* get the first point in intersectList has inbound = 2, if not, set inbound value */
  has_inbound = 0;
  /* loop through intersectList to see if there is any has inbound=1 or 2 */
  temp = intersectList;
  nintersect = length(intersectList);
  if(nintersect > 1) {
    getFirstInbound(intersectList, firstIntersect);
    if(firstIntersect->initialized) {
      has_inbound = 1;
    }
  }

  /* when has_inbound == 0, get the grid1List and grid2List */
  if( !has_inbound && nintersect > 1) {
    setInbound(intersectList, grid1List);
    getFirstInbound(intersectList, firstIntersect);
    if(firstIntersect->initialized) has_inbound = 1;
  }

  /* if has_inbound = 1, find the overlapping */
  n_out = 0;
  
  if(has_inbound) {
    maxiter1 = nintersect;
#ifdef debug_test_create_xgrid
    printf("\nNOTE from clip_2dx2d_great_circle: number of intersect is %d\n", nintersect);
    printf("\n size of grid2List is %d, size of grid1List is %d\n", length(grid2List), length(grid1List));
    printNode(intersectList, "beginning intersection list");
    printNode(grid2List, "beginning clip list");
    printNode(grid1List, "beginning subj list");
    printf("\n************************ End line_intersect_2D_3D **********************************\n\n");
#endif  
    temp1 = getNode(grid1List, *firstIntersect);
    if( temp1 == NULL) {
      double lon[10], lat[10];
      int i;
      xyz2latlon(n1_in, x1_in, y1_in, z1_in, lon, lat);
      for(i=0; i< n1_in; i++) printf("lon1 = %g, lat1 = %g\n", lon[i]*R2D, lat[i]*R2D);
      printf("\n");
      xyz2latlon(n2_in, x2_in, y2_in, z2_in, lon, lat);
      for(i=0; i< n2_in; i++) printf("lon2 = %g, lat2 = %g\n", lon[i]*R2D, lat[i]*R2D);
      printf("\n");	
	
      error_handler("firstIntersect is not in the grid1List");
    }
    addNode(polyList, *firstIntersect);
    nintersect--;
#ifdef debug_test_create_xgrid    
    printNode(polyList, "polyList at stage 1");
#endif      
	
    /* Loop over the grid1List and grid2List to find again the firstIntersect */
    curList = grid1List;
    curListNum = 0;
	
    /* Loop through curList to find the next intersection, the loop will end
       when come back to firstIntersect
    */
    copyNode(curIntersect, *firstIntersect);
    iter1 = 0;
    found1 = 0;

    while( iter1 < maxiter1 ) {
#ifdef debug_test_create_xgrid
      printf("\n----------- At iteration = %d\n\n", iter1+1 );
      printNode(curIntersect, "curIntersect at the begining of iter1");      
#endif
      /* find the curIntersect in curList and get the next intersection points */
      temp1 =  getNode(curList, *curIntersect);
      temp2 = temp1->Next;
      if( temp2 == NULL ) temp2 = curList;

      maxiter2 = length(curList);
      found2 = 0;
      iter2  = 0;
      /* Loop until find the next intersection */
      while( iter2 < maxiter2 ) {
	int temp2IsIntersect;

	temp2IsIntersect = 0;
	if( isIntersect( *temp2 ) ) { /* copy the point and switch to the grid2List */
	  struct Node *temp3;
	  
	  /* first check if temp2 is the firstIntersect */
	  if( sameNode( *temp2, *firstIntersect) ) {
	    found1 = 1; 
	    break;
	  }	  
	  
	  temp3 = temp2->Next;
	  if( temp3 == NULL) temp3 = curList;
	  if( temp3 == NULL) error_handler("creat_xgrid.c: temp3 can not be NULL");
	  found2 = 1;
	  /* if next node is inside or an intersection,
	     need to keep on curList
	  */
	  temp2IsIntersect = 1;
	  if( isIntersect(*temp3) || (temp3->isInside == 1)  ) found2 = 0;
	}
	if(found2) {
	  copyNode(curIntersect, *temp2);
	  break;
	}
	else {
	  addNode(polyList, *temp2);
#ifdef debug_test_create_xgrid    
	  printNode(polyList, "polyList at stage 2");
#endif      	  
	  if(temp2IsIntersect) { 
	    nintersect--;
	  }
	}
	temp2 = temp2->Next;
	if( temp2 == NULL ) temp2 = curList;
	iter2 ++;
      }
      if(found1) break;
    
      if( !found2 ) error_handler(" not found the next intersection ");

      /* if find the first intersection, the poly found */
      if( sameNode( *curIntersect, *firstIntersect) ) {
	found1 = 1; 
	break;
      }

      /* add curIntersect to polyList and remove it from intersectList and curList */
      addNode(polyList, *curIntersect);
#ifdef debug_test_create_xgrid    
      printNode(polyList, "polyList at stage 3");
#endif      
      nintersect--;
      
      
      /* switch curList */
      if( curListNum == 0) {
	curList = grid2List;
	curListNum = 1;
      }
      else {
	curList = grid1List;
	curListNum = 0;
      }
      iter1++;
    }
    if(!found1) error_handler("not return back to the first intersection");

    /* currently we are only clipping convex polygon to convex polygon */
    if( nintersect > 0) error_handler("After clipping, nintersect should be 0");

    /* copy the polygon to x_out, y_out, z_out */
    temp1 = polyList;
    while (temp1 != NULL) {
      getCoordinate(*temp1, x_out+n_out, y_out+n_out, z_out+n_out);
      temp1 = temp1->Next;
      n_out++;
    }
    
    /* if(n_out < 3) error_handler(" The clipped region has < 3 vertices"); */
    if( n_out < 3) n_out = 0;
#ifdef debug_test_create_xgrid
    printNode(polyList, "polyList after clipping");
#endif
  }

  /* check if grid1 is inside grid2 */
  if(n_out==0){
    /* first check number of points in grid1 is inside grid2 */
    int n, n1in2;
    /* One possible is that grid1List is inside grid2List */
#ifdef debug_test_create_xgrid
    printf("\nNOTE from clip_2dx2d_great_circle: check if grid1 is inside grid2\n");
#endif    
    n1in2 = 0;
    temp = grid1List;
    while(temp) {
      if(temp->intersect != 1) {
#ifdef debug_test_create_xgrid	
	printf("grid1->isInside = %d\n", temp->isInside);
#endif 
	if( temp->isInside == 1) n1in2++;
      }
      temp = getNextNode(temp);
    }
    if(npts1==n1in2) { /* grid1 is inside grid2 */
      n_out = npts1;
      n = 0;
      temp = grid1List;
      while( temp ) {
	getCoordinate(*temp, &x_out[n], &y_out[n], &z_out[n]);
	n++;
	temp = getNextNode(temp);
      }
    }
    if(n_out>0) return n_out;
  }
  
  /* check if grid2List is inside grid1List */
  if(n_out ==0){
    int n, n2in1;
#ifdef debug_test_create_xgrid
    printf("\nNOTE from clip_2dx2d_great_circle: check if grid2 is inside grid1\n");
#endif    
    
    temp = grid2List;
    n2in1 = 0;
    while(temp) {
      if(temp->intersect != 1) {
#ifdef debug_test_create_xgrid	
	printf("grid2->isInside = %d\n", temp->isInside);
#endif
	if( temp->isInside == 1) n2in1++;
      }
      temp = getNextNode(temp);
    }
      
    if(npts2==n2in1) { /* grid2 is inside grid1 */
      n_out = npts2;
      n = 0;
      temp = grid2List;
      while( temp ) {
	getCoordinate(*temp, &x_out[n], &y_out[n], &z_out[n]);
	n++;
	temp = getNextNode(temp);
      }
      
    }
  }  

  
  return n_out;
}


/* Intersects between the line a and the seqment s
   where both line and segment are great circle lines on the sphere represented by
   3D cartesian points.
   [sin sout] are the ends of a line segment
   returns true if the lines could be intersected, false otherwise.
   inbound means the direction of (a1,a2) go inside or outside of (q1,q2,q3)
*/

int line_intersect_2D_3D(double *a1, double *a2, double *q1, double *q2, double *q3,
			 double *intersect, double *u_a, double *u_q, int *inbound){

  /* Do this intersection by reprsenting the line a1 to a2 as a plane through the
     two line points and the origin of the sphere (0,0,0). This is the
     definition of a great circle arc.
  */
  double plane[9];
  double plane_p[2];
  double u;
  double p1[3], v1[3], v2[3];
  double c1[3], c2[3], c3[3];
  double coincident, sense, norm;
  int    i;
  int is_inter1, is_inter2;

  *inbound = 0;
  
  /* first check if any vertices are the same */
  if(samePoint(a1[0], a1[1], a1[2], q1[0], q1[1], q1[2])) {
    *u_a = 0;
    *u_q = 0;
    intersect[0] = a1[0];
    intersect[1] = a1[1];
    intersect[2] = a1[2];
#ifdef debug_test_create_xgrid
    printf("\nNOTE from line_intersect_2D_3D: u_a = %19.15f, u_q=%19.15f, inbound=%d\n", *u_a, *u_q, *inbound); 
#endif      
    return 1;
   }
   else if (samePoint(a1[0], a1[1], a1[2], q2[0], q2[1], q2[2])) {
    *u_a = 0;
    *u_q = 1;
    intersect[0] = a1[0];
    intersect[1] = a1[1];
    intersect[2] = a1[2];
#ifdef debug_test_create_xgrid
    printf("\nNOTE from line_intersect_2D_3D: u_a = %19.15f, u_q=%19.15f, inbound=%d\n", *u_a, *u_q, *inbound); 
#endif      
    return 1;
  }
   else if(samePoint(a2[0], a2[1], a2[2], q1[0], q1[1], q1[2])) {
#ifdef debug_test_create_xgrid
    printf("\nNOTE from line_intersect_2D_3D: u_a = %19.15f, u_q=%19.15f, inbound=%d\n", *u_a, *u_q, *inbound); 
#endif           
    *u_a = 1;
    *u_q = 0;
    intersect[0] = a2[0];
    intersect[1] = a2[1];
    intersect[2] = a2[2];
    return 1;
   }
   else if (samePoint(a2[0], a2[1], a2[2], q2[0], q2[1], q2[2])) {
#ifdef debug_test_create_xgrid
    printf("\nNOTE from line_intersect_2D_3D: u_a = %19.15f, u_q=%19.15f, inbound=%d\n", *u_a, *u_q, *inbound); 
#endif           
    *u_a = 1;
    *u_q = 1;
    intersect[0] = a2[0];
    intersect[1] = a2[1];
    intersect[2] = a2[2];
    return 1;
  }


  /* Load points defining plane into variable (these are supposed to be in counterclockwise order) */
  plane[0]=q1[0];
  plane[1]=q1[1];
  plane[2]=q1[2];
  plane[3]=q2[0];
  plane[4]=q2[1];
  plane[5]=q2[2];
  plane[6]=0.0;
  plane[7]=0.0;
  plane[8]=0.0;

  /* Intersect the segment with the plane */
  is_inter1 = intersect_tri_with_line(plane, a1, a2, plane_p, u_a);

  if(!is_inter1)
     return 0;

  if(fabs(*u_a) < EPSLN8) *u_a = 0;
  if(fabs(*u_a-1) < EPSLN8) *u_a = 1;

  
#ifdef debug_test_create_xgrid
    printf("\nNOTE from line_intersect_2D_3D: u_a = %19.15f\n", *u_a);
#endif


  if( (*u_a < 0) || (*u_a > 1) ) return 0;

  /* Load points defining plane into variable (these are supposed to be in counterclockwise order) */
  plane[0]=a1[0];
  plane[1]=a1[1];
  plane[2]=a1[2];
  plane[3]=a2[0];
  plane[4]=a2[1];
  plane[5]=a2[2];
  plane[6]=0.0;
  plane[7]=0.0;
  plane[8]=0.0;

  /* Intersect the segment with the plane */
  is_inter2 = intersect_tri_with_line(plane, q1, q2, plane_p, u_q);

  if(!is_inter2)
     return 0;

  if(fabs(*u_q) < EPSLN8) *u_q = 0;
  if(fabs(*u_q-1) < EPSLN8) *u_q = 1;
#ifdef debug_test_create_xgrid
    printf("\nNOTE from line_intersect_2D_3D: u_q = %19.15f\n", *u_q);
#endif
  

  if( (*u_q < 0) || (*u_q > 1) ) return 0;
  
  u =*u_a;
  
  /* The two planes are coincidental */
  vect_cross(a1, a2, c1);
  vect_cross(q1, q2, c2);
  vect_cross(c1, c2, c3);
  coincident = metric(c3);

  if(fabs(coincident) < EPSLN30) return 0;
  
  /* Calculate point of intersection */
  intersect[0]=a1[0] + u*(a2[0]-a1[0]);
  intersect[1]=a1[1] + u*(a2[1]-a1[1]);
  intersect[2]=a1[2] + u*(a2[2]-a1[2]);

  norm = metric( intersect );
  for(i = 0; i < 3; i ++) intersect[i] /= norm;

  /* when u_q =0 or u_q =1, the following could not decide the inbound value */
  if(*u_q != 0 && *u_q != 1){
  
    p1[0] = a2[0]-a1[0];
    p1[1] = a2[1]-a1[1];
    p1[2] = a2[2]-a1[2];
    v1[0] = q2[0]-q1[0];
    v1[1] = q2[1]-q1[1];
    v1[2] = q2[2]-q1[2];
    v2[0] = q3[0]-q2[0];
    v2[1] = q3[1]-q2[1];
    v2[2] = q3[2]-q2[2];

    vect_cross(v1, v2, c1);
    vect_cross(v1, p1, c2);

    sense = dot(c1, c2);
    *inbound = 1;
    if(sense > 0) *inbound = 2; /* v1 going into v2 in CCW sense */
  }
#ifdef debug_test_create_xgrid
    printf("\nNOTE from line_intersect_2D_3D: inbound=%d\n", *inbound); 
#endif      
  
  return 1;
}


/*------------------------------------------------------------------------------
  double poly_ctrlat(const double x[], const double y[], int n)
  This routine is used to calculate the latitude of the centroid 
   ---------------------------------------------------------------------------*/

double poly_ctrlat(const double x[], const double y[], int n)
{
  double ctrlat = 0.0;
  int    i;

  for (i=0;i<n;i++) {
    int ip = (i+1) % n;
    double dx = (x[ip]-x[i]);
    double dy, avg_y, hdy;
    double lat1, lat2;
    lat1 = y[ip];
    lat2 = y[i];
    dy = lat2 - lat1;
    hdy = dy*0.5;
    avg_y = (lat1+lat2)*0.5;
    if      (dx==0.0) continue;
    if(dx > M_PI)  dx = dx - 2.0*M_PI;
    if(dx < -M_PI) dx = dx + 2.0*M_PI;

    if ( fabs(hdy)< SMALL_VALUE ) /* cheap area calculation along latitude */
      ctrlat -= dx*(2*cos(avg_y) + lat2*sin(avg_y) - cos(lat1) );
    else 
      ctrlat -= dx*( (sin(hdy)/hdy)*(2*cos(avg_y) + lat2*sin(avg_y)) - cos(lat1) );
  }
  return (ctrlat*RADIUS*RADIUS);
}; /* poly_ctrlat */        

/*------------------------------------------------------------------------------
  double poly_ctrlon(const double x[], const double y[], int n, double clon)
  This routine is used to calculate the lontitude of the centroid.
   ---------------------------------------------------------------------------*/
double poly_ctrlon(const double x[], const double y[], int n, double clon)
{
  double ctrlon = 0.0;
  int    i;

  clon = clon;
  for (i=0;i<n;i++) {
    int ip = (i+1) % n;
    double phi1, phi2, dphi, lat1, lat2, dphi1, dphi2;
    double f1, f2, fac, fint;
    phi1   = x[ip];
    phi2   = x[i];
    lat1 = y[ip];
    lat2 = y[i];    
    dphi   = phi1 - phi2;
    
    if      (dphi==0.0) continue;

    f1 = 0.5*(cos(lat1)*sin(lat1)+lat1);
    f2 = 0.5*(cos(lat2)*sin(lat2)+lat2);

     /* this will make sure longitude of centroid is at 
        the same interval as the center of any grid */  
    if(dphi > M_PI)  dphi = dphi - 2.0*M_PI;
    if(dphi < -M_PI) dphi = dphi + 2.0*M_PI;
    dphi1 = phi1 - clon;
    if( dphi1 > M_PI) dphi1 -= 2.0*M_PI;
    if( dphi1 <-M_PI) dphi1 += 2.0*M_PI;
    dphi2 = phi2 -clon;
    if( dphi2 > M_PI) dphi2 -= 2.0*M_PI;
    if( dphi2 <-M_PI) dphi2 += 2.0*M_PI;    

    if(abs(dphi2 -dphi1) < M_PI) {
      ctrlon -= dphi * (dphi1*f1+dphi2*f2)/2.0;
    }
    else {
      if(dphi1 > 0.0)
	fac = M_PI;
      else
	fac = -M_PI;
      fint = f1 + (f2-f1)*(fac-dphi1)/abs(dphi);
      ctrlon -= 0.5*dphi1*(dphi1-fac)*f1 - 0.5*dphi2*(dphi2+fac)*f2
	+ 0.5*fac*(dphi1+dphi2)*fint;
	}
    
  }
  return (ctrlon*RADIUS*RADIUS);
};   /* poly_ctrlon */

/* -----------------------------------------------------------------------------
   double box_ctrlat(double ll_lon, double ll_lat, double ur_lon, double ur_lat)
   This routine is used to calculate the latitude of the centroid.
   ---------------------------------------------------------------------------*/
double box_ctrlat(double ll_lon, double ll_lat, double ur_lon, double ur_lat)
{
  double dphi = ur_lon-ll_lon;
  double ctrlat;
  
  if(dphi > M_PI)  dphi = dphi - 2.0*M_PI;
  if(dphi < -M_PI) dphi = dphi + 2.0*M_PI;
  ctrlat = dphi*(cos(ur_lat) + ur_lat*sin(ur_lat)-(cos(ll_lat) + ll_lat*sin(ll_lat)));
  return (ctrlat*RADIUS*RADIUS); 
}; /* box_ctrlat */

/*------------------------------------------------------------------------------
  double box_ctrlon(double ll_lon, double ll_lat, double ur_lon, double ur_lat, double clon)
  This routine is used to calculate the lontitude of the centroid 
   ----------------------------------------------------------------------------*/
double box_ctrlon(double ll_lon, double ll_lat, double ur_lon, double ur_lat, double clon)
{
  double phi1, phi2, dphi, lat1, lat2, dphi1, dphi2;
  double f1, f2, fac, fint;  
  double ctrlon  = 0.0;
  int i;
  clon = clon;  
  for( i =0; i<2; i++) {
    if(i == 0) {
      phi1 = ur_lon;
      phi2 = ll_lon;
      lat1 = lat2 = ll_lat;
    }
    else {
      phi1 = ll_lon;
      phi2 = ur_lon;
      lat1 = lat2 = ur_lat;
    }
    dphi   = phi1 - phi2;
    f1 = 0.5*(cos(lat1)*sin(lat1)+lat1);
    f2 = 0.5*(cos(lat2)*sin(lat2)+lat2);

    if(dphi > M_PI)  dphi = dphi - 2.0*M_PI;
    if(dphi < -M_PI) dphi = dphi + 2.0*M_PI;
    /* make sure the center is in the same grid box. */
    dphi1 = phi1 - clon;
    if( dphi1 > M_PI) dphi1 -= 2.0*M_PI;
    if( dphi1 <-M_PI) dphi1 += 2.0*M_PI;
    dphi2 = phi2 -clon;
    if( dphi2 > M_PI) dphi2 -= 2.0*M_PI;
    if( dphi2 <-M_PI) dphi2 += 2.0*M_PI;    

    if(abs(dphi2 -dphi1) < M_PI) {
      ctrlon -= dphi * (dphi1*f1+dphi2*f2)/2.0;
    }
    else {
      if(dphi1 > 0.0)
	fac = M_PI;
      else
	fac = -M_PI;
      fint = f1 + (f2-f1)*(fac-dphi1)/abs(dphi);
      ctrlon -= 0.5*dphi1*(dphi1-fac)*f1 - 0.5*dphi2*(dphi2+fac)*f2
	+ 0.5*fac*(dphi1+dphi2)*fint;
    }
  }
  return (ctrlon*RADIUS*RADIUS);    
} /* box_ctrlon */

/*******************************************************************************
  double grid_box_radius(double *x, double *y, double *z, int n);
  Find the radius of the grid box, the radius is defined the
  maximum distance between any two vertices
*******************************************************************************/ 
double grid_box_radius(const double *x, const double *y, const double *z, int n)
{
  double radius;
  int i, j;
  
  radius = 0;
  for(i=0; i<n-1; i++) {
    for(j=i+1; j<n; j++) {
      radius = max(radius, pow(x[i]-x[j],2.)+pow(y[i]-y[j],2.)+pow(z[i]-z[j],2.));
    }
  }
  
  radius = sqrt(radius);

  return (radius);
  
}; /* grid_box_radius */

/*******************************************************************************
  double dist_between_boxes(const double *x1, const double *y1, const double *z1, int n1,
			    const double *x2, const double *y2, const double *z2, int n2);
  Find the distance between any two grid boxes. The distance is defined by the maximum
  distance between any vertices of these two box
*******************************************************************************/
double dist_between_boxes(const double *x1, const double *y1, const double *z1, int n1,
			  const double *x2, const double *y2, const double *z2, int n2)
{
  double dist;
  int i, j;

  dist = 0.0;
  for(i=0; i<n1; i++) {
    for(j=0; j<n2; j++) {   
      dist = max(dist, pow(x1[i]-x2[j],2.)+pow(y1[i]-y2[j],2.)+pow(z1[i]-z2[j],2.));
    }
  }

  dist = sqrt(dist);
  return (dist);

}; /* dist_between_boxes */

/*******************************************************************************
 int inside_edge(double x0, double y0, double x1, double y1, double x, double y)
 determine a point(x,y) is inside or outside a given edge with vertex,
 (x0,y0) and (x1,y1). return 1 if inside and 0 if outside. <y1-y0, -(x1-x0)> is
 the outward edge normal from vertex <x0,y0> to <x1,y1>. <x-x0,y-y0> is the vector
 from <x0,y0> to <x,y>. 
 if Inner produce <x-x0,y-y0>*<y1-y0, -(x1-x0)> > 0, outside, otherwise inside.
 inner product value = 0 also treate as inside.
*******************************************************************************/
int inside_edge(double x0, double y0, double x1, double y1, double x, double y)
{
   const double SMALL = 1.e-12;
   double product;
   
   product = ( x-x0 )*(y1-y0) + (x0-x1)*(y-y0);
   return (product<=SMALL) ? 1:0;
   
 }; /* inside_edge */


/* The following is a test program to test subroutines in create_xgrid.c */

#ifdef test_create_xgrid

#include "create_xgrid.h"
#include <math.h>

#define D2R (M_PI/180)
#define R2D (180/M_PI)
#define MAXPOINT 1000

int main(int argc, char* argv[])
{

  double lon1_in[MAXPOINT], lat1_in[MAXPOINT];
  double lon2_in[MAXPOINT], lat2_in[MAXPOINT];
  double x1_in[MAXPOINT], y1_in[MAXPOINT], z1_in[MAXPOINT];
  double x2_in[MAXPOINT], y2_in[MAXPOINT], z2_in[MAXPOINT];
  double lon_out[20], lat_out[20];
  double x_out[20], y_out[20], z_out[20];
  int    n1_in, n2_in, n_out, i, j;
  int    nlon1=0, nlat1=0, nlon2=0, nlat2=0;
  int    n;
  int    ntest = 11;
  

  for(n=11; n<=ntest; n++) {
    
    switch (n) {
    case 1:
      /****************************************************************

       test clip_2dx2d_great_cirle case 1:
       box 1: (20,10), (20,12), (22,12), (22,10)
       box 2: (21,11), (21,14), (24,14), (24,11)
       out  : (21, 12.0018), (22, 12), (22, 11.0033), (21, 11) 

      ****************************************************************/
      n1_in = 4; n2_in = 4;  
      /* first a simple lat-lon grid box to clip another lat-lon grid box */
      lon1_in[0] = 20; lat1_in[0] = 10;
      lon1_in[1] = 20; lat1_in[1] = 12;
      lon1_in[2] = 22; lat1_in[2] = 12;
      lon1_in[3] = 22; lat1_in[3] = 10;
      lon2_in[0] = 21; lat2_in[0] = 11;
      lon2_in[1] = 21; lat2_in[1] = 14;
      lon2_in[2] = 24; lat2_in[2] = 14;
      lon2_in[3] = 24; lat2_in[3] = 11;
      break;
      
    case 2:
      /****************************************************************

        test clip_2dx2d_great_cirle case 2: two identical box
        box 1: (20,10), (20,12), (22,12), (22,10)
        box 2: (20,10), (20,12), (22,12), (22,10)
        out  : (20,10), (20,12), (22,12), (22,10)

      ****************************************************************/      
      lon1_in[0] = 20; lat1_in[0] = 10;
      lon1_in[1] = 20; lat1_in[1] = 12;
      lon1_in[2] = 22; lat1_in[2] = 12;
      lon1_in[3] = 22; lat1_in[3] = 10;  

      for(i=0; i<n2_in; i++) {
	lon2_in[i] = lon1_in[i];
	lat2_in[i] = lat1_in[i];
      }
      break;

    case 3:
      /****************************************************************

       test clip_2dx2d_great_cirle case 3: one cubic sphere grid close to the pole with lat-lon grid.
       box 1: (251.7, 88.98), (148.3, 88.98), (57.81, 88.72), (342.2, 88.72) 
       box 2: (150, 88), (150, 90), (152.5, 90), (152.5, 88)
       out  : (152.5, 89.0642), (150, 89.0165), (0, 90) 

      ****************************************************************/
      n1_in = 4; n2_in = 4;  
      /* first a simple lat-lon grid box to clip another lat-lon grid box */
      lon1_in[0] = 251.7; lat1_in[0] = 88.98;
      lon1_in[1] = 148.3; lat1_in[1] = 88.98;
      lon1_in[2] = 57.81; lat1_in[2] = 88.72;
      lon1_in[3] = 342.2; lat1_in[3] = 88.72;

      lon2_in[0] = 150; lat2_in[0] = 88;
      lon2_in[1] = 150; lat2_in[1] = 90;
      lon2_in[2] = 152.5; lat2_in[2] = 90;
      lon2_in[3] = 152.5; lat2_in[3] = 88;
      /*
      for(i=0; i<4; i++) {
	lon2_in[i] = lon1_in[i];
	lat2_in[i] = lat1_in[i];
      }
      */
      break;

    case 4:
      /****************************************************************

       test clip_2dx2d_great_cirle case 4: One box contains the pole
       box 1: (-160, 88.5354), (152.011, 87.8123) , (102.985, 88.4008), (20, 89.8047) 
       box 2: (145,88), (145,90), (150,90), (150,88)
       out  : (145.916, 88.0011), (145, 88.0249), (0, 90), (150, 88) 

      ****************************************************************/
      n1_in = 4; n2_in = 4;  
      /* first a simple lat-lon grid box to clip another lat-lon grid box */

      lon1_in[0] = -160;  lat1_in[0] = 88.5354;
      lon1_in[1] = 152.011; lat1_in[1] = 87.8123;
      lon1_in[2] = 102.985; lat1_in[2] = 88.4008;
      lon1_in[3] = 20; lat1_in[3] = 89.8047;  

      lon2_in[0] = 145; lat2_in[0] = 88;
      lon2_in[1] = 145; lat2_in[1] = 90;
      lon2_in[2] = 150; lat2_in[2] = 90;
      lon2_in[3] = 150; lat2_in[3] = 88;
      break;

    case 5:
      /****************************************************************

       test clip_2dx2d_great_cirle case 5: One tripolar grid around the pole with lat-lon grid.
       box 1: (-202.6, 87.95), (-280, 89.56), (-100, 90), (-190, 88)
       box 2: (21,11), (21,14), (24,14), (24,11)
       out  : (150, 88.7006), (145,  88.9507), (0, 90) 

      ****************************************************************/
      n1_in = 4; n2_in = 4;  
      /* first a simple lat-lon grid box to clip another lat-lon grid box */

      lon1_in[0] = -202.6;  lat1_in[0] = 87.95;
      lon1_in[1] = -280.;   lat1_in[1] = 89.56;
      lon1_in[2] = -100.0; lat1_in[2] = 90;
      lon1_in[3] = -190.; lat1_in[3] = 88;  

      lon2_in[0] = 145; lat2_in[0] = 88;
      lon2_in[1] = 145; lat2_in[1] = 90;
      lon2_in[2] = 150; lat2_in[2] = 90;
      lon2_in[3] = 150; lat2_in[3] = 88;
      break; 

    case 6:
      /****************************************************************

       test clip_2dx2d_great_cirle case 6: One cubic sphere grid arounc the pole with one tripolar grid box
                                       around the pole.
       box 1: (-160, 88.5354), (152.011, 87.8123) , (102.985, 88.4008), (20, 89.8047) 
       box 2: (-202.6, 87.95), (-280, 89.56), (-100, 90), (-190, 88)
       out  : (170, 88.309), (157.082, 88.0005), (83.714, 89.559), (80, 89.6094), (0, 90), (200, 88.5354) 


      ****************************************************************/
      n1_in = 4; n2_in = 4;  
      /* first a simple lat-lon grid box to clip another lat-lon grid box */

      lon1_in[0] = -160;  lat1_in[0] = 88.5354;
      lon1_in[1] = 152.011; lat1_in[1] = 87.8123;
      lon1_in[2] = 102.985; lat1_in[2] = 88.4008;
      lon1_in[3] = 20; lat1_in[3] = 89.8047;  

      lon2_in[0] = -202.6;  lat2_in[0] = 87.95;
      lon2_in[1] = -280.;   lat2_in[1] = 89.56;
      lon2_in[2] = -100.0;  lat2_in[2] = 90;
      lon2_in[3] = -190.;   lat2_in[3] = 88;
      break;

    case 7:
      /****************************************************************

       test clip_2dx2d_great_cirle case 7: One small grid box inside a big grid box.
       box 1: (20,10), (20,12), (22,12), (22,10)
       box 2: (18,8), (18,14), (24,14), (24,8)
       out  : (20,10), (20,12), (22,12), (22,10)

      ****************************************************************/
      n1_in = 4; n2_in = 4;  
      /* first a simple lat-lon grid box to clip another lat-lon grid box */
      lon1_in[0] = 20; lat1_in[0] = 10;
      lon1_in[1] = 20; lat1_in[1] = 12;
      lon1_in[2] = 22; lat1_in[2] = 12;
      lon1_in[3] = 22; lat1_in[3] = 10;
      lon2_in[0] = 18; lat2_in[0] = 8;
      lon2_in[1] = 18; lat2_in[1] = 14;
      lon2_in[2] = 24; lat2_in[2] = 14;
      lon2_in[3] = 24; lat2_in[3] = 8;
      break;

    case 8:
      /****************************************************************

       test clip_2dx2d_great_cirle case 8: Cubic sphere grid at tile = 1, point (i=25,j=1)
          with N45 at (i=141,j=23)
       box 1: 
       box 2: 
       out  : None

      ****************************************************************/
      n1_in = 4; n2_in = 4;  
      /* first a simple lat-lo
	 n grid box to clip another lat-lon grid box */
      lon1_in[0] = 350.0; lat1_in[0] = -45;
      lon1_in[1] = 350.0; lat1_in[1] = -43.43;
      lon1_in[2] = 352.1; lat1_in[2] = -43.41;
      lon1_in[3] = 352.1; lat1_in[3] = -44.98;
      lon2_in[0] = 350.0;   lat2_in[0] = -46;
      lon2_in[1] = 350.0;   lat2_in[1] = -44;
      lon2_in[2] = 352.5; lat2_in[2] = -44;
      lon2_in[3] = 352.5; lat2_in[3] = -46;
      break;      

    case 9:      
      /****************************************************************

       test clip_2dx2d_great_cirle case 9: Cubic sphere grid at tile = 1, point (i=1,j=1)
          with N45 at (i=51,j=61)
       box 1: 
       box 2: 
       out  : None

      ****************************************************************/
      n1_in = 4; n2_in = 4;  

      lon1_in[0] = 305.0; lat1_in[0] = -35.26;
      lon1_in[1] = 305.0; lat1_in[1] = -33.80;
      lon1_in[2] = 306.6; lat1_in[2] = -34.51;
      lon1_in[3] = 306.6; lat1_in[3] = -35.99;
      lon2_in[0] = 125;   lat2_in[0] = 32;
      lon2_in[1] = 125;   lat2_in[1] = 34;
      lon2_in[2] = 127.5; lat2_in[2] = 34;
      lon2_in[3] = 127.5; lat2_in[3] = 32;
      break;      

    case 10:      
      /****************************************************************

       test clip_2dx2d_great_cirle case 10: Cubic sphere grid at tile = 3, point (i=24,j=1)
          with N45 at (i=51,j=46)
       box 1: 
       box 2: 
       out  : None

      ****************************************************************/
      n1_in = 4; n2_in = 4;  

      lon1_in[0] = 125.0; lat1_in[0] = 1.46935;
      lon1_in[1] = 126.573; lat1_in[1] = 1.5091;
      lon1_in[2] = 126.573; lat1_in[2] = 0;
      lon1_in[3] = 125.0; lat1_in[3] = 0;
      lon2_in[0] = 125;   lat2_in[0] = 0;
      lon2_in[1] = 125;   lat2_in[1] = 2;
      lon2_in[2] = 127.5; lat2_in[2] = 2;
      lon2_in[3] = 127.5; lat2_in[3] = 0;
      break;      

    case 11:      
      /****************************************************************

       test clip_2dx2d_great_cirle case 10: Cubic sphere grid at tile = 3, point (i=24,j=1)
          with N45 at (i=51,j=46)
       box 1: 
       box 2: 
       out  : 

      ****************************************************************/
      nlon1 = 1;
      nlat1 = 1;
      nlon2 = 1;
      nlat2 = 1;
      n1_in = (nlon1+1)*(nlat1+1);
      n2_in = (nlon2+1)*(nlat2+1);

      lon1_in[0] = 350.0; lat1_in[0] = 90.00;
      lon1_in[1] = 170.0; lat1_in[1] = 87.92;
      lon1_in[2] = 260.0; lat1_in[2] = 87.92;
      lon1_in[3] = 215.0;  lat1_in[3] = 87.06;
      
/*       lon1_in[0] = 35.0; lat1_in[0] = 87.06; */
/*       lon1_in[1] = 80.0; lat1_in[1] = 87.92; */
/*       lon1_in[2] = 125.0; lat1_in[2] = 87.06; */
/*       lon1_in[3] = 350.0; lat1_in[3] = 87.92; */
/*       lon1_in[4] = 350.0; lat1_in[4] = 90.00; */
/*       lon1_in[5] = 170.0; lat1_in[5] = 87.92; */
/*       lon1_in[6] = 305.0; lat1_in[6] = 87.06; */
/*       lon1_in[7] = 260.0; lat1_in[7] = 87.92; */
/*       lon1_in[8] = 215.0;  lat1_in[8] = 87.06; */

      lon2_in[0] = 167.5; lat2_in[0] = 88;
      lon2_in[1] = 170;   lat2_in[1] = 88;
      lon2_in[2] = 167.5; lat2_in[2] = 90;
      lon2_in[3] = 170;   lat2_in[3] = 90;
      
/*       nlon1 = 3; */
/*       nlat1 = 2; */
/*       nlon2 = 1; */
/*       nlat2 = 1; */
/*       n1_in = (nlon1+1)*(nlat1+1); */
/*       n2_in = (nlon2+1)*(nlat2+1); */

/*       lon1_in[0] = 35.00;     lat1_in[0] = -59.90; */
/*       lon1_in[1] = 37.64;     lat1_in[1] = -58.69; */
/*       lon1_in[2] = 40.07;     lat1_in[2] = -57.44; */
/*       lon1_in[3] = 42.32;     lat1_in[3] = -56.15; */
/*       lon1_in[4] = 32.36;     lat1_in[4] = -58.69; */
/*       lon1_in[5] = 35.00;     lat1_in[5] = -57.56; */
/*       lon1_in[6] = 37.45;     lat1_in[6] = -56.39; */
/*       lon1_in[7] = 39.74;     lat1_in[7] = -55.18; */
/*       lon1_in[8] = 29.93;     lat1_in[8] = -57.44; */
/*       lon1_in[9] = 32.55;     lat1_in[9] = -56.39; */
/*       lon1_in[10] = 35.00;     lat1_in[10] = -55.29; */
/*       lon1_in[11] = 37.30;     lat1_in[11] = -54.16; */
/*       lon2_in[0] = 35;   lat2_in[0] = -58; */
/*       lon2_in[1] = 37.5; lat2_in[1] = -58; */
/*       lon2_in[2] = 35;   lat2_in[2] = -56; */
/*       lon2_in[3] = 37.5; lat2_in[3] = -56; */

/*       nlon1 = 1; */
/*       nlat1 = 1; */
/*       nlon2 = 1; */
/*       nlat2 = 1; */
/*       n1_in = (nlon1+1)*(nlat1+1); */
/*       n2_in = (nlon2+1)*(nlat2+1); */

/*       lon1_in[0] = 305;     lat1_in[0] = -35.26; */
/*       lon1_in[1] = 306;     lat1_in[1] = -35.99; */
/*       lon1_in[2] = 305;     lat1_in[2] = -33.80; */
/*       lon1_in[3] = 306;     lat1_in[3] = -34.51; */
/*       lon2_in[0] = 305;   lat2_in[0] = -34; */
/*       lon2_in[1] = 307.5; lat2_in[1] = -34; */
/*       lon2_in[2] = 305;   lat2_in[2] = -32; */
/*       lon2_in[3] = 307.5; lat2_in[3] = -32; */

       nlon1 = 2; 
       nlat1 = 2; 
       nlon2 = 1; 
       nlat2 = 1;
      n1_in = (nlon1+1)*(nlat1+1);
      n2_in = (nlon2+1)*(nlat2+1);
       
      lon1_in[0] = 111.3; lat1_in[0] = 1.591;
      lon1_in[1] = 109.7; lat1_in[1] = 2.926;
      lon1_in[2] = 108.2; lat1_in[2] = 4.256;
      lon1_in[3] = 110.0; lat1_in[3] = 0.000;
      lon1_in[4] = 108.4; lat1_in[4] = 1.335;
      lon1_in[5] = 106.8; lat1_in[5] = 2.668;
      lon1_in[6] = 108.7; lat1_in[6] = -1.591;
      lon1_in[7] = 107.1; lat1_in[7] = -0.256;
      lon1_in[8] = 105.5;  lat1_in[8] = 1.078;

      lon2_in[0] = 107.5; lat2_in[0] = 0;
      lon2_in[1] = 110;   lat2_in[1] = 0;
      lon2_in[2] = 107.5; lat2_in[2] = 2;
      lon2_in[3] = 110;   lat2_in[3] = 2;
      
      break;      
      
    case 12:
      /****************************************************************

       test : create_xgrid_great_circle
       box 1: (20,10), (20,12), (22,12), (22,10)
       box 2: (21,11), (21,14), (24,14), (24,11)
       out  : (21, 12.0018), (22, 12), (22, 11.0033), (21, 11) 

      ****************************************************************/
      nlon1 = 2;
      nlat1 = 2;
      nlon2 = 3;
      nlat2 = 3;
      n1_in = (nlon1+1)*(nlat1+1);
      n2_in = (nlon2+1)*(nlat2+1);
      
      /* first a simple lat-lon grid box to clip another lat-lon grid box */
      for(j=0; j<=nlat1; j++) for(i=0; i<=nlon1; i++){
	lon1_in[j*(nlon1+1)+i] = 20.0 + (i-1)*2.0;
	lat1_in[j*(nlon1+1)+i] = 10.0 + (j-1)*2.0;
      }
       for(j=0; j<=nlat2; j++) for(i=0; i<=nlon2; i++){
	lon2_in[j*(nlon2+1)+i] = 19.0 + (i-1)*2.0;
	lat2_in[j*(nlon2+1)+i] = 9.0 + (j-1)*2.0;
      }
	
      break;

    case 13:

      nlon1 = 1;
      nlat1 = 1;
      nlon2 = 1;
      nlat2 = 1;
      n1_in = (nlon1+1)*(nlat1+1);
      n2_in = (nlon2+1)*(nlat2+1);

/*       lon1_in[0] = ; lat1_in[0] = ; */
/*       lon1_in[1] = ; lat1_in[1] = ; */
/*       lon1_in[2] = ; lat1_in[2] = ; */
/*       lon1_in[3] = ; lat1_in[3] = ; */
/*       lon2_in[0] = ; lat2_in[0] = ; */
/*       lon2_in[1] = ; lat2_in[1] = ; */
/*       lon2_in[2] = ; lat2_in[2] = ; */
/*       lon2_in[3] = ; lat2_in[3] = ;     */
      
/*       lon1_in[0] = 1.35536; lat1_in[0] = 1.16251; */
/*       lon1_in[1] = 1.36805; lat1_in[1] = 1.15369; */
/*       lon1_in[2] = 1.37843; lat1_in[2] = 1.16729; */
/*       lon1_in[3] = 1.39048; lat1_in[3] = 1.15826; */
/*       lon2_in[0] = 1.34611; lat2_in[0] = 1.16372; */
/*       lon2_in[1] = 1.35616; lat2_in[1] = 1.15802;    */
/*       lon2_in[2] = 1.35143; lat2_in[2] = 1.16509; */
/*       lon2_in[3] = 1.36042; lat2_in[3] = 1.15913; */

/*       lon1_in[0] = 12.508065121288551; lat1_in[0] = -87.445883646793547; */
/*       lon1_in[1] = 325.425637772; lat1_in[1] = -86.481216821859505; */
/*       lon1_in[2] = 97.5; lat1_in[2] = -89.802136057677174; */
/*       lon1_in[3] = 277.5; lat1_in[3] = -87.615232005344637; */

/*       for(j=0; j<=nlat2; j++) for(i=0; i<=nlon2; i++) { */
/* 	lon2_in[j*(nlon2+1)+i] = -280.0 + i*1.0; */
/* 	lat2_in[j*(nlon2+1)+i] = -90.0 + j*8.0; */
/*       } */
      lon1_in[0] = 120.369397984526174; lat1_in[0] = 16.751543427495864;
      lon1_in[1] = 119.999999999999986; lat1_in[1] = 16.751871929590038;
      lon1_in[2] = 120.369397846883501; lat1_in[2] = 16.397797979598028;
      lon1_in[3] = 119.999999999999986; lat1_in[3] = 16.398120477217255;
      lon2_in[0] = 120.369415056522087; lat2_in[0] = 16.752176828509153;
      lon2_in[1] = 119.999999999999986; lat2_in[1] = 16.752505523196167;
      lon2_in[2] = 120.369415056522087; lat2_in[2] = 16.397797949548146;
      lon2_in[3] = 119.999999999999986; lat2_in[3] = 16.398120477217255;      
      
      break;

      
    default:
      error_handler("test_create_xgrid: incorrect case number");
    }

    /* convert to radian */

    for(i=0; i<n1_in; i++) {
      lon1_in[i] *= D2R; lat1_in[i] *=D2R;
    }
    for(i=0; i<n2_in; i++) {
      lon2_in[i] *= D2R; lat2_in[i] *=D2R;
    }
    
  
    printf("\n*********************************************************\n");
    printf("\n               Case %d                                    \n", n);
    printf("\n*********************************************************\n");


    if( n > 10 ) {
      int nxgrid;
      int *i1, *j1, *i2, *j2;
      double *xarea, *xclon, *xclat, *mask1;

      mask1 = (double *)malloc(nlon1*nlat1*sizeof(double));
      i1    = (int    *)malloc(MAXXGRID*sizeof(int));
      j1    = (int    *)malloc(MAXXGRID*sizeof(int));
      i2    = (int    *)malloc(MAXXGRID*sizeof(int));
      j2    = (int    *)malloc(MAXXGRID*sizeof(int));
      xarea = (double *)malloc(MAXXGRID*sizeof(double));
      xclon = (double *)malloc(MAXXGRID*sizeof(double));
      xclat = (double *)malloc(MAXXGRID*sizeof(double));

      for(i=0; i<nlon1*nlat1; i++) mask1[i] = 1.0;
      
      nxgrid = create_xgrid_great_circle(&nlon1, &nlat1, &nlon2, &nlat2, lon1_in, lat1_in,
					 lon2_in, lat2_in, mask1, i1, j1, i2, j2,
					 xarea, xclon, xclat);
      printf("\n*********************************************************\n");
      printf("\n     First input grid box longitude, latitude   \n \n");
      for(i=0; i<n1_in; i++) printf(" %g,  %g \n", lon1_in[i]*R2D, lat1_in[i]*R2D);
  
      printf("\n     Second input grid box longitude, latitude \n \n");
      for(i=0; i<n2_in; i++) printf(" %g,  %g \n", lon2_in[i]*R2D, lat2_in[i]*R2D);

      printf("\n  Number of exchange grid is %d\n", nxgrid);
      for(i=0; i<nxgrid; i++) {
	printf("(i1,j1)=(%d,%d), (i2,j2)=(%d, %d), xgrid_area=%g, xgrid_clon=%g, xgrid_clat=%g\n",
	       i1[i], j1[i], i2[i], j2[i], xarea[i], xclon[i], xclat[i]);
      }

      /* comparing the area sum of exchange grid and grid1 area */
      {
	double *x1, *y1, *z1, *area1;
	double area_sum;
	int    i;
	area_sum = 0.0;
	
	for(i=0; i<nxgrid; i++) {
	  area_sum+= xarea[i];
	}

	area1 = (double *)malloc((nlon1)*(nlat1)*sizeof(double));
	get_grid_great_circle_area_(&nlon1, &nlat1, lon1_in, lat1_in, area1);      

	printf("xgrid area sum is %g, grid 1 area is %g\n", area_sum, area1[0]);
      }
	
      printf("\n");
      free(i1);
      free(i2);
      free(j1);
      free(j2);
      free(xarea);
      free(xclon);
      free(xclat);
      free(mask1);      
    }
    else {
      latlon2xyz(n1_in, lon1_in, lat1_in, x1_in, y1_in, z1_in);
      latlon2xyz(n2_in, lon2_in, lat2_in, x2_in, y2_in, z2_in);
    
      n_out = clip_2dx2d_great_circle(x1_in, y1_in, z1_in, 4, x2_in, y2_in, z2_in, n2_in,
				      x_out, y_out,  z_out);
      xyz2latlon(n_out, x_out, y_out, z_out, lon_out, lat_out);

      printf("\n*********************************************************\n");
      printf("\n     First input grid box longitude, latitude   \n \n");
      for(i=0; i<n1_in; i++) printf(" %g,  %g \n", lon1_in[i]*R2D, lat1_in[i]*R2D);
  
      printf("\n     Second input grid box longitude, latitude \n \n");
      for(i=0; i<n2_in; i++) printf(" %g,  %g \n", lon2_in[i]*R2D, lat2_in[i]*R2D);
  
      printf("\n     output clip grid box longitude, latitude for case 1 \n \n");
      for(i=0; i<n_out; i++) printf(" %g,  %g \n", lon_out[i]*R2D, lat_out[i]*R2D);
      printf("\n");
    }
  }
}


#endif
