/*
  Copyright 2011 NOAA Geophysical Fluid Dynamics Lab, Princeton, NJ
  This program is distributed under the terms of the GNU General Public
  License. See the file COPYING contained in this directory
*/
#ifndef INTERP_H_
#define INTERP_H_
/*********************************************************************
                     interp.h
   This header files contains defition of some interpolation routine  (1-D or 2-D).
   contact: Zhi.Liang@noaa.gov
*********************************************************************/
void cubic_spline_sp(int size1, int size2, const double *grid1, const double *grid2, const double *data1,
                  double *data2 );
void cubic_spline(int size1, int size2, const double *grid1, const double *grid2, const double *data1,
		  double *data2, double yp1, double ypn  );
void conserve_interp(int nx_src, int ny_src, int nx_dst, int ny_dst, const double *x_src,
		     const double *y_src, const double *x_dst, const double *y_dst,
		     const double *mask_src, const double *data_src, double *data_dst );
void conserve_interp_great_circle(int nx_src, int ny_src, int nx_dst, int ny_dst, const double *x_src,
		     const double *y_src, const double *x_dst, const double *y_dst,
		     const double *mask_src, const double *data_src, double *data_dst );
void linear_vertical_interp(int nx, int ny, int nk1, int nk2, const double *grid1, const double *grid2,
			    double *data1, double *data2);
#endif
