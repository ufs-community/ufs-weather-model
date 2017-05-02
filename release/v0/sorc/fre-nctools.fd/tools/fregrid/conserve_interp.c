#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include <math.h>
#include "constant.h"
#include "globals.h"
#include "create_xgrid.h"
#include "mosaic_util.h"
#include "conserve_interp.h"
#include "fregrid_util.h"
#include "mpp.h"
#include "mpp_io.h"
#include "read_mosaic.h"

#define  AREA_RATIO (1.e-3)
#define  MAXVAL (1.e20)
#define  TOLERANCE  (1.e-10)
/*******************************************************************************
  void setup_conserve_interp
  Setup the interpolation weight for conservative interpolation 
*******************************************************************************/
void setup_conserve_interp(int ntiles_in, const Grid_config *grid_in, int ntiles_out,
			   Grid_config *grid_out, Interp_config *interp, unsigned int opcode) 
{
  int    n, m, i, ii, jj, nx_in, ny_in, nx_out, ny_out, tile;
  size_t nxgrid, nxgrid2, nxgrid_prev;
  int    *i_in=NULL, *j_in=NULL, *i_out=NULL, *j_out=NULL;
  int   *tmp_t_in=NULL, *tmp_i_in=NULL, *tmp_j_in=NULL, *tmp_i_out=NULL, *tmp_j_out=NULL;
  double *tmp_di_in, *tmp_dj_in;
  double *xgrid_area=NULL, *tmp_area=NULL, *xgrid_clon=NULL, *xgrid_clat=NULL;
  
  double garea;
  typedef struct{
    double *area;
    double *clon;
    double *clat;
  } CellStruct;
  CellStruct *cell_in;
  
  garea = 4*M_PI*RADIUS*RADIUS;
  
  if( opcode & READ) {
    for(n=0; n<ntiles_out; n++) {
      if( interp[n].file_exist ) { /* reading from file */
	int *t_in, *ind;
	int fid, vid;
	
	nxgrid     = read_mosaic_xgrid_size(interp[n].remap_file);
	i_in       = (int    *)malloc(nxgrid   * sizeof(int   ));
	j_in       = (int    *)malloc(nxgrid   * sizeof(int   ));
	i_out      = (int    *)malloc(nxgrid   * sizeof(int   ));
	j_out      = (int    *)malloc(nxgrid   * sizeof(int   ));
	xgrid_area = (double *)malloc(nxgrid   * sizeof(double));
	if(opcode & CONSERVE_ORDER2) {
	  xgrid_clon = (double *)malloc(nxgrid   * sizeof(double));
	  xgrid_clat = (double *)malloc(nxgrid   * sizeof(double));
	}    
	t_in       = (int    *)malloc(nxgrid*sizeof(int   ));
	ind        = (int    *)malloc(nxgrid*sizeof(int   ));
	if(opcode & CONSERVE_ORDER1)
	  read_mosaic_xgrid_order1(interp[n].remap_file, i_in, j_in, i_out, j_out, xgrid_area);
	else
	  read_mosaic_xgrid_order2(interp[n].remap_file, i_in, j_in, i_out, j_out, xgrid_area, xgrid_clon, xgrid_clat);
	
	/*--- rescale the xgrid area */
	for(i=0; i<nxgrid; i++) xgrid_area[i] *= garea;
	fid = mpp_open(interp[n].remap_file, MPP_READ);
	vid = mpp_get_varid(fid, "tile1");
      	mpp_get_var_value(fid, vid, t_in);
	mpp_close(fid);
	/*distribute the exchange grid on each pe according to target grid index*/
	interp[n].nxgrid = 0;
	for(i=0; i<nxgrid; i++) {
	  if( i_out[i] <= grid_out[n].iec && i_out[i] >= grid_out[n].isc &&
	      j_out[i] <= grid_out[n].jec && j_out[i] >= grid_out[n].jsc )
	    ind[interp[n].nxgrid++] = i;
	}
	interp[n].i_in   = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
	interp[n].j_in   = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
	interp[n].i_out  = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
	interp[n].j_out  = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
	interp[n].area   = (double *)malloc(interp[n].nxgrid*sizeof(double));
	interp[n].t_in   = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));

	for(i=0; i< interp[n].nxgrid; i++) {
	  interp[n].i_in [i] = i_in [ind[i]];
	  interp[n].j_in [i] = j_in [ind[i]];
	  interp[n].t_in [i] = t_in [ind[i]] - 1;
	  interp[n].i_out[i] = i_out[ind[i]] - grid_out[n].isc;
	  interp[n].j_out[i] = j_out[ind[i]] - grid_out[n].jsc;
	  interp[n].area [i] = xgrid_area[ind[i]];
     	}
	if(opcode & CONSERVE_ORDER2) {
	  interp[n].di_in   = (double *)malloc(interp[n].nxgrid*sizeof(double));
	  interp[n].dj_in   = (double *)malloc(interp[n].nxgrid*sizeof(double));
	  for(i=0; i< interp[n].nxgrid; i++) {
	    interp[n].di_in[i] = xgrid_clon[ind[i]];
	    interp[n].dj_in[i] = xgrid_clat[ind[i]];
	  }
	}	
	free(t_in);
	free(ind);
      }
    }
    if(mpp_pe() == mpp_root_pe())printf("NOTE: Finish reading index and weight for conservative interpolation from file.\n");
  }
  else {
    i_in       = (int    *)malloc(MAXXGRID   * sizeof(int   ));
    j_in       = (int    *)malloc(MAXXGRID   * sizeof(int   ));
    i_out      = (int    *)malloc(MAXXGRID   * sizeof(int   ));
    j_out      = (int    *)malloc(MAXXGRID   * sizeof(int   ));
    xgrid_area = (double *)malloc(MAXXGRID   * sizeof(double));
    xgrid_clon = (double *)malloc(MAXXGRID   * sizeof(double));
    xgrid_clat = (double *)malloc(MAXXGRID   * sizeof(double));;
    cell_in    = (CellStruct *)malloc(ntiles_in * sizeof(CellStruct));
    for(m=0; m<ntiles_in; m++) {
      nx_in = grid_in[m].nx;
      ny_in = grid_in[m].ny;
      cell_in[m].area = (double *)malloc(nx_in*ny_in*sizeof(double));
      cell_in[m].clon = (double *)malloc(nx_in*ny_in*sizeof(double));
      cell_in[m].clat = (double *)malloc(nx_in*ny_in*sizeof(double));
      for(n=0; n<nx_in*ny_in; n++) {
	cell_in[m].area[n] = 0;
	cell_in[m].clon[n] = 0;
        cell_in[m].clat[n] = 0;
      }
    }
    for(n=0; n<ntiles_out; n++) {
      nx_out    = grid_out[n].nxc;
      ny_out    = grid_out[n].nyc;      
      interp[n].nxgrid = 0;
      for(m=0; m<ntiles_in; m++) {
	double *mask;
	double y_min, y_max, yy;
	int jstart, jend, ny_now, j;

        nx_in = grid_in[m].nx;
	ny_in = grid_in[m].ny;

	mask = (double *)malloc(nx_in*ny_in*sizeof(double));
	for(i=0; i<nx_in*ny_in; i++) mask[i] = 1.0; 

	if(opcode & GREAT_CIRCLE) {
	  nxgrid = create_xgrid_great_circle(&nx_in, &ny_in, &nx_out, &ny_out, grid_in[m].lonc,
					     grid_in[m].latc,  grid_out[n].lonc,  grid_out[n].latc,
					     mask, i_in, j_in, i_out, j_out, xgrid_area, xgrid_clon, xgrid_clat);
	  }
	else {
	  y_min = minval_double((nx_out+1)*(ny_out+1), grid_out[n].latc);
	  y_max = maxval_double((nx_out+1)*(ny_out+1), grid_out[n].latc);
	  jstart = ny_in; jend = -1;
	  for(j=0; j<=ny_in; j++) for(i=0; i<=nx_in; i++) {
	    yy = grid_in[m].latc[j*(nx_in+1)+i];
           if( yy > y_min ) {
               if(j < jstart ) jstart = j;
            }
            if( yy < y_max ) {
               if(j > jend ) jend = j;
            }

	  }
	  jstart = max(0, jstart-1);
	  jend   = min(ny_in-1, jend+1);	
	  ny_now = jend-jstart+1;
	
	  if(opcode & CONSERVE_ORDER1) {
	    nxgrid = create_xgrid_2dx2d_order1(&nx_in, &ny_now, &nx_out, &ny_out, grid_in[m].lonc+jstart*(nx_in+1),
					       grid_in[m].latc+jstart*(nx_in+1),  grid_out[n].lonc,  grid_out[n].latc,
					       mask, i_in, j_in, i_out, j_out, xgrid_area);
	    for(i=0; i<nxgrid; i++) j_in[i] += jstart;
	  }
	  else if(opcode & CONSERVE_ORDER2) {
	    int g_nxgrid;
	    int    *g_i_in, *g_j_in;
	    double *g_area, *g_clon, *g_clat;
	  
	    nxgrid = create_xgrid_2dx2d_order2(&nx_in, &ny_now, &nx_out, &ny_out, grid_in[m].lonc+jstart*(nx_in+1),
					       grid_in[m].latc+jstart*(nx_in+1),  grid_out[n].lonc,  grid_out[n].latc,
					       mask, i_in, j_in, i_out, j_out, xgrid_area, xgrid_clon, xgrid_clat);
	    for(i=0; i<nxgrid; i++) j_in[i] += jstart;
	  
	    /* For the purpose of bitiwise reproducing, the following operation is needed. */
	    g_nxgrid = nxgrid;
	    mpp_sum_int(1, &g_nxgrid);
	    if(g_nxgrid > 0) {
	      g_i_in = (int    *)malloc(g_nxgrid*sizeof(int   ));
	      g_j_in = (int    *)malloc(g_nxgrid*sizeof(int   ));			   
	      g_area = (double *)malloc(g_nxgrid*sizeof(double));
	      g_clon = (double *)malloc(g_nxgrid*sizeof(double));
	      g_clat = (double *)malloc(g_nxgrid*sizeof(double));
	      mpp_gather_field_int   (nxgrid, i_in,       g_i_in);
	      mpp_gather_field_int   (nxgrid, j_in,       g_j_in);
	      mpp_gather_field_double(nxgrid, xgrid_area, g_area);
	      mpp_gather_field_double(nxgrid, xgrid_clon, g_clon);
	      mpp_gather_field_double(nxgrid, xgrid_clat, g_clat);
	      for(i=0; i<g_nxgrid; i++) {
		ii = g_j_in[i]*nx_in+g_i_in[i];
		cell_in[m].area[ii] += g_area[i];
		cell_in[m].clon[ii] += g_clon[i];
		cell_in[m].clat[ii] += g_clat[i];
	      }
	      free(g_i_in);
	      free(g_j_in);
	      free(g_area);
	      free(g_clon);
	      free(g_clat);
	    }
	  }
	  else
	    mpp_error("conserve_interp: interp_method should be CONSERVE_ORDER1 or CONSERVE_ORDER2");
	}

      	free(mask);
	if(nxgrid > 0) {
	  nxgrid_prev = interp[n].nxgrid;
	  interp[n].nxgrid += nxgrid;
	  if(nxgrid_prev == 0 ) {
	    interp[n].i_in   = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
	    interp[n].j_in   = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
	    interp[n].i_out  = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
	    interp[n].j_out  = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
	    interp[n].area   = (double *)malloc(interp[n].nxgrid*sizeof(double));
	    interp[n].t_in   = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
	    for(i=0; i<interp[n].nxgrid; i++) {
	      interp[n].t_in [i] = m;
	      interp[n].i_in [i] = i_in [i];
	      interp[n].j_in [i] = j_in [i];
	      interp[n].i_out[i] = i_out[i];
	      interp[n].j_out[i] = j_out[i];
	      interp[n].area[i]  = xgrid_area[i];
	    }
	    if(opcode & CONSERVE_ORDER2) {
	      interp[n].di_in   = (double *)malloc(interp[n].nxgrid*sizeof(double));
	      interp[n].dj_in   = (double *)malloc(interp[n].nxgrid*sizeof(double));
	      for(i=0; i<interp[n].nxgrid; i++) {
		interp[n].di_in [i] = xgrid_clon[i]/xgrid_area[i];
		interp[n].dj_in [i] = xgrid_clat[i]/xgrid_area[i];
	      }
	    }
	  }
	  else {
	    tmp_i_in  = interp[n].i_in;
	    tmp_j_in  = interp[n].j_in;
	    tmp_i_out = interp[n].i_out;
	    tmp_j_out = interp[n].j_out;
	    tmp_area  = interp[n].area;
	    tmp_t_in  = interp[n].t_in;
	    interp[n].i_in   = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
	    interp[n].j_in   = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
	    interp[n].i_out  = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
	    interp[n].j_out  = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
	    interp[n].area   = (double *)malloc(interp[n].nxgrid*sizeof(double));
	    interp[n].t_in   = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
	    for(i=0; i<nxgrid_prev; i++) {
	      interp[n].t_in [i] = tmp_t_in [i];
	      interp[n].i_in [i] = tmp_i_in [i];
	      interp[n].j_in [i] = tmp_j_in [i];
	      interp[n].i_out[i] = tmp_i_out[i];
	      interp[n].j_out[i] = tmp_j_out[i];
	      interp[n].area [i] = tmp_area [i];
	    }
	    for(i=0; i<nxgrid; i++) {
	      ii = i + nxgrid_prev;
	      interp[n].t_in [ii] = m;
	      interp[n].i_in [ii] = i_in [i];
	      interp[n].j_in [ii] = j_in [i];
	      interp[n].i_out[ii] = i_out[i];
	      interp[n].j_out[ii] = j_out[i];
	      interp[n].area [ii] = xgrid_area[i];
	    }
	    if(opcode & CONSERVE_ORDER2) {
	      tmp_di_in  = interp[n].di_in;
	      tmp_dj_in  = interp[n].dj_in;
	      interp[n].di_in   = (double *)malloc(interp[n].nxgrid*sizeof(double));
	      interp[n].dj_in   = (double *)malloc(interp[n].nxgrid*sizeof(double));
	      for(i=0; i<nxgrid_prev; i++) { 
		interp[n].di_in [i] = tmp_di_in [i];
		interp[n].dj_in [i] = tmp_dj_in [i];
	      }
	      for(i=0; i<nxgrid; i++) {
		ii = i + nxgrid_prev;
		jj = j_in [i]*nx_in+i_in [i];
		interp[n].di_in [ii] = xgrid_clon[i]/xgrid_area[i];
		interp[n].dj_in [ii] = xgrid_clat[i]/xgrid_area[i];
	      }
	      free(tmp_di_in);
	      free(tmp_dj_in);
	    }
	    free(tmp_t_in);
	    free(tmp_i_in);
	    free(tmp_j_in);
	    free(tmp_i_out);
	    free(tmp_j_out);
	    free(tmp_area);
	  }
	}  /* if(nxgrid>0) */
      }
    }
    if(opcode & CONSERVE_ORDER2) {
      /* subtrack the grid_in clon and clat to get the distance between xgrid and grid_in */
      for(n=0; n<ntiles_in; n++) {
	double x1_in[50], y1_in[50], lon_in_avg, clon, clat;
	int    j, n0, n1, n2, n3, n1_in;
	/* calcualte cell area */
     	nx_in = grid_in[n].nx;
	ny_in = grid_in[n].ny;
	for(j=0; j<ny_in; j++) for(i=0; i<nx_in; i++) {
	  ii = j*nx_in + i;
	  if(cell_in[n].area[ii] > 0) {
	    if( fabs(cell_in[n].area[ii]-grid_in[n].cell_area[ii])/grid_in[n].cell_area[ii] < AREA_RATIO ) {
	      cell_in[n].clon[ii] /= cell_in[n].area[ii];
	      cell_in[n].clat[ii] /= cell_in[n].area[ii];
	    }
	    else {
	      n0 = j*(nx_in+1)+i;       n1 = j*(nx_in+1)+i+1;
	      n2 = (j+1)*(nx_in+1)+i+1; n3 = (j+1)*(nx_in+1)+i;
	      x1_in[0] = grid_in[n].lonc[n0]; y1_in[0] = grid_in[n].latc[n0];
	      x1_in[1] = grid_in[n].lonc[n1]; y1_in[1] = grid_in[n].latc[n1];
	      x1_in[2] = grid_in[n].lonc[n2]; y1_in[2] = grid_in[n].latc[n2];
	      x1_in[3] = grid_in[n].lonc[n3]; y1_in[3] = grid_in[n].latc[n3];
	      n1_in = fix_lon(x1_in, y1_in, 4, M_PI);
	      lon_in_avg = avgval_double(n1_in, x1_in);
              clon = poly_ctrlon(x1_in, y1_in, n1_in, lon_in_avg);
	      clat = poly_ctrlat (x1_in, y1_in, n1_in );
	      cell_in[n].clon[ii] = clon/grid_in[n].cell_area[ii];
	      cell_in[n].clat[ii] = clat/grid_in[n].cell_area[ii];
	    }
	  } 
	}
      }
      for(n=0; n<ntiles_out; n++) {
	for(i=0; i<interp[n].nxgrid; i++) {
	  tile = interp[n].t_in[i];
	  ii   = interp[n].j_in[i] * grid_in[tile].nx + interp[n].i_in[i];
          interp[n].di_in[i] -= cell_in[tile].clon[ii];
	  interp[n].dj_in[i] -= cell_in[tile].clat[ii];
	}
      }

      /* free the memory */
      for(n=0; n<ntiles_in; n++) {
	free(cell_in[n].area);
	free(cell_in[n].clon);
	free(cell_in[n].clat);
      }
      free(cell_in);
    }
    if( opcode & WRITE) { /* write out remapping information */
      for(n=0; n<ntiles_out; n++) {
	int nxgrid;
	
	nxgrid = interp[n].nxgrid;
	mpp_sum_int(1, &nxgrid);
	if(nxgrid > 0) {
	  size_t start[4], nwrite[4];
	  int    fid, dim_string, dim_ncells, dim_two, dims[4];
	  int    id_xgrid_area, id_tile1_dist;
	  int    id_tile1_cell, id_tile2_cell, id_tile1;
	  int    *gdata_int, *ldata_int;	  
	  double *gdata_dbl;
	  
	  fid = mpp_open( interp[n].remap_file, MPP_WRITE);
	  dim_string = mpp_def_dim(fid, "string", STRING);
	  dim_ncells = mpp_def_dim(fid, "ncells", nxgrid);
	  dim_two    = mpp_def_dim(fid, "two", 2);
	  dims[0] = dim_ncells; dims[1] = dim_two;
	  id_tile1      = mpp_def_var(fid, "tile1",      NC_INT, 1, &dim_ncells, 1,
				      "standard_name", "tile_number_in_mosaic1");
	  id_tile1_cell = mpp_def_var(fid, "tile1_cell", NC_INT, 2, dims, 1,
				      "standard_name", "parent_cell_indices_in_mosaic1");
	  id_tile2_cell = mpp_def_var(fid, "tile2_cell", NC_INT, 2, dims, 1,
				      "standard_name", "parent_cell_indices_in_mosaic2");
	  id_xgrid_area = mpp_def_var(fid, "xgrid_area", NC_DOUBLE, 1, &dim_ncells, 2,
				      "standard_name", "exchange_grid_area", "units", "m2");
	  if(opcode & CONSERVE_ORDER2) id_tile1_dist = mpp_def_var(fid, "tile1_distance", NC_DOUBLE, 2, dims, 1,
								   "standard_name", "distance_from_parent1_cell_centroid");
	  mpp_end_def(fid);
	  for(i=0; i<4; i++) {
	    start[i] = 0; nwrite[i] = 1;
	  }
	  nwrite[0] = nxgrid;
	  gdata_int = (int *)malloc(nxgrid*sizeof(int));
	  if(interp[n].nxgrid>0) ldata_int = (int *)malloc(interp[n].nxgrid*sizeof(int));
          mpp_gather_field_int(interp[n].nxgrid, interp[n].t_in, gdata_int);
	  for(i=0; i<nxgrid; i++) gdata_int[i]++;
	  mpp_put_var_value(fid, id_tile1, gdata_int);

	  mpp_gather_field_int(interp[n].nxgrid, interp[n].i_in, gdata_int);
	  for(i=0; i<nxgrid; i++) gdata_int[i]++;
	  mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, gdata_int);

          for(i=0; i<interp[n].nxgrid; i++) ldata_int[i] = interp[n].i_out[i] + grid_out[n].isc + 1; 
	  mpp_gather_field_int(interp[n].nxgrid, ldata_int, gdata_int);
	  mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, gdata_int);

	  mpp_gather_field_int(interp[n].nxgrid, interp[n].j_in, gdata_int);
	  for(i=0; i<nxgrid; i++) gdata_int[i]++;
	  start[1] = 1;
	  mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, gdata_int);

          for(i=0; i<interp[n].nxgrid; i++) ldata_int[i] = interp[n].j_out[i] + grid_out[n].jsc + 1; 	  
	  mpp_gather_field_int(interp[n].nxgrid, ldata_int, gdata_int);
	  mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, gdata_int);

	  free(gdata_int);
	  if(interp[n].nxgrid>0)free(ldata_int);
	  
	  gdata_dbl = (double *)malloc(nxgrid*sizeof(double));
	  mpp_gather_field_double(interp[n].nxgrid, interp[n].area, gdata_dbl);
	  mpp_put_var_value(fid, id_xgrid_area, gdata_dbl);
	  
	  if(opcode & CONSERVE_ORDER2) {
	    start[1] = 0;
	    mpp_gather_field_double(interp[n].nxgrid, interp[n].di_in, gdata_dbl);
	    mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, gdata_dbl);
	    start[1] = 1;
	    mpp_gather_field_double(interp[n].nxgrid, interp[n].dj_in, gdata_dbl);
	    mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, gdata_dbl);
	  }
	  
	  free(gdata_dbl);
	  mpp_close(fid);
	}
      }
    }
    if(mpp_pe() == mpp_root_pe())printf("NOTE: done calculating index and weight for conservative interpolation\n");
  }

  /* check the input area match exchange grid area */
  if(opcode & CHECK_CONSERVE) {
    int nx1, ny1, max_i, max_j, i, j;
    double max_ratio, ratio_change;
    double *area2;
    
    /* sum over exchange grid to get the area of grid_in */
    nx1  = grid_out[0].nxc;
    ny1  = grid_out[0].nyc;

    area2 = (double *)malloc(nx1*ny1*sizeof(double));
    
    for(n=0; n<ntiles_out; n++) {
      for(i=0; i<nx1*ny1; i++) area2[i] = 0;
      for(i=0; i<interp[n].nxgrid; i++) {
	ii = interp[n].j_out[i]*nx1 + interp[n].i_out[i];
	area2[ii] +=  interp[n].area[i];
      }
      max_ratio = 0;
      max_i = 0;
      max_j = 0;
      /* comparing area1 and area2 */
      for(j=0; j<ny1; j++) for(i=0; i<nx1; i++) {
	ii = j*nx1+i;
	ratio_change = fabs(grid_out[n].cell_area[ii]-area2[ii])/grid_out[n].cell_area[ii];
	if(ratio_change > max_ratio) {
	  max_ratio = ratio_change;
	  max_i = i;
	  max_j = j;
	}
	if( ratio_change > 1.e-4 ) {
	  printf("(i,j)=(%d,%d), change = %g, area1=%g, area2=%g\n", i, j, ratio_change, grid_out[n].cell_area[ii],area2[ii]);
	}
      }
      ii = max_j*nx1+max_i;
      printf("The maximum ratio change at (%d,%d) = %g, area1=%g, area2=%g\n", max_i, max_j, max_ratio, grid_out[n].cell_area[ii],area2[ii]);
      
    }
    
    free(area2);
    
  }
      
  free(i_in);
  free(j_in);
  free(i_out);
  free(j_out);
  free(xgrid_area);
  if(xgrid_clon) free(xgrid_clon);
  if(xgrid_clat) free(xgrid_clat);
  
}; /* setup_conserve_interp */


/*******************************************************************************
 void do_scalar_conserve_interp( )
 doing conservative interpolation
*******************************************************************************/
void do_scalar_conserve_interp(Interp_config *interp, int varid, int ntiles_in, const Grid_config *grid_in,
			       int ntiles_out, const Grid_config *grid_out, const Field_config *field_in,
			       Field_config *field_out, unsigned int opcode, int nz)
{
  int nx1, ny1, nx2, ny2, i1, j1, i2, j2, tile, n, m, i, j, n1, n2;
  int k, n0;
  int has_missing, halo, interp_method;
  int weight_exist;
  int cell_measures, cell_methods;
  double area, missing, di, dj, area_missing;
  double *out_area;
  double gsum_out;
  int monotonic;
  Monotone_config *monotone_data;
  
  gsum_out = 0;
  interp_method = field_in->var[varid].interp_method;
  halo = 0;
  monotonic = 0;
  if(interp_method == CONSERVE_ORDER2) {
    halo = 1;
    monotonic = opcode & MONOTONIC;
  }
    
  area_missing = field_in->var[varid].area_missing;
  has_missing = field_in->var[varid].has_missing;
  weight_exist = grid_in[0].weight_exist;
  cell_measures = field_in->var[varid].cell_measures;
  cell_methods = field_in->var[varid].cell_methods;

  missing = -MAXVAL;
  if(has_missing) missing = field_in->var[varid].missing;
  
  if( nz>1 && has_missing ) mpp_error("conserve_interp: has_missing should be false when nz > 1");
  if( nz>1 && cell_measures ) mpp_error("conserve_interp: cell_measures should be false when nz > 1");
  if( nz>1 && cell_methods == CELL_METHODS_SUM ) mpp_error("conserve_interp: cell_methods should not be sum when nz > 1");  
  /*  if( nz>1 && monotonic ) mpp_error("conserve_interp: monotonic should be false when nz > 1"); */
  
  if(monotonic) monotone_data = (Monotone_config *)malloc(ntiles_in*sizeof(Monotone_config));
  
  for(m=0; m<ntiles_out; m++) {
    nx2 = grid_out[m].nxc;
    ny2 = grid_out[m].nyc;
    out_area = (double *)malloc(nx2*ny2*nz*sizeof(double));
    
    for(i=0; i<nx2*ny2*nz; i++) field_out[m].data[i] = 0.0;

    for(i=0; i<nx2*ny2; i++) out_area[i] = 0.0;
    if(interp_method == CONSERVE_ORDER1) {
      if(has_missing) {
	for(n=0; n<interp[m].nxgrid; n++) {
	  i2   = interp[m].i_out[n];
	  j2   = interp[m].j_out[n];
	  i1   = interp[m].i_in [n];
	  j1   = interp[m].j_in [n];    
	  tile = interp[m].t_in [n];
	  area = interp[m].area [n];	  
	  nx1  = grid_in[tile].nx;
	  ny1  = grid_in[tile].ny;
          if(weight_exist) area *= grid_in[tile].weight[j1*nx1+i1];
	  n1 = j1*nx1+i1;
	  n0 = j2*nx2+i2;
	  
	  if( field_in[tile].data[n1] != missing ) {
            if( cell_measures )
	      area *= (field_in[tile].area[n1]/grid_in[tile].cell_area[n1]);
	    else if( cell_methods == CELL_METHODS_SUM )
	      area /= grid_in[tile].cell_area[n1];
  	    field_out[m].data[n0] += (field_in[tile].data[n1]*area);
            out_area[n0] += area;
          }
          else if(cell_measures) {
            if( field_in[tile].area[n1] != area_missing ){
              area *= (field_in[tile].area[n1]/grid_in[tile].cell_area[n1]);
              out_area[n0] += area;
            }
          }
        }
      }
      else {
	for(n=0; n<interp[m].nxgrid; n++) {
	  i2   = interp[m].i_out[n];
	  j2   = interp[m].j_out[n];
	  i1   = interp[m].i_in [n];
	  j1   = interp[m].j_in [n];    
	  tile = interp[m].t_in [n];
	  area = interp[m].area [n];
	  nx1  = grid_in[tile].nx;
	  ny1  = grid_in[tile].ny;
	  if(weight_exist) area *= grid_in[tile].weight[j1*nx1+i1];
	  for(k=0; k<nz; k++) {
	    n1 = k*nx1*ny1 + j1*nx1+i1;
	    n0 = k*nx2*ny2 + j2*nx2+i2;
            if( cell_measures )
	      area *= (field_in[tile].area[n1]/grid_in[tile].cell_area[n1]);
	    else if(  cell_methods == CELL_METHODS_SUM )
	      area /= grid_in[tile].cell_area[n1];
  	    field_out[m].data[n0] += (field_in[tile].data[n1]*area);
	    out_area[n0] += area;
	  }
	}
      } 
    }
    else if(monotonic) {
      int ii, jj;
      double f_bar;
      double *xdata;
      for(n=0; n<ntiles_in; n++) {
	nx1 =  grid_in[n].nx;
	ny1 =  grid_in[n].ny;
	monotone_data[n].f_bar_max = (double *)malloc(nx1*ny1*sizeof(double));
	monotone_data[n].f_bar_min = (double *)malloc(nx1*ny1*sizeof(double));
	monotone_data[n].f_max     = (double *)malloc(nx1*ny1*sizeof(double));
	monotone_data[n].f_min     = (double *)malloc(nx1*ny1*sizeof(double));
	for(j=0; j<ny1; j++) for(i=0; i<nx1; i++) {
	  n1 = j*nx1+i;
	    
	  monotone_data[n].f_bar_max[n1] = -MAXVAL;
	  monotone_data[n].f_bar_min[n1] = MAXVAL;
	  monotone_data[n].f_max[n1]     = -MAXVAL;
	  monotone_data[n].f_min[n1]     = MAXVAL;
	  n1 = j*nx1+i;
	  for(jj=j-1; jj<=j+1; jj++) for(ii=i-1; ii<=i+1; ii++) {
	    n2 = (jj+1)*(nx1+2)+ii+1;
	    if( field_in[n].data[n2] != missing ) {
	      if( field_in[n].data[n2] > monotone_data[n].f_bar_max[n1] ) monotone_data[n].f_bar_max[n1] = field_in[n].data[n2];
	      if( field_in[n].data[n2] < monotone_data[n].f_bar_min[n1] ) monotone_data[n].f_bar_min[n1] = field_in[n].data[n2];
	    }
	  }
	}
      }      
      
      xdata = (double *)malloc(interp[m].nxgrid*sizeof(double));
      for(n=0; n<interp[m].nxgrid; n++) {
	i1   = interp[m].i_in [n];
	j1   = interp[m].j_in [n];
	di   = interp[m].di_in[n];
	dj   = interp[m].dj_in[n];
	tile = interp[m].t_in [n];
	n1 = j1*nx1+i1;
        n2 = (j1+1)*(nx1+2)+i1+1;
	if( field_in[tile].data[n2] != missing ) {
	  if( field_in[tile].grad_mask[n1] ) { /* use zero gradient */
	    xdata[n] = field_in[tile].data[n2];
	  }
	  else {
	    xdata[n] = field_in[tile].data[n2]+field_in[tile].grad_x[n1]*di+field_in[tile].grad_y[n1]*dj;
	  }
	  if(monotonic) {
	    if( xdata[n] > monotone_data[tile].f_max[n1]) monotone_data[tile].f_max[n1] = xdata[n];
	    if( xdata[n] < monotone_data[tile].f_min[n1]) monotone_data[tile].f_min[n1] = xdata[n];
	  }
	}
	else
	  xdata[n] = missing;
      }

      /* get the global f_max and f_min */
      if(mpp_npes() >1) {
	for(n=0; n<ntiles_in; n++) {
	  mpp_min_double(nx1*ny1, monotone_data[n].f_min);
	  mpp_max_double(nx1*ny1, monotone_data[n].f_max);
	}
      }

      /* adjust the exchange grid cell data to make it monotonic */
      for(n=0; n<interp[m].nxgrid; n++) {
	i1   = interp[m].i_in [n];
	j1   = interp[m].j_in [n];
	tile = interp[m].t_in [n];
	n1 = j1*nx1+i1;
	n2 = (j1+1)*(nx1+2)+i1+1;
	f_bar = field_in[tile].data[n2];
	if(xdata[n] == missing) continue;
	  
	if( monotone_data[tile].f_max[n1] > monotone_data[tile].f_bar_max[n1] ) {
	  /* z1l: Due to truncation error, we might get xdata[n] > f_bar_max[n1]. So
	     we allow some tolerance. What is the suitable tolerance? */
	  xdata[n] = f_bar + ((xdata[n]-f_bar)/(monotone_data[tile].f_max[n1]-f_bar))
	    * (monotone_data[tile].f_bar_max[n1]-f_bar);
	  if( xdata[n] > monotone_data[tile].f_bar_max[n1]) {
	    if(xdata[n] - monotone_data[tile].f_bar_max[n1] < TOLERANCE ) xdata[n] = monotone_data[tile].f_bar_max[n1];
	    if( xdata[n] > monotone_data[tile].f_bar_max[n1]) {
	      printf(" n = %d, n1 = %d, xdata = %f, f_bar_max=%f\n", n, n1, xdata[n], monotone_data[tile].f_bar_max[n1]);
	      mpp_error(" xdata is greater than f_bar_max ");
	    }
	  }
	}
	else if( monotone_data[tile].f_min[n1] < monotone_data[tile].f_bar_min[n1] ) {
	  /* z1l: Due to truncation error, we might get xdata[n] < f_bar_min[n1]. So
	     we allow some tolerance. What is the suitable tolerance? */
	  xdata[n] = f_bar + ((xdata[n]-f_bar)/(monotone_data[tile].f_min[n1]-f_bar)) * (monotone_data[tile].f_bar_min[n1]-f_bar);
	  if( xdata[n] < monotone_data[tile].f_bar_min[n1]) {
	    if(monotone_data[tile].f_bar_min[n1] - xdata[n]< TOLERANCE ) xdata[n] = monotone_data[tile].f_bar_min[n1];
	    if( xdata[n] < monotone_data[tile].f_bar_min[n1]) {
	      printf(" n = %d, n1 = %d, xdata = %f, f_bar_min=%f\n", n, n1, xdata[n], monotone_data[tile].f_bar_min[n1]);
	      mpp_error(" xdata is less than f_bar_min ");
	    }
	  }	     
	}
      }
      for(n=0; n<ntiles_in; n++) {	
	free(monotone_data[n].f_bar_max);
	free(monotone_data[n].f_bar_min);
	free(monotone_data[n].f_max);
	free(monotone_data[n].f_min);
      }

      /* remap onto destination grid */
      for(n=0; n<interp[m].nxgrid; n++) {
	i2   = interp[m].i_out[n];
	j2   = interp[m].j_out[n];
	i1   = interp[m].i_in [n];
	j1   = interp[m].j_in [n];
	tile = interp[m].t_in [n];
	area = interp[m].area [n];
	if(xdata[n] == missing) continue;
	if(weight_exist) area *= grid_in[tile].weight[j1*nx1+i1];
	n1 = j1*nx1+i1;
	n0 = j2*nx2+i2;
	if( cell_measures )
	  area *= (field_in[tile].area[n1]/grid_in[tile].cell_area[n1]);
	else if( cell_methods == CELL_METHODS_SUM )
	  area /= grid_in[tile].cell_area[n1];	    
	field_out[m].data[n0] += xdata[n]*area;
	out_area[n0] += area;
      }
      free(xdata);
    }
    else {
      if(has_missing) {
	for(n=0; n<interp[m].nxgrid; n++) {
	  i2   = interp[m].i_out[n];
	  j2   = interp[m].j_out[n];
	  i1   = interp[m].i_in [n];
	  j1   = interp[m].j_in [n];
	  di   = interp[m].di_in[n];
	  dj   = interp[m].dj_in[n];
	  tile = interp[m].t_in [n];
	  area = interp[m].area [n];
	  nx1  = grid_in[tile].nx;
	  ny1  = grid_in[tile].ny;

	  if(weight_exist) area *= grid_in[tile].weight[j1*nx1+i1];
	  n2 = (j1+1)*(nx1+2)+i1+1;
	  n0 = j2*nx2+i2;
	  
	  if( field_in[tile].data[n2] != missing ) {
	    n1 = j1*nx1+i1;
	    n0 = j2*nx2+i2;
            if( cell_measures )
	      area *= (field_in[tile].area[n1]/grid_in[tile].cell_area[n1]);
	    else if( cell_methods == CELL_METHODS_SUM )
	      area /= grid_in[tile].cell_area[n1];	    
	    if(field_in[tile].grad_mask[n1]) { /* use zero gradient */
	      field_out[m].data[n0] += field_in[tile].data[n2]*area;
	    }
	    else {
	      field_out[m].data[n0] += (field_in[tile].data[n2]+field_in[tile].grad_x[n1]*di
					+field_in[tile].grad_y[n1]*dj)*area;
	    }
	    out_area[n0] += area;
	  }
	}
      }
      else {
	for(n=0; n<interp[m].nxgrid; n++) {
	  i2   = interp[m].i_out[n];
	  j2   = interp[m].j_out[n];
	  i1   = interp[m].i_in [n];
	  j1   = interp[m].j_in [n];
	  di   = interp[m].di_in[n];
	  dj   = interp[m].dj_in[n];
	  tile = interp[m].t_in [n];
	  area = interp[m].area [n];

	  nx1  = grid_in[tile].nx;
	  ny1  = grid_in[tile].ny;
	  if(weight_exist) area *= grid_in[tile].weight[j1*nx1+i1];
	  for(k=0; k<nz; k++) {
	    n0 = k*nx2*ny2 + j2*nx2+i2;
	    n1 = k*nx1*ny1+j1*nx1+i1;
	    n2 = k*(nx1+2)*(ny1+2)+(j1+1)*(nx1+2)+i1+1;
	    if( cell_measures )
	      area *= (field_in[tile].area[n1]/grid_in[tile].cell_area[n1]);
	    else if( cell_methods == CELL_METHODS_SUM )
	      area /= grid_in[tile].cell_area[n1];	  
	    field_out[m].data[n0] += (field_in[tile].data[n2]+field_in[tile].grad_x[n1]*di
				      +field_in[tile].grad_y[n1]*dj)*area;
	    out_area[n0] += area;	    
	  }
	}
      }
    }

    if(opcode & CHECK_CONSERVE) {
      for(i=0; i<nx2*ny2*nz; i++) {
	if(out_area[i] > 0) gsum_out += field_out[m].data[i];
      }
    }

    if( cell_measures || ( !(opcode & TARGET) && !(cell_methods == CELL_METHODS_SUM)) ) {
      for(i=0; i<nx2*ny2*nz; i++) {
	if(out_area[i] > 0)
	  field_out[m].data[i] /= out_area[i];
	else
	  field_out[m].data[i] = missing;      
      }
    }
    else if( opcode & TARGET ) {
      for(i=0; i<nx2*ny2; i++) {
	if(out_area[i] > 0)
	  for(k=0; k<nz; k++) field_out[m].data[k*nx2*ny2+i] /= grid_out[m].cell_area[i];
	else
	  for(k=0; k<nz; k++) field_out[m].data[k*nx2*ny2+i] = missing;
      }
    }
    else if ( cell_methods == CELL_METHODS_SUM ) {
      for(i=0; i<nx2*ny2*nz; i++) {
        if(out_area[i] == 0)
          for(k=0; k<nz; k++) field_out[m].data[k*nx2*ny2+i] = missing;
      }
    }

    free(out_area);
  }
    

  /* conservation check if needed */
  if(opcode & CHECK_CONSERVE) {
    double gsum_in, dd;
    gsum_in = 0;
    for(n=0; n<ntiles_in; n++) {

      nx1  = grid_in[n].nx;
      ny1  = grid_in[n].ny;


      if( cell_measures ) {
        for(j=0; j<ny1; j++) for(i=0; i<nx1; i++) {
  	  dd = field_in[n].data[(j+halo)*(nx1+2*halo)+i+halo];
	  if(dd != missing) gsum_in += dd*field_in[n].area[j*nx1+i];
        }      
      }
      else if ( cell_methods == CELL_METHODS_SUM ) {
        for(j=0; j<ny1; j++) for(i=0; i<nx1; i++) {
	  dd = field_in[n].data[(j+halo)*(nx1+2*halo)+i+halo];
	  if(dd != missing) gsum_in += dd;
        }
      }     
      else {
        for(k=0; k<nz; k++) for(j=0; j<ny1; j++) for(i=0; i<nx1; i++) {
	  dd = field_in[n].data[k*(nx1+2*halo)*(ny1+2*halo)+(j+halo)*(nx1+2*halo)+i+halo];
	  if(dd != missing) gsum_in += dd*grid_in[n].cell_area[j*nx1+i];
        }
      }
    }
    mpp_sum_double(1, &gsum_out);
    
    if(mpp_pe() == mpp_root_pe()) printf("the flux(data*area) sum of %s: input = %g, output = %g, diff = %g. \n",
					 field_in->var[varid].name, gsum_in, gsum_out, gsum_out-gsum_in);
    
  }
  
  
}; /* do_scalar_conserve_interp */


/*******************************************************************************
 void do_vector_conserve_interp( )
 doing conservative interpolation
*******************************************************************************/
void do_vector_conserve_interp(Interp_config *interp, int varid, int ntiles_in, const Grid_config *grid_in, int ntiles_out, 
                               const Grid_config *grid_out, const Field_config *u_in,  const Field_config *v_in,
                               Field_config *u_out, Field_config *v_out, unsigned int opcode)  
{
  int          nx1, ny1, nx2, ny2, i1, j1, i2, j2, tile, n, m, i;
  double       area, missing, tmp_x, tmp_y;
  double       *out_area;

  missing = u_in->var[varid].missing;  
  /* first rotate input data */
  for(n = 0; n < ntiles_in; n++) {
    if(grid_in[n].rotate) {
      nx1 = grid_in[n].nx;
      ny1 = grid_in[n].ny;
      for(i=0; i<nx1*ny1; i++) {
	tmp_x = u_in[n].data[i];
	tmp_y = v_in[n].data[i];
	if( tmp_x != missing && tmp_y != missing) {
	  u_in[n].data[i] = tmp_x * grid_in[n].cosrot[i] - tmp_y * grid_in[n].sinrot[i];
	  v_in[n].data[i] = tmp_x * grid_in[n].sinrot[i] + tmp_y * grid_in[n].cosrot[i];
	}
      }
    }
  }
  
  for(m=0; m<ntiles_out; m++) {
    nx2 = grid_out[m].nxc;
    ny2 = grid_out[m].nyc;
    out_area = (double *)malloc(nx2*ny2*sizeof(double));
    
    for(i=0; i<nx2*ny2; i++) {
      u_out[m].data[i] = 0.0;
      v_out[m].data[i] = 0.0;
    }
    for(i=0; i<nx2*ny2; i++) out_area[i] = 0.0;
    
    for(n=0; n<interp[m].nxgrid; n++) {
      i2   = interp[m].i_out[n];
      j2   = interp[m].j_out[n];
      i1   = interp[m].i_in [n];
      j1   = interp[m].j_in [n];    
      tile = interp[m].t_in [n];
      area = interp[m].area [n];
      nx1  = grid_in[tile].nx;
      ny1  = grid_in[tile].ny;
      tmp_x = u_in[tile].data[j1*nx1+i1];
      tmp_y = v_in[tile].data[j1*nx1+i1];
      if( tmp_x != missing && tmp_y != missing ) {
	u_out[m].data[j2*nx2+i2] += u_in[tile].data[j1*nx1+i1]*area;
	v_out[m].data[j2*nx2+i2] += v_in[tile].data[j1*nx1+i1]*area;
	out_area[j2*nx2+i2] += area;
      }
    }
    if(opcode & TARGET) {
      for(i=0; i<nx2*ny2; i++) {
	if(out_area[i] > 0) {
	  u_out[m].data[i] /= grid_out[m].area[i];
	  v_out[m].data[i] /= grid_out[m].area[i];
	}
	else {
	  u_out[m].data[i] = missing;
	  v_out[m].data[i] = missing;
	}
      }
    }
    else {
      for(i=0; i<nx2*ny2; i++) {
	if(out_area[i] > 0) {
	  u_out[m].data[i] /= out_area[i];
	  v_out[m].data[i] /= out_area[i];
	}
	else {
	  u_out[m].data[i] = missing;
	  v_out[m].data[i] = missing;
	}
      }      
    }
    /* rotate the data if needed */
    if(grid_out[m].rotate) {
      for(i=0; i<nx2*ny2; i++) {
	tmp_x = u_out[m].data[i];
	tmp_y = v_out[m].data[i];
	if( tmp_x != missing && tmp_y != missing) {
	  u_out[m].data[i] =  tmp_x * grid_out[m].cosrot[i] + tmp_y * grid_out[m].sinrot[i];
	  v_out[m].data[i] = -tmp_x * grid_out[m].sinrot[i] + tmp_y * grid_out[m].cosrot[i];
	}
      }
    }
    free(out_area);
  }
  
}; /* do_vector_conserve_interp */
