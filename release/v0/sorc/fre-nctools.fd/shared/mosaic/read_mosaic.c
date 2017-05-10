#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "read_mosaic.h"
#include "constant.h"
#include "mosaic_util.h"
#ifdef use_netCDF
#include <netcdf.h>
#endif
/*********************************************************************
    void netcdf_error( int status )
    status is the returning value of netcdf call. this routine will
    handle the error when status is not NC_NOERR.
********************************************************************/
void handle_netcdf_error(const char *msg, int status )
{
  char errmsg[512];

  sprintf( errmsg, "%s: %s", msg, (char *)nc_strerror(status) );
  error_handler(errmsg);

}; /* handle_netcdf_error */

/***************************************************************************
  void get_file_dir(const char *file, char *dir)
  get the directory where file is located. The dir will be the complate path
  before the last "/". If no "/" exist in file, the path will be current ".".
***************************************************************************/
void get_file_dir(const char *file, char *dir)
{
  int len;
  char *strptr = NULL;

  /* get the diretory */
 
  strptr = strrchr(file, '/');
  if(strptr) {
    len = strptr - file;
    strncpy(dir, file, len);
  }
  else {
    len = 1;
    strcpy(dir, ".");
  }
  dir[len] = 0;

}; /* get_file_dir */


int field_exist(const char* file, const char *name)
{
  int ncid, varid, status, existed;
  char msg[512];  
#ifdef use_netCDF
  status = nc_open(file, NC_NOWRITE, &ncid);
  if(status != NC_NOERR) {
    sprintf(msg, "field_exist: in opening file %s", file);
    handle_netcdf_error(msg, status);
  }
  status = nc_inq_varid(ncid, name, &varid);  
  if(status == NC_NOERR)
    existed = 1;
  else
    existed = 0;
    
  status = nc_close(ncid);
  if(status != NC_NOERR) {
    sprintf(msg, "field_exist: in closing file %s.", file);
    handle_netcdf_error(msg, status);
  }

  return existed;
#else
  error_handler("read_mosaic: Add flag -Duse_netCDF when compiling");

#endif
  return 0; 
}; /* field_exist */

int get_dimlen(const char* file, const char *name)
{
  int ncid, dimid, status, len;
  size_t size;
  char msg[512];

  len = 0;
#ifdef use_netCDF  
  status = nc_open(file, NC_NOWRITE, &ncid);
  if(status != NC_NOERR) {
    sprintf(msg, "in opening file %s", file);
    handle_netcdf_error(msg, status);
  }
  
  status = nc_inq_dimid(ncid, name, &dimid);
  if(status != NC_NOERR) {
    sprintf(msg, "in getting dimid of %s from file %s.", name, file);
    handle_netcdf_error(msg, status);
  }
  
  status = nc_inq_dimlen(ncid, dimid, &size);
  if(status != NC_NOERR) {
    sprintf(msg, "in getting dimension size of %s from file %s.", name, file);
    handle_netcdf_error(msg, status);
  }
  status = nc_close(ncid);
  if(status != NC_NOERR) {
    sprintf(msg, "in closing file %s.", file);
    handle_netcdf_error(msg, status);
  }
  
  len = size;
  if(status != NC_NOERR) {
    sprintf(msg, "in closing file %s", file);
    handle_netcdf_error(msg, status);
  }
#else
  error_handler("read_mosaic: Add flag -Duse_netCDF when compiling");
#endif
  
  return len;
  
}; /* get_dimlen */

/*******************************************************************************
   void get_string_data(const char *file, const char *name, char *data)
   get string data of field with "name" from "file".
******************************************************************************/
void get_string_data(const char *file, const char *name, char *data)
{
  int ncid, varid, status;
  char msg[512];

#ifdef use_netCDF    
  status = nc_open(file, NC_NOWRITE, &ncid);
  if(status != NC_NOERR) {
    sprintf(msg, "in opening file %s", file);
    handle_netcdf_error(msg, status);
  }
  status = nc_inq_varid(ncid, name, &varid);
  if(status != NC_NOERR) {
    sprintf(msg, "in getting varid of %s from file %s.", name, file);
    handle_netcdf_error(msg, status);
  }     
  status = nc_get_var_text(ncid, varid, data);
  if(status != NC_NOERR) {
    sprintf(msg, "in getting data of %s from file %s.", name, file);
    handle_netcdf_error(msg, status);
  }
  status = nc_close(ncid);
  if(status != NC_NOERR) {
    sprintf(msg, "in closing file %s.", file);
    handle_netcdf_error(msg, status);
  }  
#else
  error_handler("read_mosaic: Add flag -Duse_netCDF when compiling");
#endif
  
}; /* get_string_data */

/*******************************************************************************
   void get_string_data_level(const char *file, const char *name, const size_t *start, const size_t *nread, char *data)
   get string data of field with "name" from "file".
******************************************************************************/
void get_string_data_level(const char *file, const char *name, char *data, const int *level)
{
  int ncid, varid, status, i;
  size_t start[4], nread[4];
  char msg[512];

#ifdef use_netCDF  
  status = nc_open(file, NC_NOWRITE, &ncid);
  if(status != NC_NOERR) {
    sprintf(msg, "in opening file %s", file);
    handle_netcdf_error(msg, status);
  }
  status = nc_inq_varid(ncid, name, &varid);
  if(status != NC_NOERR) {
    sprintf(msg, "in getting varid of %s from file %s.", name, file);
    handle_netcdf_error(msg, status);
  }
  for(i=0; i<4; i++) {
    start[i] = 0; nread[i] = 1;
  }
  start[0] = *level; nread[1] = STRING;
  status = nc_get_vara_text(ncid, varid, start, nread, data);
  if(status != NC_NOERR) {
    sprintf(msg, "in getting data of %s from file %s.", name, file);
    handle_netcdf_error(msg, status);
  }
  status = nc_close(ncid);
  if(status != NC_NOERR) {
    sprintf(msg, "in closing file %s.", file);
    handle_netcdf_error(msg, status);
  }  
#else
  error_handler("read_mosaic: Add flag -Duse_netCDF when compiling");
#endif
  
}; /* get_string_data_level */


/*******************************************************************************
   void get_var_data(const char *file, const char *name, double *data)
   get var data of field with "name" from "file".
******************************************************************************/
void get_var_data(const char *file, const char *name, void *data)
{

  int ncid, varid, status;  
  nc_type vartype;
  char msg[512];

#ifdef use_netCDF    
  status = nc_open(file, NC_NOWRITE, &ncid);
  if(status != NC_NOERR) {
    sprintf(msg, "in opening file %s", file);
    handle_netcdf_error(msg, status);
  }
  status = nc_inq_varid(ncid, name, &varid);
  if(status != NC_NOERR) {
    sprintf(msg, "in getting varid of %s from file %s.", name, file);
    handle_netcdf_error(msg, status);
  }

  status = nc_inq_vartype(ncid, varid, &vartype);
  if(status != NC_NOERR) {
    sprintf(msg, "get_var_data: in getting vartype of of %s in file %s ", name, file);
    handle_netcdf_error(msg, status);
  }

  switch (vartype) {
  case NC_DOUBLE:case NC_FLOAT:
#ifdef OVERLOAD_R4
  status = nc_get_var_float(ncid, varid, data);
#else
  status = nc_get_var_double(ncid, varid, data);
#endif
  break;
  case NC_INT:
    status = nc_get_var_int(ncid, varid, data);
    break;
  default:
    sprintf(msg, "get_var_data: field %s in file %s has an invalid type, "
            "the type should be NC_DOUBLE, NC_FLOAT or NC_INT", name, file);
    error_handler(msg);
  }
  if(status != NC_NOERR) {
    sprintf(msg, "in getting data of %s from file %s.", name, file);
    handle_netcdf_error(msg, status);
  }
  status = nc_close(ncid);
  if(status != NC_NOERR) {
    sprintf(msg, "in closing file %s.", file);
    handle_netcdf_error(msg, status);
  }  
#else
  error_handler("read_mosaic: Add flag -Duse_netCDF when compiling");
#endif
  
}; /* get_var_data */

/*******************************************************************************
   void get_var_data(const char *file, const char *name, double *data)
   get var data of field with "name" from "file".
******************************************************************************/
void get_var_data_region(const char *file, const char *name, const size_t *start, const size_t *nread, void *data)
{

  int ncid, varid, status;  
  nc_type vartype;
  char msg[512];

#ifdef use_netCDF    
  status = nc_open(file, NC_NOWRITE, &ncid);
  if(status != NC_NOERR) {
    sprintf(msg, "get_var_data_region: in opening file %s", file);
    handle_netcdf_error(msg, status);
  }
  status = nc_inq_varid(ncid, name, &varid);
  if(status != NC_NOERR) {
    sprintf(msg, "in getting varid of %s from file %s.", name, file);
    handle_netcdf_error(msg, status);
  }

  status = nc_inq_vartype(ncid, varid, &vartype);
  if(status != NC_NOERR) {
    sprintf(msg, "get_var_data_region: in getting vartype of of %s in file %s ", name, file);
    handle_netcdf_error(msg, status);
  }

  switch (vartype) {
  case NC_DOUBLE:case NC_FLOAT:
#ifdef OVERLOAD_R4
    status = nc_get_vara_float(ncid, varid, start, nread, data);
#else      
    status = nc_get_vara_double(ncid, varid, start, nread, data);
#endif
    break;
  case NC_INT:
    status = nc_get_vara_int(ncid, varid, start, nread, data);
    break;
  default:
    sprintf(msg, "get_var_data_region: field %s in file %s has an invalid type, "
            "the type should be NC_DOUBLE, NC_FLOAT or NC_INT", name, file);
    error_handler(msg);
  }

  if(status != NC_NOERR) {
    sprintf(msg, "get_var_data_region: in getting data of %s from file %s.", name, file);
    handle_netcdf_error(msg, status);
  }
  status = nc_close(ncid);
  if(status != NC_NOERR) {
    sprintf(msg, "get_var_data_region: in closing file %s.", file);
    handle_netcdf_error(msg, status);
  }  
#else
  error_handler("read_mosaic: Add flag -Duse_netCDF when compiling");
#endif
  
}; /* get_var_data_region */

/******************************************************************************
   void get_var_text_att(const char *file, const char *name, const char *attname, char *att)
   get text attribute of field 'name' from 'file
******************************************************************************/
void get_var_text_att(const char *file, const char *name, const char *attname, char *att)
{
  int ncid, varid, status;  
  char msg[512];

#ifdef use_netCDF    
  status = nc_open(file, NC_NOWRITE, &ncid);
  if(status != NC_NOERR) {
    sprintf(msg, "in opening file %s", file);
    handle_netcdf_error(msg, status);
  }
  status = nc_inq_varid(ncid, name, &varid);
  if(status != NC_NOERR) {
    sprintf(msg, "in getting varid of %s from file %s.", name, file);
    handle_netcdf_error(msg, status);
  }     
  status = nc_get_att_text(ncid, varid, attname, att);
  if(status != NC_NOERR) {
    sprintf(msg, "in getting attribute %s of %s from file %s.", attname, name, file);
    handle_netcdf_error(msg, status);
  }
  status = nc_close(ncid);
  if(status != NC_NOERR) {
    sprintf(msg, "in closing file %s.", file);
    handle_netcdf_error(msg, status);
  }  
#else
  error_handler("read_mosaic: Add flag -Duse_netCDF when compiling");
#endif
  
}; /* get_var_text_att */

/***********************************************************************
  return number of overlapping cells.
***********************************************************************/
#ifndef __AIX
int read_mosaic_xgrid_size_( const char *xgrid_file )
{
  return read_mosaic_xgrid_size(xgrid_file);
}
#endif

int read_mosaic_xgrid_size( const char *xgrid_file )
{
  int ncells;
  
  ncells = get_dimlen(xgrid_file, "ncells");
  return ncells;
}

#ifdef OVERLOAD_R4
float get_global_area(void) 
{
  float garea;
#else
double get_global_area(void) 
{
  double garea;
#endif
  garea = 4*M_PI*RADIUS*RADIUS;

  return garea;
};

#ifndef __AIX
#ifdef OVERLOAD_R4
float get_global_area_(void)
{
  float garea;
#else
double get_global_area_(void)
{
  double garea;
#endif
  garea = 4*M_PI*RADIUS*RADIUS;

  return garea;
};
#endif


/****************************************************************************/
#ifndef __AIX
#ifdef OVERLOAD_R4
void read_mosaic_xgrid_order1_(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, float *area )
#else
void read_mosaic_xgrid_order1_(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, double *area )
#endif
{
  read_mosaic_xgrid_order1(xgrid_file, i1, j1, i2, j2, area);
  
};
#endif

#ifdef OVERLOAD_R4
void read_mosaic_xgrid_order1(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, float *area )
#else
void read_mosaic_xgrid_order1(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, double *area )
#endif
{
  int    ncells, n;
  int    *tile1_cell, *tile2_cell;
#ifdef OVERLOAD_R4  
  float garea;
#else
  double garea;
#endif
  
  ncells = get_dimlen(xgrid_file, "ncells");

  tile1_cell       = (int *)malloc(ncells*2*sizeof(int));
  tile2_cell       = (int *)malloc(ncells*2*sizeof(int));
  get_var_data(xgrid_file, "tile1_cell", tile1_cell);
  get_var_data(xgrid_file, "tile2_cell", tile2_cell);

  get_var_data(xgrid_file, "xgrid_area", area);

  garea = 4*M_PI*RADIUS*RADIUS;
  
  for(n=0; n<ncells; n++) {
    i1[n] = tile1_cell[n*2] - 1;
    j1[n] = tile1_cell[n*2+1] - 1;
    i2[n] = tile2_cell[n*2] - 1;
    j2[n] = tile2_cell[n*2+1] - 1;
    area[n] /= garea; /* rescale the exchange grid area to unit earth area */
  }

  free(tile1_cell);
  free(tile2_cell);
  
}; /* read_mosaic_xgrid_order1 */


#ifndef __AIX
#ifdef OVERLOAD_R4
void read_mosaic_xgrid_order1_region_(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, float *area, int *isc, int *iec )
#else
void read_mosaic_xgrid_order1_region_(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, double *area, int *isc, int *iec )
#endif
{
  read_mosaic_xgrid_order1_region(xgrid_file, i1, j1, i2, j2, area, isc, iec);
  
};
#endif

#ifdef OVERLOAD_R4
void read_mosaic_xgrid_order1_region(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, float *area, int *isc, int *iec )
#else
void read_mosaic_xgrid_order1_region(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, double *area, int *isc, int *iec )
#endif
{
  int    ncells, n, i;
  int    *tile1_cell, *tile2_cell;
  size_t start[4], nread[4];
#ifdef OVERLOAD_R4  
  float garea;
#else
  double garea;
#endif
  
  ncells = *iec-*isc+1;

  tile1_cell       = (int *)malloc(ncells*2*sizeof(int));
  tile2_cell       = (int *)malloc(ncells*2*sizeof(int));
  for(i=0; i<4; i++) {
    start[i] = 0; nread[i] = 1;
  }

  start[0] = *isc;
  nread[0] = ncells;
  nread[1] = 2;

  get_var_data_region(xgrid_file, "tile1_cell", start, nread, tile1_cell);
  get_var_data_region(xgrid_file, "tile2_cell", start, nread, tile2_cell);

  nread[1] = 1;
  
  get_var_data_region(xgrid_file, "xgrid_area", start, nread, area);

  garea = 4*M_PI*RADIUS*RADIUS;
  
  for(n=0; n<ncells; n++) {
    i1[n] = tile1_cell[n*2] - 1;
    j1[n] = tile1_cell[n*2+1] - 1;
    i2[n] = tile2_cell[n*2] - 1;
    j2[n] = tile2_cell[n*2+1] - 1;
    area[n] /= garea; /* rescale the exchange grid area to unit earth area */
  }

  free(tile1_cell);
  free(tile2_cell);
  
}; /* read_mosaic_xgrid_order1 */

/* NOTE: di, dj is for tile1, */
/****************************************************************************/
#ifndef __AIX
#ifdef OVERLOAD_R4
void read_mosaic_xgrid_order2_(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, float *area, float *di, float *dj )
#else
void read_mosaic_xgrid_order2_(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, double *area, double *di, double *dj )
#endif
{
  read_mosaic_xgrid_order2(xgrid_file, i1, j1, i2, j2, area, di, dj);
  
};
#endif
#ifdef OVERLOAD_R4
void read_mosaic_xgrid_order2(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, float *area, float *di, float *dj )
#else
void read_mosaic_xgrid_order2(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, double *area, double *di, double *dj )
#endif

{
  int    ncells, n;
  int    *tile1_cell, *tile2_cell;
  double *tile1_distance;
#ifdef OVERLOAD_R4
  float garea;
#else
  double garea;
#endif
  ncells = get_dimlen(xgrid_file, "ncells");

  tile1_cell       = (int    *)malloc(ncells*2*sizeof(int   ));
  tile2_cell       = (int    *)malloc(ncells*2*sizeof(int   ));
  tile1_distance   = (double *)malloc(ncells*2*sizeof(double));
  get_var_data(xgrid_file, "tile1_cell", tile1_cell);
  get_var_data(xgrid_file, "tile2_cell", tile2_cell);
  get_var_data(xgrid_file, "xgrid_area", area);
  get_var_data(xgrid_file, "tile1_distance", tile1_distance);

  garea = 4*M_PI*RADIUS*RADIUS;
  
  for(n=0; n<ncells; n++) {
    i1[n] = tile1_cell[n*2] - 1;
    j1[n] = tile1_cell[n*2+1] - 1;
    i2[n] = tile2_cell[n*2] - 1;
    j2[n] = tile2_cell[n*2+1] - 1;
    di[n] = tile1_distance[n*2];
    dj[n] = tile1_distance[n*2+1];
    area[n] /= garea; /* rescale the exchange grid area to unit earth area */
  }

  free(tile1_cell);
  free(tile2_cell);
  free(tile1_distance);
  
}; /* read_mosaic_xgrid_order2 */

/******************************************************************************
  int read_mosaic_ntiles(const char *mosaic_file)
  return number tiles in mosaic_file
******************************************************************************/
#ifndef __AIX
int read_mosaic_ntiles_(const char *mosaic_file)
{
  return read_mosaic_ntiles(mosaic_file);
}
#endif
int read_mosaic_ntiles(const char *mosaic_file)
{

  int ntiles;

  ntiles = get_dimlen(mosaic_file, "ntiles");

  return ntiles;
  
}; /* read_mosaic_ntiles */

/******************************************************************************
  int read_mosaic_ncontacts(const char *mosaic_file)
  return number of contacts in mosaic_file
******************************************************************************/
#ifndef __AIX
int read_mosaic_ncontacts_(const char *mosaic_file)
{
  return read_mosaic_ncontacts(mosaic_file);
}
#endif
int read_mosaic_ncontacts(const char *mosaic_file)
{

  int ncontacts;

  if(field_exist(mosaic_file, "contacts") )
    ncontacts = get_dimlen(mosaic_file, "ncontact");
  else
    ncontacts = 0;
  
  return ncontacts;
  
}; /* read_mosaic_ncontacts */


/*****************************************************************************
  void read_mosaic_grid_sizes(const char *mosaic_file, int *nx, int *ny)
  read mosaic grid size of each tile, currently we are assuming the refinement is 2.
  We assume the grid files are located at the same directory as mosaic_file.
*****************************************************************************/
#ifndef __AIX
void read_mosaic_grid_sizes_(const char *mosaic_file, int *nx, int *ny)
{
  read_mosaic_grid_sizes(mosaic_file, nx, ny);
}
#endif
void read_mosaic_grid_sizes(const char *mosaic_file, int *nx, int *ny)
{
  int ntiles, n;
  char gridfile[STRING], tilefile[2*STRING];
  char dir[STRING];
  const int x_refine = 2, y_refine = 2;

  get_file_dir(mosaic_file, dir);  
  ntiles = get_dimlen(mosaic_file, "ntiles");
  for(n = 0; n < ntiles; n++) {
    get_string_data_level(mosaic_file, "gridfiles", gridfile, &n);
    sprintf(tilefile, "%s/%s", dir, gridfile);
    nx[n] = get_dimlen(tilefile, "nx");
    ny[n] = get_dimlen(tilefile, "ny");
    if(nx[n]%x_refine != 0) error_handler("Error from read_mosaic_grid_sizes: nx is not divided by x_refine");
    if(ny[n]%y_refine != 0) error_handler("Error from read_mosaic_grid_sizes: ny is not divided by y_refine");
    nx[n] /= x_refine;
    ny[n] /= y_refine;
  }
  
}; /* read_mosaic_grid_sizes */
  

/******************************************************************************
  void read_mosaic_contact(const char *mosaic_file)
  read mosaic contact information
******************************************************************************/
#ifndef __AIX
void read_mosaic_contact_(const char *mosaic_file, int *tile1, int *tile2, int *istart1, int *iend1,
			 int *jstart1, int *jend1, int *istart2, int *iend2, int *jstart2, int *jend2)
{
  read_mosaic_contact(mosaic_file, tile1, tile2, istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2);
}
#endif

/* transfer the index from supergrid to model grid 
   return 0 if istart = iend
   otherwise return 1.
*/

int transfer_to_model_index(int istart_in, int iend_in, int *istart_out, int *iend_out, int refine_ratio)
{

   int type;

   if( istart_in == iend_in ) {
      type = 0;
      istart_out[0] = (istart_in+1)/refine_ratio-1;
      iend_out[0]  = istart_out[0];
   }   
   else {
      type = 1;
      if( iend_in > istart_in ) {
        istart_out[0] = istart_in - 1;
        iend_out[0]   = iend_in - refine_ratio;
      }
      else {
        istart_out[0] = istart_in - refine_ratio;
        iend_out[0]   = iend_in - 1;
      }

      if( istart_out[0]%refine_ratio || iend_out[0]%refine_ratio)
         error_handler("Error from read_mosaic: mismatch between refine_ratio and istart_in/iend_in");
      istart_out[0] /= refine_ratio;
      iend_out[0]   /= refine_ratio;
   }      

   return type;

}


void read_mosaic_contact(const char *mosaic_file, int *tile1, int *tile2, int *istart1, int *iend1,
			 int *jstart1, int *jend1, int *istart2, int *iend2, int *jstart2, int *jend2)
{
  char contacts[STRING];
  char **gridtiles;
#define MAXVAR 40
  char pstring[MAXVAR][STRING];
  int ntiles, ncontacts, n, m, l, found;
  unsigned int nstr;
  const int x_refine = 2, y_refine = 2;
  int i1_type, j1_type, i2_type, j2_type;  

  ntiles = get_dimlen(mosaic_file, "ntiles");
  gridtiles = (char **)malloc(ntiles*sizeof(char *));
  for(n=0; n<ntiles; n++) {
    gridtiles[n] = (char *)malloc(STRING*sizeof(char));
    get_string_data_level(mosaic_file, "gridtiles", gridtiles[n], &n);
  }
    
  ncontacts = get_dimlen(mosaic_file, "ncontact"); 
  for(n = 0; n < ncontacts; n++) {
    get_string_data_level(mosaic_file, "contacts", contacts, &n);
    /* parse the string contacts to get tile number */
    tokenize( contacts, ":", STRING, MAXVAR, (char *)pstring, &nstr);
    if(nstr != 4) error_handler("Error from read_mosaic: number of elements "
				 "in contact seperated by :/:: should be 4");
    found = 0;
    for(m=0; m<ntiles; m++) {
      if(strcmp(gridtiles[m], pstring[1]) == 0) { /*found the tile name */
	found = 1;
	tile1[n] = m+1;
	break;
      }
    }
    if(!found) error_handler("error from read_mosaic: the first tile name specified "
			     "in contact is not found in tile list");
    found = 0;
    for(m=0; m<ntiles; m++) {
      if(strcmp(gridtiles[m], pstring[3]) == 0) { /*found the tile name */
	found = 1;
	tile2[n] = m+1;
	break;
      }
    }
    if(!found) error_handler("error from read_mosaic: the second tile name specified "
			     "in contact is not found in tile list");    
    get_string_data_level(mosaic_file, "contact_index", contacts, &n);
    /* parse the string to get contact index */
    tokenize( contacts, ":,", STRING, MAXVAR, (char *)pstring, &nstr);
    if(nstr != 8) error_handler("Error from read_mosaic: number of elements "
				 "in contact_index seperated by :/, should be 8");
    /* make sure the string is only composed of numbers */
    for(m=0; m<nstr; m++) for(l=0; l<strlen(pstring[m]); l++) {
      if(pstring[m][l] > '9' ||  pstring[m][l] < '0' ) {
	error_handler("Error from read_mosaic: some of the character in "
		      "contact_indices except token is not digit number");
      }
    }
    istart1[n] = atoi(pstring[0]);
    iend1[n]   = atoi(pstring[1]);
    jstart1[n] = atoi(pstring[2]);
    jend1[n]   = atoi(pstring[3]);
    istart2[n] = atoi(pstring[4]);
    iend2[n]   = atoi(pstring[5]);
    jstart2[n] = atoi(pstring[6]);
    jend2[n]   = atoi(pstring[7]);
    i1_type = transfer_to_model_index(istart1[n], iend1[n], istart1+n, iend1+n, x_refine);
    j1_type = transfer_to_model_index(jstart1[n], jend1[n], jstart1+n, jend1+n, y_refine);
    i2_type = transfer_to_model_index(istart2[n], iend2[n], istart2+n, iend2+n, x_refine);
    j2_type = transfer_to_model_index(jstart2[n], jend2[n], jstart2+n, jend2+n, y_refine);
    if( i1_type == 0 && j1_type == 0 )
      error_handler("Error from read_mosaic_contact:istart1==iend1 and jstart1==jend1");
    if( i2_type == 0 && j2_type == 0 )
      error_handler("Error from read_mosaic_contact:istart2==iend2 and jstart2==jend2");
    if( i1_type + j1_type != i2_type + j2_type )
      error_handler("Error from read_mosaic_contact: It is not a line or overlap contact");

  }

  for(m=0; m<ntiles; m++) {
    free(gridtiles[m]);
  }

  free(gridtiles);
 

}; /* read_mosaic_contact */


/******************************************************************************
  void read_mosaic_grid_data(const char *mosaic_file, const char *name, int nx, int ny,
                             double *data, int level, int ioff, int joff)
  read mosaic grid information onto model grid. We assume the refinement is 2 right now.
  We may remove this restriction in the future. nx and ny are model grid size. level
  is the tile number. ioff and joff to indicate grid location. ioff =0 and joff = 0
  for C-cell. ioff=0 and joff=1 for E-cell, ioff=1 and joff=0 for N-cell,
  ioff=1 and joff=1 for T-cell
******************************************************************************/
void read_mosaic_grid_data(const char *mosaic_file, const char *name, int nx, int ny,
                           double *data, int level, int ioff, int joff)
{
  char   tilefile[STRING], gridfile[STRING], dir[STRING];
  double *tmp;
  int    ni, nj, nxp, nyp, i, j;

  get_file_dir(mosaic_file, dir);
  
  get_string_data_level(mosaic_file, "gridfiles", gridfile, &level);
  sprintf(tilefile, "%s/%s", dir, gridfile);
  
  ni = get_dimlen(tilefile, "nx");
  nj = get_dimlen(tilefile, "ny");

  if( ni != nx*2 || nj != ny*2) error_handler("supergrid size should be double of the model grid size");
  tmp = (double *)malloc((ni+1)*(nj+1)*sizeof(double));
  get_var_data( tilefile, name, tmp);
  nxp = nx + 1 - ioff;
  nyp = ny + 1 - joff;
  for(j=0; j<nyp; j++) for(i=0; i<nxp; i++) data[j*nxp+i] = tmp[(2*j+joff)*(ni+1)+2*i+ioff];
  free(tmp);
   
}; /* read_mosaic_grid_data */


