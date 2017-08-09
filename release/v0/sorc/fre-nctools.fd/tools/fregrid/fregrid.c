/*
  This program remaps (scalar or vector) data from the input grid
  to the output grid

 AUTHOR: Zhi Liang (Zhi.Liang@noaa.gov)
          NOAA Geophysical Fluid Dynamics Lab, Princeton, NJ
 
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  For the full text of the GNU General Public License,
  write to: Free Software Foundation, Inc.,
            675 Mass Ave, Cambridge, MA 02139, USA.  
-----------------------------------------------------------------------
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include "globals.h"
#include "constant.h"
#include "read_mosaic.h"
#include "mpp_io.h"
#include "mpp.h"
#include "mosaic_util.h"
#include "conserve_interp.h"
#include "bilinear_interp.h"
#include "fregrid_util.h"

char *usage[] = {
  "",
  "  fregrid --input_mosaic input_mosaic --input_file input_file                         ",
  "          [--scalar_field scalar_fld] [--u_field u_fld]  [--v_field v_fld]            ",
  "          [--output_mosaic output_mosaic] [--lonBegin #decimal] [--lonEnd #decimal]   ",
  "          [--latBegin #decimal] [--latEnd #decimal] [--nlon #integer]                 ",
  "          [--nlat #integer] [--KlevelBegin #integer] [--KlevelEnd #integer]           ",
  "          [--LstepBegin #integer] [--LstepEnd #integer]                               ",
  "          [--output_file output_file] [--input_dir input_dir]                         ",
  "          [--output_dir output_dir] [--remap_file remap_file]                         ",
  "          [--interp_method method] [--grid_type grid_type] [--test_case test_case]    ",
  "          [--symmetry] [--target_grid] [--finer_step #] [--fill_missing]              ",
  "          [--center_y] [--check_conserve] [--weight_file weight_file]                 ",
  "          [--weight_field --weight_field] [--dst_vgrid dst_vgrid]                     ",
  "          [--extrapolate] [--stop_crit #] [--standard_dimension]                      ",
  "          [--associated_file_dir dir]                                                 ",
  "                                                                                      ",
  "fregrid remaps data (scalar or vector) from input_mosaic onto                         ",
  "output_mosaic.  Note that the target grid also could be specified                     ",
  "through lonBegin, lonEnd, latBegin, latEnd, nlon and nlat. Currently                  ",
  "only T-cell scalar regridding and AGRID vector regridding is                          ",
  "available.  Bilinear interpolation is implemented only for cubic grid                 ",
  "vector interpolation. The interpolation algorithm used is controlled                  ",
  "by --interp_method with default 'conserve_order1'. Currently the                      ",
  "'conserve_order1', 'conserve_order2' and 'bilinear' remapping schemes                 ",
  "are implemented. 'bilinear' is only used to remap data from cubic grid                ",
  "to latlon grid. Alternative schemes can be added if needed. fregrid                   ",
  "expects NetCDF format input. scalar_field and/or u_field/v_field must                 ",
  "be specified. u_fld and v_fld must be paired together.                                ",
  "                                                                                      ",
  "fregrid takes the following flags:                                                    ",
  "                                                                                      ",
  "REQUIRED:                                                                             ",
  "                                                                                      ",
  "--input_mosaic  input_mosaic  specify the input mosaic information. This file         ",
  "                              contains list of tile files which specify the grid      ",
  "                              information for each tile.                              ",
  "                                                                                      ",
  "OPTIONAL FLAGS                                                                        ",
  "                                                                                      ",  
  "--input_file    input_file    specify the input file name. The suffix '.nc' can be    ",
  "                              omitted. The suffix 'tile#' should not present for      ",
  "                              multiple-tile files. The number of files must be 1 for  ",
  "                              scalar regridding and can be 1 or 2 for vector          ",
  "                              regridding. File path should not be includes.           ",
  "                                                                                      ",
  "--scalar_field    scalar_fld  specify the scalar field name to be regridded. The      ",
  "                              multiple entry field names are seperated by comma.      ",
  "                                                                                      ",    
  "--u_field         u_fld       specify the vector field u-componentname to be          ",
  "                              regridded. The multiple entry field names are seperated ",
  "                              by comma. u_field must be paired together with v_field. ",
  "                                                                                      ",    
  "--v_field         v_fld       specify the vector field v-componentname to be          ",
  "                              regridded. The multiple entry field names are seperated ",
  "                              by comma. v_field must be paired together with u_field. ",
  "                                                                                      ",  
  "--output_mosaic output_mosaic specify the output mosaic information. This file        ",
  "                              contains list of tile files which specify the grid      ",
  "                              information for each tile. If output_mosaic is not      ",
  "                              specified, nlon and nlat must be specified.             ",
  "                                                                                      ",
  "--lonBegin  #decimal          specify the starting longitude(in degree) of the        ",
  "                              geographical region of the target grid on which the     ",
  "                              output is desired. The default value is 0.              ",
  "                                                                                      ",
  "--lonEnd   #decimal           specify the ending longitude(in degree) of the          ",
  "                              geographical region of the target grid on which the     ",
  "                              output is desired. The default value is 360.            ",
  "                                                                                      ",
  "--latBegin  #decimal          specify the starting latitude(in degree) of the         ",
  "                              geographical region of the target grid on which the     ",
  "                              output is desired. The default value is -90.            ",
  "                                                                                      ",
  "--latEnd   #decimal           specify the ending latitude(in degree) of the           ",
  "                              geographical region of the target grid on which the     ",
  "                              output is desired. The default value is 90.             ",  
  "                                                                                      ",
  "--nlon #integer               specify number of grid box cells in x-direction for a   ",
  "                              regular lat-lon grid.                                   ",
  "                                                                                      ",
  "--nlat #integer               specify number of grid box cells in y-direction for a   ",
  "                              regular lat-lon grid.                                   ",
  "                                                                                      ",  
  "--KlevelBegin #integer        specify begin index of the k-level (depth axis) that    ",
  "                              to be regridded.                                        ",
  "                                                                                      ",
  "--KlevelEnd #integer          specify end index of the k-level (depth axis) that      ",
  "                              to be regridded.                                        ",  
  "                                                                                      ",
  "--LstepBegin #integer         specify the begin index of L-step (time axis) that      ",
  "                              to be regridded.                                        ",
  "                                                                                      ",
  "--LstepEnd #integer           specify the end index of L-step (time axis) that        ",
  "                              to be regridded.                                        ",  
  "                                                                                      ",  
  "--output_file   output_file   specify the output file name. If not presented,         ",
  "                              output_file will take the value of input_file. The      ",
  "                              suffix '.nc' can be omitted. The suffix 'tile#' should  ",
  "                              not present for multiple-tile files. The number of      ",
  "                              files must be 1 for scalar regridding and can be 1 or 2 ",
  "                              for vector regridding. File path should not be includes.",  
  "                                                                                      ",
  "--input_dir     input_dir     specify the path that stores input_file. If not         ",
  "                              presented, the input file is assumed to be stored in    ",
  "                              current diretory.                                       ",
  "                                                                                      ",
  "--output_dir   output_dir     specify the path that will store output file. If not    ",
  "                              presented, the output file will be stored in current    ",
  "                              diretory.                                               ",
  "                                                                                      ",
  "--remap_file   remap_file     specify the file name that saves remapping information. ",
  "                              If remap_file is specified and the file does not exist, ",
  "                              remapping information will be calculated ans stored in  ",
  "                              remap_file. If remap_file is specified and the file     ",
  "                              exists, remapping information will be read from         ",
  "                              remap_file.                                             ",
  "                                                                                      ",
  "--interp_method interp_method specify the remapping algorithm to be used. Default is  ",
  "                              'conserve_order1'. Currently only 'conserve_order1',    ",
  "                              'conserve_order2', 'conserve_order2_monotonic' and      ",
  "                              'bilinear' remapping scheme are   ",
  "                              implemented in this tool. The bilinear scheme can only  ",
  "                              be used to remap data from cubic grid to regular latlon ",
  "                              grid. When interp_method is 'bilinear', nlon and nlat   ",
  "                              must be specified and the output data in y-direction    ",
  "                              will be located at the center of cell or bound of the   ",
  "                              cell depending on the setting of y_center.              ",
  "                                                                                      ",
  "--test_case test_case         specify the test function to be used for testing.       ",
  "                                                                                      ",
  "--grid_type     grid_type     specify the vector field grid location. default is      ",
  "                              AGRID and only AGRID is implemented yet.                ",
  "                                                                                      ",
  "--symmetry                    indicate the grid is symmetry or not.                   ",
  "                                                                                      ",
  "--target_grid                 use taget grid cell area instead of calculating based on",
  "                              exchange grid area. default is off.                     ",
  "                                                                                      ",
  "---finer_step #integer        This is used only for bilinear interpolation. Set       ",
  "                              finer_step to a positive integer to reduce noise in     ",
  "                              interpolation and get a relatively smooth output. The   ",
  "                              default value is 0. When finer_step is greater than 0,  ",
  "                              fregrid will first remap data from source grid onto a   ",
  "                              finer grid with resolution that is power of 2 of        ",
  "                              destination grid resolution using bilinear              ",
  "                              interpolation, then using volume averaging to remap     ",
  "                              data from finer grid onto destination grid.             ",
  "                                                                                      ",
  "--center_y                    output latitude will locate at cell center, i.e., the   ",
  "                              starting latitude will be -89 when nlat = 90. when      ",
  "                              center_y is not set, starting latitude will be -90. for ",
  "                              bilinear interpolation. For conservative interpolation, ",
  "                              center_y is assumed.                                    ",
  "                                                                                      ",
  "--check_conserve              check the conservation of conservative interpolation.   ",
  "                              The area sum will be printed out for input and output   ",
  "                              mosaic.                                                 ",
  "                                                                                      ",
  "--monotonic                   When specified, use monotonic interpolation when        ",
  "                              interp_method is 'conserve_order2'.                     ",
  "                                                                                      ",
  "--weight_file                 Specify the filename that store weight_field. The       ",
  "                              suffix '.tile#.nc' should not present for multiple-tile ",
  "                              files. weight_field is used to adjust the source weight.",
  "                              Normally it could be area fraction. When weight_field   ",
  "                              is specified, the weight_file will default to be        ",
  "                              input_file if weight_file is not specified.             ",
  "                                                                                      ",
  "--weight_field                Specify the name of weight field in weight_file         ",
  "                                                                                      ",
  "--dst_vgrid                   specify the destination vertical grid file. Data will   ",
  "                              be remapped onto the destination vertical grid.         ",
  "                              When --dst_vgrid is specified, --extrapolate is         ",
  "                              assumed to be specified.                                ",
  "                                                                                      ",
  "--extrapolate                 Will extrapolate data onto masked points when specified.",
  "                                                                                      ",
  "--stop_crit #                 The stopping criteria when extrapping data onto missing ",
  "                              points. Default is 0.005                                ",
  "                                                                                      ",
  "--standard_dimension          When specified, the dimension and field name for        ",
  "                              longitude and latitude axis will be 'lon' and 'lat'.    ",
  "                              'lon_bnd' and 'lat_bnd' will be longitude and latitude  ",
  "                              bound name. The dimension of lon_bounds is (2,nlon) and ",
  "                              the dimension of lat_bounds is (2,nlat).                "
  "                                                                                      ",
  "--associated_file_dir dir     Specify the path of the associated files                ",
  "                                                                                      ",
  "--debug                       Will print out memory usage and running time            ",
  "                                                                                      ",
  "--nthreads #                  Specify number of OpenMP threads.                       ",
  "                                                                                      ",
  "  Example 1: Remap C48 data onto N45 grid.                                            ",          
  "             (use GFDL-CM3 data as example)                                           ",
  "   fregrid --input_mosaic C48_mosaic.nc --input_dir input_dir --input_file input_file ",
  "           --scalar_field temp,salt --nlon 144 --nlat 90                              ",
  "                                                                                      ",
  "  Example 2: Remap data onto cm2m ocean grid with extrapolation                       ",
  "             and vertical interpolation.                                              ",
  "   fregrid --input_mosaic levitus_mosaic.nc --input_dir inputdir                      ",
  "           --input_file WOA09_ann_theta.nc --scalar_field POTENTIAL_TEMP              ",
  "           --output_file WOA09_ann_theta_cm2g_extrap.nc                               ",
  "           --output_mosaic cm2m_ocean_mosaic.nc --extrapolate                         ",
  "           --dst_vgrid inputdir/cm2m_ocean_vgrid.nc                                   ",
  "                                                                                      ",
  NULL};
#define EPSLN10  (1.e-10)
const double D2R = M_PI/180.;
char tagname[] = "$Name: bronx-10_performance_z1l $";

int main(int argc, char* argv[])
{
  unsigned int opcode = 0;
  char    *mosaic_in=NULL;            /* input mosaic file name */
  char    *mosaic_out=NULL;           /* input mosaic file name */
  char    *dir_in=NULL;               /* input file location */
  char    *dir_out=NULL;              /* output file location */
  int     ntiles_in = 0;              /* number of tiles in input mosaic */
  int     ntiles_out = 0;             /* number of tiles in output mosaic */
  unsigned int nfiles     = 0;        /* number of input file */
  unsigned int nfiles_out = 0;             /* number of output file */
  char    input_file [NFILE][STRING];
  char    output_file[NFILE][STRING];
  char    scalar_name[NVAR] [STRING];
  char    scalar_name_remap[NVAR][STRING];
  char    u_name     [NVAR] [STRING];
  char    v_name     [NVAR] [STRING];
  char    *test_case = NULL;
  double  test_param = 1;
  char    *associated_file_dir = NULL;
  int     check_conserve = 0; /* 0 means no check */
  double  lonbegin = 0, lonend = 360;
  double  latbegin = -90, latend = 90;			  
  int     nlon = 0, nlat = 0;
  int     kbegin = 0, kend = -1; 
  int     lbegin = 0, lend = -1;
  char    *remap_file = NULL;
  char    interp_method[STRING] = "conserve_order1";
  int     y_at_center = 0;
  int     grid_type = AGRID;
  unsigned int nscalar=0, nvector=0, nvector2=0;
  unsigned int nscalar_orig;
  int     option_index, c, i, n, m, l;
  char    entry[MAXSTRING];  /* should be long enough */
  char    txt[STRING];
  char    history[MAXATT];
  int     fill_missing = 0;
  int     extrapolate = 0;
  int     vertical_interp = 0;
  char    *dst_vgrid = NULL;
  double  stop_crit=0.005;
  unsigned int  finer_step = 0;
  int     debug = 0;
  int     great_circle_algorithm_in, great_circle_algorithm_out;
  
  char          wt_file_obj[512];
  char          *weight_file=NULL;
  char          *weight_field = NULL;

  VGrid_config   vgrid_in;            /* store vertical input grid */
  VGrid_config   vgrid_out;           /* store vertical output grid */
  Grid_config   *grid_in    = NULL;   /* store input grid  */
  Grid_config   *grid_out   = NULL;   /* store output grid */
  Field_config  *scalar_in  = NULL;   /* store input scalar data */
  Field_config  *scalar_out = NULL;   /* store output scalar data */
  Field_config  *u_in       = NULL;   /* store input vector u-component */
  Field_config  *v_in       = NULL;   /* store input vector v-component */
  Field_config  *u_out      = NULL;   /* store input vector u-component */
  Field_config  *v_out      = NULL;   /* store input vector v-component */
  File_config   *file_in    = NULL;   /* store input file information */
  File_config   *file_out   = NULL;   /* store output file information */
  File_config   *file2_in   = NULL;   /* store input file information */
  File_config   *file2_out  = NULL;   /* store output file information */
  Bound_config  *bound_T    = NULL;   /* store halo update information for T-cell*/
  Interp_config *interp     = NULL;   /* store remapping information */
  int save_weight_only      = 0;
  int nthreads = 1;

  double time_get_in_grid=0, time_get_out_grid=0, time_get_input=0;
  double time_setup_interp=0, time_do_interp=0, time_write=0;
  clock_t time_start, time_end;
  
  int errflg = (argc == 1);
  int fid;
  
  static struct option long_options[] = {
    {"input_mosaic",     required_argument, NULL, 'a'},
    {"output_mosaic",    required_argument, NULL, 'b'},
    {"input_dir",        required_argument, NULL, 'c'},
    {"output_dir",       required_argument, NULL, 'd'},
    {"input_file",       required_argument, NULL, 'e'},
    {"output_file",      required_argument, NULL, 'f'},
    {"remap_file",       required_argument, NULL, 'g'},
    {"test_case",        required_argument, NULL, 'i'},
    {"interp_method",    required_argument, NULL, 'j'},
    {"test_parameter",   required_argument, NULL, 'k'},
    {"symmetry",         no_argument,       NULL, 'l'},
    {"grid_type",        required_argument, NULL, 'm'},
    {"target_grid",      no_argument,       NULL, 'n'},
    {"finer_step",       required_argument, NULL, 'o'},
    {"fill_missing",     no_argument,       NULL, 'p'},
    {"nlon",             required_argument, NULL, 'q'},
    {"nlat",             required_argument, NULL, 'r'},
    {"scalar_field",     required_argument, NULL, 's'},
    {"check_conserve",   no_argument,       NULL, 't'},
    {"u_field",          required_argument, NULL, 'u'},
    {"v_field",          required_argument, NULL, 'v'},
    {"center_y",         no_argument,       NULL, 'y'},
    {"lonBegin",         required_argument, NULL, 'A'},
    {"lonEnd",           required_argument, NULL, 'B'},
    {"latBegin",         required_argument, NULL, 'C'},
    {"latEnd",           required_argument, NULL, 'D'},
    {"KlevelBegin",      required_argument, NULL, 'E'},
    {"KlevelEnd",        required_argument, NULL, 'F'},
    {"LstepBegin",       required_argument, NULL, 'G'},
    {"LstepEnd",         required_argument, NULL, 'H'},
    {"weight_file",      required_argument, NULL, 'I'},
    {"weight_field",     required_argument, NULL, 'J'},
    {"extrapolate",      no_argument,       NULL, 'L'},
    {"dst_vgrid",        required_argument, NULL, 'M'},
    {"stop_crit",        required_argument, NULL, 'N'},
    {"standard_dimension", no_argument,     NULL, 'O'},
    {"debug",             no_argument,     NULL, 'P'},
    {"nthreads",         required_argument, NULL, 'Q'},
    {"associated_file_dir", required_argument, NULL, 'R'},
    {"help",             no_argument,       NULL, 'h'},
    {0, 0, 0, 0},
  };  
  
  /* start parallel */
  mpp_init(&argc, &argv);
  mpp_domain_init();
  
  while ((c = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
    switch (c) {
    case 'a':
      mosaic_in  = optarg;
      break;
    case 'b':
      mosaic_out = optarg;
      break;
    case 'c':
      dir_in = optarg;
      break;
    case 'd':
      dir_out = optarg;
      break;
    case 'e':
      if(strlen(optarg) >= MAXSTRING) mpp_error("fregrid: the entry is not long for option -e");
      strcpy(entry, optarg);
      tokenize(entry, ",", STRING, NFILE, (char *)input_file, &nfiles);
      break;
    case 'f':
      if(strlen(optarg) >= MAXSTRING)  mpp_error("fregrid: the entry is not long for option -f");      
      strcpy(entry, optarg);
      tokenize(entry, ",", STRING, NFILE, (char *)output_file, &nfiles_out);
      break;
    case 'g':
      remap_file = optarg;
      break;
    case 's':
      if(strlen(optarg) >= MAXSTRING) mpp_error("fregrid: the entry is not long for option -s");      
      strcpy(entry, optarg);
      tokenize(entry, ",", STRING, NVAR, (char *)scalar_name, &nscalar);
      break;
    case 'u':
      if(strlen(optarg) >= MAXSTRING) mpp_error("fregrid: the entry is not long for option -u");      
      strcpy(entry, optarg);
      tokenize(entry, ",", STRING, NVAR, (char *)u_name, &nvector);
      break;        
    case 'v':
      if(strlen(optarg) >= MAXSTRING) mpp_error("fregrid: the entry is not long for option -v");      
      strcpy(entry, optarg);
      tokenize(entry, ",", STRING, NVAR, (char *)v_name, &nvector2);
      break;      
    case 'j':
      strcpy(interp_method, optarg);
      break;
    case 'i':
      test_case = optarg;
      break;
    case 'k':
      test_param = atof(optarg);
      break;      
    case 'l':
      opcode |= SYMMETRY;
      break;
    case 'm':
      if(strcmp(optarg, "AGRID") == 0)
	grid_type = AGRID;
      else if(strcmp(optarg, "BGRID") == 0)
	grid_type = BGRID;
      else
	mpp_error("fregrid: only AGRID and BGRID vector regridding are implmented, contact developer");
      break;
    case 'n':
      opcode |= TARGET;
      break;
    case 'o':
      finer_step = atoi(optarg);
      break;
    case 'p':
      fill_missing = 1;
      break;
    case 'q':
      nlon = atoi(optarg);
      break;
    case 'r':
      nlat = atoi(optarg);
      break;
    case 't':
      check_conserve = 1;
      break;
    case 'y':
      y_at_center = 1;
      break;
    case 'A':
      lonbegin = atof(optarg);
      break;
    case 'B':
      lonend = atof(optarg);
      break;
    case 'C':
      latbegin = atof(optarg);
      break;
    case 'D':
      latend = atof(optarg);
      break;
    case 'E':
      kbegin = atoi(optarg);
      break;
    case 'F':
      kend = atoi(optarg);
      break;
    case 'G':
      lbegin = atoi(optarg);
      break;
    case 'H':
      lend = atoi(optarg);
      break;
    case 'I':
      weight_file = optarg;
      break;
    case 'J':
      weight_field = optarg;
      break;
    case 'L':
      extrapolate = 1;
      break;
    case 'M':
      dst_vgrid = optarg;
      vertical_interp = 1;
      break;
    case 'N':
      stop_crit = atof(optarg);
      break;
    case 'O':
      opcode |= STANDARD_DIMENSION;
      break;
    case 'P':
      debug = 1;
      break;  
    case 'Q':
      nthreads = atoi(optarg);
      break;
    case 'R':
      associated_file_dir = optarg;
      break;
    case '?':
      errflg++;
      break;
    }
  }

  if (errflg) {
    char **u = usage;
    while (*u) { fprintf(stderr, "%s\n", *u); u++; }
    exit(2);
  }      
  /* check the arguments */
  if( !mosaic_in  ) mpp_error("fregrid: input_mosaic is not specified");
  if( !mosaic_out ) {
    if(nlon == 0 || nlat ==0 ) mpp_error("fregrid: when output_mosaic is not specified, nlon and nlat should be specified");
    if(lonend <= lonbegin) mpp_error("fregrid: when output_mosaic is not specified, lonEnd should be larger than lonBegin");
    if(latend <= latbegin) mpp_error("fregrid: when output_mosaic is not specified, latEnd should be larger than latBegin");
  }
  else {
    if(nlon !=0 || nlat != 0) mpp_error("fregrid: when output_mosaic is specified, nlon and nlat should not be specified");
  }

  if(!strcmp(interp_method, "conserve_order1") ) {
    if(mpp_pe() == mpp_root_pe())printf("****fregrid: first order conservative scheme will be used for regridding.\n");
    opcode |= CONSERVE_ORDER1;
  }
  else if(!strcmp(interp_method, "conserve_order2") ) {
    if(mpp_pe() == mpp_root_pe())printf("****fregrid: second order conservative scheme will be used for regridding.\n");
    opcode |= CONSERVE_ORDER2;
  }
  else if(!strcmp(interp_method, "conserve_order2_monotonic") ) {
    if(mpp_pe() == mpp_root_pe())printf("****fregrid: second order monotonic conservative scheme will be used for regridding.\n");
    opcode |= CONSERVE_ORDER2;
    opcode |= MONOTONIC;
  }
  else if(!strcmp(interp_method, "bilinear") ) {
    if(mpp_pe() == mpp_root_pe())printf("****fregrid: bilinear remapping scheme will be used for regridding.\n");  
    opcode |= BILINEAR;
  }
  else
    mpp_error("fregrid: interp_method must be 'conserve_order1', 'conserve_order2', 'conserve_order2_monotonic'  or 'bilinear'");
      
  if( nfiles == 0) {
    if(nvector > 0 || nscalar > 0 || nvector2 > 0)
      mpp_error("fregrid: when --input_file is not specified, --scalar_field, --u_field and --v_field should also not be specified");
    if(!remap_file) mpp_error("fregrid: when --input_file is not specified, remap_file must be specified to save weight information");
    save_weight_only = 1;
    if(mpp_pe()==mpp_root_pe())printf("NOTE: No input file specified in this run, no data file will be regridded "
				      "and only weight information is calculated.\n");
  }
  else if( nfiles == 1 || nfiles ==2) {
    if( nvector != nvector2 ) mpp_error("fregrid: number of fields specified in u_field must be the same as specified in v_field");
    if( nscalar+nvector==0 ) mpp_error("fregrid: both scalar_field and vector_field are not specified");
    /* when nvector =2 and nscalar=0, nfiles can be 2 otherwise nfiles must be 1 */
    if( nscalar && nfiles != 1 )
      mpp_error("fregrid: when scalar_field is specified, number of files must be 1");
    if( nfiles_out == 0 ) {
      for(i=0; i<nfiles; i++) strcpy(output_file[i], input_file[i]);
    }
    else if (nfiles_out != nfiles )
      mpp_error("fregrid:number of input file is not equal to number of output file");
  }
  else
    mpp_error("fregrid: number of input file should be 1 or 2");

  if(kbegin != 0 || kend != -1) { /* at least one of kbegin and kend is set */
    if(kbegin < 1 || kend < kbegin) mpp_error("fregrid:KlevelBegin should be a positive integer and no larger "
					      "than KlevelEnd when you want pick certain klevel");
  }
  if(lbegin != 0 || lend != -1) { /* at least one of lbegin and lend is set */
     if(lbegin < 1 || lend < lbegin) mpp_error("fregrid:LstepBegin should be a positive integer and no larger "
					      "than LstepEnd when you want pick certain Lstep");
  }

  if(weight_field) {
    if(nvector >0) mpp_error("fregrid: weight_field should not be specified for vector interpolation, contact developer");
    if(!weight_file) {
      
      if(nfiles==0) mpp_error("fregrid: weight_field is specified, but both weight_file and input_file are not specified");
      if(dir_in)
	sprintf(wt_file_obj, "%s/%s", dir_in, input_file[0]);
      else
	sprintf(wt_file_obj, "./%s", input_file[0]);
      weight_file = wt_file_obj;
    }
  }

  if(nvector > 0) {
    opcode |= VECTOR;
    if(grid_type == AGRID)
      opcode |= AGRID;
    else if(grid_type == BGRID)
      opcode |= BGRID;
  }
  
  /* define history to be the history in the grid file */
  strcpy(history,argv[0]);

  for(i=1;i<argc;i++) {
    strcat(history, " ");
    if(strlen(argv[i]) > MAXENTRY) { /* limit the size of each entry, here we are assume the only entry that is longer than
					MAXENTRY= 256 is the option --scalar_field --u_field and v_field */
      if(strcmp(argv[i-1], "--scalar_field") && strcmp(argv[i-1], "--u_field") && strcmp(argv[i-1], "--v_field") )
	mpp_error("fregrid: the entry ( is not scalar_field, u_field, v_field ) is too long, need to increase parameter MAXENTRY");
      strcat(history, "(**please see the field list in this file**)" );
    }
    else
      strcat(history, argv[i]);
  }
  
{
  int base_cpu;

#if defined(_OPENMP)
  omp_set_num_threads(nthreads);
  base_cpu = get_cpu_affinity();
#pragma omp parallel
  set_cpu_affinity(base_cpu+omp_get_thread_num() );
#endif

}

  /* get the mosaic information of input and output mosaic*/
  fid = mpp_open(mosaic_in, MPP_READ);
  ntiles_in = mpp_get_dimlen(fid, "ntiles");
  mpp_close(fid);

  /* second order conservative interpolation is only avail for the cubic sphere input grid */
  if( ntiles_in != 6 && (opcode & CONSERVE_ORDER2) ) 
    mpp_error("fregrid: when the input grid is not cubic sphere grid, interp_method can not be conserve_order2");

  if(mosaic_out) {
    fid = mpp_open(mosaic_out, MPP_READ);
    ntiles_out = mpp_get_dimlen(fid, "ntiles");
    mpp_close(fid);
  }
  else
    ntiles_out = 1;

  if(test_case) {
    if(nfiles != 1) mpp_error("fregrid: when test_case is specified, nfiles should be 1");
    sprintf(output_file[0], "%s.%s.output", test_case, interp_method);
  }

  if(check_conserve) opcode |= CHECK_CONSERVE;

  if( opcode & STANDARD_DIMENSION ) printf("fregrid: --standard_dimension is set\n");
  
  if( opcode & BILINEAR ) {
    int ncontact;
    ncontact = read_mosaic_ncontacts(mosaic_in);
    if( nlon == 0 || nlat == 0) mpp_error("fregrid: when interp_method is bilinear, nlon and nlat should be specified");
    if(ntiles_in != 6) mpp_error("fregrid: when interp_method is bilinear, the input mosaic should be 6 tile cubic grid");
    if(ncontact !=12)  mpp_error("fregrid: when interp_method is bilinear, the input mosaic should be 12 contact cubic grid");
    if(mpp_npes() > 1) mpp_error("fregrid: parallel is not implemented for bilinear remapping");
  }
  else 
    y_at_center = 1;

  if(extrapolate) opcode |= EXTRAPOLATE;
  
  /* memory allocation for data structure */
  grid_in   = (Grid_config *)malloc(ntiles_in *sizeof(Grid_config));
  grid_out  = (Grid_config *)malloc(ntiles_out*sizeof(Grid_config));
  bound_T   = (Bound_config *)malloc(ntiles_in *sizeof(Bound_config));
  interp    = (Interp_config *)malloc(ntiles_out*sizeof(Interp_config));

  if(debug) {
    print_mem_usage("Before calling get_input_grid");
    time_start = clock();
  }
  get_input_grid( ntiles_in, grid_in, bound_T, mosaic_in, opcode, &great_circle_algorithm_in, save_weight_only );
  if(debug) {
    time_end = clock();
    time_get_in_grid = 1.0*(time_end - time_start)/CLOCKS_PER_SEC;
    print_mem_usage("After calling get_input_grid");
    time_start = clock();
  }
  if(mosaic_out) 
    get_output_grid_from_mosaic( ntiles_out, grid_out, mosaic_out, opcode, &great_circle_algorithm_out );
  else {
    great_circle_algorithm_out = 0;
    get_output_grid_by_size(ntiles_out, grid_out, lonbegin, lonend, latbegin, latend,
			    nlon, nlat, finer_step, y_at_center, opcode);
  }
  if(debug) {
    time_end = clock();
    time_get_out_grid = 1.0*(time_end - time_start)/CLOCKS_PER_SEC;
    print_mem_usage("After calling get_output_grid");
  }
  /* find out if great_circle algorithm is used in the input grid or output grid */
  
  if( great_circle_algorithm_in == 0 && great_circle_algorithm_out == 0 )
    opcode |= LEGACY_CLIP;
  else {
    opcode |= GREAT_CIRCLE;
    /* currently only first-order conservative is implemented */
    if( !(opcode & CONSERVE_ORDER1) )
      mpp_error("fregrid: when clip_method is 'conserve_great_circle', interp_methos need to be 'conserve_order1', contact developer");
  }

  /* get the grid cell_area */
  get_input_output_cell_area(ntiles_in, grid_in, ntiles_out, grid_out, opcode);
  if(debug) print_mem_usage("After get_input_output_cell_area");  
  /* currently extrapolate are limited to ntiles = 1. extrapolate are limited to lat-lon input grid */
  if( extrapolate ) {
    int i, j, ind0, ind1, ind2;
    if(test_case ) mpp_error("fregrid: extrapolate is limited to test_case is false");
    if(ntiles_in != 1) mpp_error("fregrid: extrapolate is limited to ntile_in = 1");
    /* check if the grid is lat-lon grid */
    for(j=1; j<=grid_in[0].ny; j++) for(i=1; i<=grid_in[0].nx; i++) {
      ind0 = j*(grid_in[0].nx+2)+i;
      ind1 = j*(grid_in[0].nx+2)+1;
      ind2 = 1*(grid_in[0].nx+2)+i;
      if(fabs( grid_in[0].lont[ind0]-grid_in[0].lont[ind2] ) > EPSLN10 ||
	 fabs( grid_in[0].latt[ind0]-grid_in[0].latt[ind1] ) > EPSLN10 )
	mpp_error("fregrid: extrapolate is limited to lat-lon grid");
	  
    }
  }

  /* when vertical_interp is set, extrapolate must be set */
  if( vertical_interp) extrapolate = 1;
  /* vertical_interp and extrapolate is not supported for vector interpolation */
  if( nvector > 0) {
    if(vertical_interp) mpp_error("fregrid: vertical_interp is not supported for vector fields");
    if(extrapolate) mpp_error("fregrid: extrapolate is not supported for vector fields");
  }
  
  if(remap_file) set_remap_file(ntiles_out, mosaic_out, remap_file, interp, &opcode, save_weight_only);  
  
  if(!save_weight_only) {
    file_in   = (File_config *)malloc(ntiles_in *sizeof(File_config));
    file_out  = (File_config *)malloc(ntiles_out*sizeof(File_config));
 
    if(nfiles == 2) {
      file2_in   = (File_config *)malloc(ntiles_in *sizeof(File_config));
      file2_out  = (File_config *)malloc(ntiles_out*sizeof(File_config));
    }
    if(nscalar > 0) {
      scalar_in  = (Field_config *)malloc(ntiles_in *sizeof(Field_config));
      scalar_out = (Field_config *)malloc(ntiles_out *sizeof(Field_config));
    }
    if(nvector > 0) {
      u_in  = (Field_config *)malloc(ntiles_in *sizeof(Field_config));
      u_out = (Field_config *)malloc(ntiles_out *sizeof(Field_config));    
      v_in  = (Field_config *)malloc(ntiles_in *sizeof(Field_config));
      v_out = (Field_config *)malloc(ntiles_out *sizeof(Field_config));
    }
  
    set_mosaic_data_file(ntiles_in, mosaic_in, dir_in, file_in,  input_file[0]);
    set_mosaic_data_file(ntiles_out, mosaic_out, dir_out, file_out, output_file[0]);

    vgrid_out.nz = 0;
    vgrid_in.nz = 0;
    if(vertical_interp) {
      get_output_vgrid(&vgrid_out, dst_vgrid);
      get_input_vgrid(&vgrid_in, file_in[0].name, scalar_name[0]);
      setup_vertical_interp(&vgrid_in, &vgrid_out);
    }
  
    if(nfiles == 2) {
      set_mosaic_data_file(ntiles_in, mosaic_in, dir_in, file2_in,  input_file[1]);
      set_mosaic_data_file(ntiles_out, mosaic_out, dir_out, file2_out, output_file[1]);    
    }

    for(n=0; n<ntiles_in; n++) file_in[n].fid = mpp_open(file_in[n].name, MPP_READ);

    nscalar_orig = nscalar;
    /* filter field with interp_method = "none"  */
    nscalar = 0;
    for(n=0; n<nscalar_orig; n++) {
      int vid;
      char remap_method[STRING];

      vid = mpp_get_varid(file_in[0].fid, scalar_name[n]);
      if(mpp_var_att_exist(file_in[0].fid, vid, "interp_method")) {
        mpp_get_var_att(file_in[0].fid, vid, "interp_method", remap_method);
        if(strcmp(remap_method, "none") && strcmp(remap_method, "NONE") && strcmp(remap_method, "None")) {
          strcpy(scalar_name_remap[nscalar], scalar_name[n]);
          nscalar++;
        }
      }
      else {
        strcpy(scalar_name_remap[nscalar], scalar_name[n]);
          nscalar++;
      }
    }
    set_field_struct ( ntiles_in,   scalar_in,   nscalar, scalar_name_remap[0], file_in);
    set_field_struct ( ntiles_out,  scalar_out,  nscalar, scalar_name_remap[0], file_out);
    set_field_struct ( ntiles_in,   u_in,        nvector, u_name[0], file_in);
    set_field_struct ( ntiles_out,  u_out,       nvector, u_name[0], file_out);
    if(nfiles == 1) {
      set_field_struct ( ntiles_in,   v_in,        nvector, v_name[0], file_in);
      set_field_struct ( ntiles_out,  v_out,       nvector, v_name[0], file_out);
    }
    else {
      set_field_struct ( ntiles_in,   v_in,        nvector, v_name[0], file2_in);
      set_field_struct ( ntiles_out,  v_out,       nvector, v_name[0], file2_out);
    }

    get_input_metadata(ntiles_in, nfiles, file_in, file2_in, scalar_in, u_in, v_in, grid_in,
		       kbegin, kend, lbegin, lend, opcode, associated_file_dir);

    set_weight_inf( ntiles_in, grid_in, weight_file, weight_field, file_in->has_cell_measure_att);
    
    set_output_metadata(ntiles_in, nfiles, file_in, file2_in, scalar_in, u_in, v_in,
			ntiles_out, file_out, file2_out, scalar_out, u_out, v_out, grid_out, &vgrid_out, history, tagname, opcode);

    if(debug) print_mem_usage("After set_output_metadata");
    /* when the interp_method specified through command line is CONSERVE_ORDER1, but the interp_method in the source file
       field attribute is CONSERVE_ORDER2, need to modify the interp_method value */
    if(opcode & CONSERVE_ORDER1) {
      for(l=0; l<nscalar; l++) {
	if(scalar_out->var[l].interp_method == CONSERVE_ORDER2) {
	  if(mpp_pe() == mpp_root_pe())printf("NOTE from fregrid: even though the interp_method specified through command line is "
					      "conserve_order1, the interp_method is reset to conserve_order2 because some fields in "
					      "the source data have interp_method attribute value conserve_order2");
	  opcode = opcode & ~CONSERVE_ORDER1;
	  opcode |= CONSERVE_ORDER2;
	  break;
	}
      }
    }
    if(opcode & CONSERVE_ORDER1) {
      for(l=0; l<nvector; l++) {
	if(u_out->var[l].interp_method == CONSERVE_ORDER2) {
	  if(mpp_pe() == mpp_root_pe())printf("NOTE from fregrid: even though the interp_method specified through command line is "
					      "conserve_order1, the interp_method is reset to conserve_order2 because some fields in "
					      "the source data have interp_method attribute value conserve_order2");
	  opcode = opcode & ~CONSERVE_ORDER1;
	  opcode |= CONSERVE_ORDER2;
	  break;
	}
      }
    }    
  }

  /* preparing for the interpolation, if remapping information exist, read it from remap_file,
     otherwise create the remapping information and write it to remap_file
  */

  if(debug) time_start = clock();
   if( opcode & BILINEAR ) /* bilinear interpolation from cubic to lalon */
     setup_bilinear_interp(ntiles_in, grid_in, ntiles_out, grid_out, interp, opcode );
   else
     setup_conserve_interp(ntiles_in, grid_in, ntiles_out, grid_out, interp, opcode);
   if(debug) {
     time_end = clock();
     time_setup_interp = 1.0*(time_end - time_start)/CLOCKS_PER_SEC;
     print_mem_usage("After setup interp");
   }
   if(debug) {
      print_time("get_input_grid", time_get_in_grid);
      print_time("get_output_grid", time_get_out_grid);
      print_time("setup_interp", time_setup_interp);
   }
   if(save_weight_only) {
     if(mpp_pe() == mpp_root_pe() ) {
       printf("NOTE: Successfully running fregrid and the following files which store weight information are generated.\n");
       for(n=0; n<ntiles_out; n++) {
	 printf("****%s\n", interp[n].remap_file);
       }
     }
     mpp_end();
     return 0;     
   }
  
   if(nscalar > 0) {
     get_field_attribute(ntiles_in, scalar_in);
     copy_field_attribute(ntiles_out, scalar_in, scalar_out);
   }
   
   if(nvector > 0) {
     get_field_attribute(ntiles_in, u_in);
     get_field_attribute(ntiles_in, v_in);
     copy_field_attribute(ntiles_out, u_in, u_out);
     copy_field_attribute(ntiles_out, v_in, v_out);
   }


  
   /* set time step to 1, only test scalar field now, nz need to be 1 */
   if(test_case) {
     if(nscalar != 1 || nvector != 0) mpp_error("fregrid: when test_case is specified, nscalar must be 1 and nvector must be 0");
     if(scalar_in->var->nz != 1) mpp_error("fregrid: when test_case is specified, number of vertical level must be 1");
     file_in->nt = 1;
     file_out->nt = 1;
   }
   
  /* Then doing the regridding */
  for(m=0; m<file_in->nt; m++) {
    int memsize, level_z, level_n, level_t;

    write_output_time(ntiles_out, file_out, m);
    if(nfiles > 1) write_output_time(ntiles_out, file2_out, m);
    
    /* first interp scalar variable */
    for(l=0; l<nscalar; l++) {
      if( !scalar_in->var[l].has_taxis && m>0) continue;
      level_t = m + scalar_in->var[l].lstart;
      /*--- to reduce memory usage, we are only do remapping for on horizontal level one time */
      for(level_n =0; level_n < scalar_in->var[l].nn; level_n++) {
	if(extrapolate) {
	  get_input_data(ntiles_in, scalar_in, grid_in, bound_T, l, -1, level_n, level_t, extrapolate, stop_crit);
	  allocate_field_data(ntiles_out, scalar_out, grid_out, scalar_in->var[l].nz);
	  if( opcode & BILINEAR ) 
	    do_scalar_bilinear_interp(interp, l, ntiles_in, grid_in, grid_out, scalar_in, scalar_out, finer_step, fill_missing);
	  else
	    do_scalar_conserve_interp(interp, l, ntiles_in, grid_in, ntiles_out, grid_out, scalar_in, scalar_out, opcode, scalar_in->var[l].nz);
          if(vertical_interp) do_vertical_interp(&vgrid_in, &vgrid_out, grid_out, scalar_out, l);
	  write_field_data(ntiles_out, scalar_out, grid_out, l, -1, level_n, m);
	  if(scalar_out->var[l].interp_method == CONSERVE_ORDER2) {
	    for(n=0; n<ntiles_in; n++) {
	      free(scalar_in[n].grad_x);
	      free(scalar_in[n].grad_y);
	    }
	  }
	  for(n=0; n<ntiles_in; n++) free(scalar_in[n].data);
	  for(n=0; n<ntiles_out; n++) free(scalar_out[n].data);
	}
	else {
	  for(level_z=scalar_in->var[l].kstart; level_z <= scalar_in->var[l].kend; level_z++)
	    {	    
              if(debug) time_start = clock();
              if(test_case)
		get_test_input_data(test_case, test_param, ntiles_in, scalar_in, grid_in, bound_T, opcode);
	      else
		get_input_data(ntiles_in, scalar_in, grid_in, bound_T, l, level_z, level_n, level_t, extrapolate, stop_crit);
              if(debug) {
	        time_end = clock();
		time_get_input += 1.0*(time_end - time_start)/CLOCKS_PER_SEC;
	      }

	      allocate_field_data(ntiles_out, scalar_out, grid_out, 1);
	      if(debug) time_start = clock();
	      if( opcode & BILINEAR ) 
		do_scalar_bilinear_interp(interp, l, ntiles_in, grid_in, grid_out, scalar_in, scalar_out, finer_step, fill_missing);
	      else
		do_scalar_conserve_interp(interp, l, ntiles_in, grid_in, ntiles_out, grid_out, scalar_in, scalar_out, opcode,1);
              if(debug) {
		time_end = clock();
		time_do_interp += 1.0*(time_end - time_start)/CLOCKS_PER_SEC;
	      }

	      if(debug) time_start = clock();
	      write_field_data(ntiles_out, scalar_out, grid_out, l, level_z, level_n, m);
	      if(debug) {
		time_end = clock();
	        time_write += 1.0*(time_end - time_start)/CLOCKS_PER_SEC;
	      }
	      if(scalar_out->var[l].interp_method == CONSERVE_ORDER2) {
		for(n=0; n<ntiles_in; n++) {
		  free(scalar_in[n].grad_x);
		  free(scalar_in[n].grad_y);
		  free(scalar_in[n].grad_mask);
		}
	      }
	      for(n=0; n<ntiles_in; n++) free(scalar_in[n].data);
	      for(n=0; n<ntiles_out; n++) free(scalar_out[n].data);
	    }
	}
      }
    }
   if(debug) print_mem_usage("After do interp");
    /* then interp vector field */
    for(l=0; l<nvector; l++) {
      if( !u_in[n].var[l].has_taxis && m>0) continue;
      level_t = m + u_in->var[l].lstart;
      get_input_data(ntiles_in, u_in, grid_in, bound_T, l, level_z, level_n, level_t, extrapolate, stop_crit);
      get_input_data(ntiles_in, v_in, grid_in, bound_T, l, level_z, level_n, level_t, extrapolate, stop_crit);
      allocate_field_data(ntiles_out, u_out, grid_out, u_in[n].var[l].nz);
      allocate_field_data(ntiles_out, v_out, grid_out, u_in[n].var[l].nz);
      if( opcode & BILINEAR )
	do_vector_bilinear_interp(interp, l, ntiles_in, grid_in, ntiles_out, grid_out, u_in, v_in, u_out, v_out, finer_step, fill_missing);
      else
	do_vector_conserve_interp(interp, l, ntiles_in, grid_in, ntiles_out, grid_out, u_in, v_in, u_out, v_out, opcode);
      
      write_field_data(ntiles_out, u_out, grid_out, l, level_z, level_n, m);
      write_field_data(ntiles_out, v_out, grid_out, l, level_z, level_n, m);
      for(n=0; n<ntiles_in; n++) {
	free(u_in[n].data);
	free(v_in[n].data);
      }
      for(n=0; n<ntiles_out; n++) {
	free(u_out[n].data);
	free(v_out[n].data);
      }      
    }
  }

  if(debug) {
    print_time("get_input", time_get_input);
    print_time("do_interp", time_do_interp);
    print_time("write_data", time_write);
  }
  
  if(mpp_pe() == mpp_root_pe() ) {
    printf("Successfully running fregrid and the following output file are generated.\n");
    for(n=0; n<ntiles_out; n++) {
      mpp_close(file_out[n].fid);
      printf("****%s\n", file_out[n].name);
      if( nfiles > 1 ) {
	mpp_close(file2_out[n].fid);
	printf("****%s\n", file2_out[n].name);
      }
    }
  }
      
  mpp_end();
  return 0;
  
} /* end of main */
  

  
  
