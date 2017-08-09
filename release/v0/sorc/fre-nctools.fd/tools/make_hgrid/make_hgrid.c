/*
  Copyright 2011 NOAA Geophysical Fluid Dynamics Lab, Princeton, NJ
  This program is distributed under the terms of the GNU General Public
  License. See the file COPYING contained in this directory

  This program generates various types of horizontal grids in netCDF file format

 AUTHOR: Zhi Liang (Zhi.Liang@noaa.gov)
          NOAA Geophysical Fluid Dynamics Lab, Princeton, NJ
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <netcdf.h>
#include "create_hgrid.h"
#include "mpp.h"
#include "mpp_domain.h"
#include "mpp_io.h"
#include "tool_util.h"
#include "mosaic_util.h"

const int MAXBOUNDS = 100;
const int STRINGLEN = 255;
#define REGULAR_LONLAT_GRID    1
#define TRIPOLAR_GRID          2
#define FROM_FILE              3
#define SIMPLE_CARTESIAN_GRID  4
#define SPECTRAL_GRID          5
#define CONFORMAL_CUBIC_GRID   6
#define GNOMONIC_ED            7
#define F_PLANE_GRID           8
#define BETA_PLANE_GRID        9
#define MISSING_VALUE           (-9999.)
int my_grid_type = 0;

char *usage[] = {
  "",
  "                                                                                 ",
  "                                                                                 ",
  "                    Usage of make_hgrid                                          ",
  "                                                                                 ",
  "   make_hgrid --grid_type grid_type --my_grid_file my_grid_file                  ",
  "                  --nxbnds nxbnds --nybnds nybnds                                ",
  "                  --xbnds x(1),...,x(nxbnds) --ybnds y(1),...,y(nybnds)          ",
  "                  --nlon nlon(1),...nlon(nxbnds-1)                               ",
  "                  --nlat nlat(1),...nlat(nybnds-1)                               ",
  "                  --dlon dlon(1),...dlon(nxbnds)                                 ",
  "                  --dlat dlat(1),...dlat(nybnds)                                 ",
  "                  --lat_join lat_join --num_lon num_lon --nratio nratio          ",
  "                  --simple_dx simple_dx --simple_dy simple_dy                    ",
  "                  --grid_name gridname --center center --verbose --shift_fac #   ",
  "                  --do_schmidt --stretch_fac # --target_lon # --target_lat #     ",
  "                  --nest_grid --parent_tile # --refine_ratio # --halo #          ",
  "                  --istart_nest # --iend_nest # --jstart_nest # --jend_nest #    ",
  "                  --great_circle_algorithm --out_halo #                          ",
  "                                                                                 ",
  "   This program can generate different types of horizontal grid. The             ",
  "   output data is on supergrid ( model grid size x refinement(=2) ).  For        ",
  "   'cubic_grid', six grid files which contain the grid information for           ",
  "   each tile will be generate, otherwise one file will be generated you          ",
  "   can specify the grid type through --grid_type. The value of grid_type         ",
  "   can be 'from_file', 'spectral_grid', 'spherical_grid',                        ",
  "   'conformal_cubic_grid', 'gnomonic_ed' or 'simple_cartesian_grid', with        ",
  "   default value 'spherical_grid'. --nlon and --nlat must be specified to        ",
  "   indicate supergrid size ( for cubic_grid, --nlat is not needed because        ",
  "   nlat has the same value as nlon.  Besides --nlon, --nlat, other               ",
  "   optional and requirement arguments for each type are,                         ",
  "                                                                                 ",
  "   Two algorithms are provided for creating 'tripolar_grid' and                  ",
  "   'regular_lonlat_grid'. Monotonic cubic interpolation and legacy               ",
  "   algorithm using cosine function to configure grid cell location.  The         ",
  "   monotonic cubic interpolation is developed by Russell Fiedler and the         ",
  "   detail of the alrogithm is explained in                                       ",
  "   http://www.mathworks.com/moler/interp.pdf. More details of the                ",
  "   algorithm is available in                                                     ",
  "                                                                                 ",
  "      F. N. Fritsch and R. E. Carlson, Monotone Piecewise Cubic                  ",
  "      Interpolation, SIAM Journal on Numerical Analysis, 17 (1980),              ",
  "      pp. 238-246.                                                               ",
  "                                                                                 ",
  "      D. Kahaner, C. Moler, and S. Nash, Numerical Methods and Software,         ",
  "      Prentice Hall, Englewood Cliffs, NJ, 1989.                                 ",
  "                                                                                 ",
  "   --nlon --nlat need to be specified to use the monotonic cubic                 ",
  "   algorithm.  The legacy algorithm was developed by Ron Pacanowski. This        ",
  "   algorithm uses a cosine function to do the interpolation.                     ",
  "                                                                                 ",
  "   --dlon and --dlat need to be specified to use the legacy                      ",
  "   algorithm. Monotonic cubic interpolation is a higher order                    ",
  "   interpolation and will produce smoother grid distance. It is strongly         ",
  "   suggested to use the monotonic cubic interpolation by specifying              ",
  "   argument --nlon and --nlat.                                                   ",
  "                                                                                 ",
  "   1. 'from_file':              --my_grid_file must be specified. The grid       ",
  "                                specified in my_grid_file should be super grid   ",
  "                                vertex.                                          ",
  "   2. 'spectral_grid':          no other optional or required arguments.         ",
  "   3. 'regular_lonlat_grid':    --nxbnds, --nybnds --xbnds, --ybnds, must be     ",
  "                                specified to define the grid bounds.             ",
  "   4. 'tripolar_grid':          --nxbnds, --nybnds, --xbnds, --ybnds, must be    ",
  "                                specified to define the grid bounds. --lat_join  ",
  "                                is optional with default value 65.               ",
  "   5  'conformal_cubic_grid':   --nratio is optional argument.                   ",
  "   6  'gnomonic_ed'          :  equal distance gnomonic cubic grid.              ",
  "   7. 'simple_cartesian_grid':  --xbnds, --ybnds must be specified to define     ",
  "                                the grid bounds location and grid size. number   ",
  "                                of bounds must be 2 in both and x and            ",
  "                                y-direction. --simple_dx and --simple_dy must be ",
  "                                specified to specify uniform cell length.        ",
  "   8  'f_plane_grid':           For setting geometric fractors according         ",
  "                                to f-plane. f_plane_latitude need to be specified",
  "   9  'beta_plane_grid':        For setting geometric fractors according         ",
  "                                to  beta plane. f_plane_latitude need to be      ",
  "                                specified                                        ",
  "                                                                                 ",
  "   make_hgrid take the following flags                                           ",
  "                                                                                 ",
  "   --grid_type grid_type      specify type of topography. See above for          ",
  "                              grid type option.                                  "
  "                                                                                 ",
  "   --my_grid_file file        when this flag is present, the program will read   ",
  "                              grid information from 'my_grid_file'. The file     ",
  "                              format can be ascii file or netcdf file. Multiple  ",
  "                              file entry are allowed but the number should be    ",
  "                              less than MAXBOUNDS.                                ",
  "                                                                                 ",
  "   --nxbnds nxbnds            Specify number of zonal regions for varying        ",
  "                              resolution.                                        ",
  "                                                                                 ",
  "   --nybnds nybnds            Specify number of meridinal regions for varying    ",
  "                              resolution.                                        ",
  "                                                                                 ",
  "   --xbnds x(1),.,x(nxbnds)   Specify boundaries for defining zonal regions of   ",
  "                              varying resolution. When --tripolar is present,    ",
  "                              x also defines the longitude of the two new poles. ",
  "                              nxbnds must be 2 and lon_start = x(1),             ",
  "                              lon_end = x(nxbnds) are longitude of the two       ",
  "                              new poles.                                         ",
  "                                                                                 ",
  "   --ybnds y(1),.,y(nybnds)   Specify boundaries for defining meridional         ",
  "                              regions of varying resolution                      ",
  "                                                                                 ",
  "   --nlon nlon(1),..,nlon(nxbnds-1) Number of model grid points(supergrid) for   ",
  "                                    each zonal regions of varying resolution.    ",
  "                                                                                 ", 
  "   --nlat nlat(1),..,nlat(nybnds-1) Number of model grid points(supergid) for    ",
  "                                    each meridinal regions of varying resolution.",
  "                                                                                 ",
  "   --dlon dlon(1),..,dlon(nxbnds)   nominal resolution of zonal regions          ",
  "                                                                                 ",
  "   --dlat dlat(1),..,dlat(nybnds)   nominal resolution of meridional regions     ",
  "                                                                                 ",
  "   --lat_join lat_join        Specify latitude for joining spherical and rotated ",
  "                              bipolar grid. Default value is 65 degree.          ",
  "                                                                                 ",
  "   --nratio nratio            Speicify the refinement ratio when calculating     ",
  "                              cell length and area of supergrid.                 ",
  "                                                                                 ",
  "   --simple_dx dimple_dx      Specify the uniform cell length in x-direction for ",
  "                              simple cartesian grid.                             ",
  "                                                                                 ",
  "   --simple_dy dimple_dy      Specify the uniform cell length in y-direction for ",
  "                              simple cartesian grid.                             ",  
  "                                                                                 ",
  "   --grid_name grid_name      Specify the grid name. The output grid file name   ",
  "                              will be grid_name.nc if there is one tile and      ",
  "                              grid_name.tile#.nc if there is more than one tile. ",
  "                              The default value will be horizontal_grid.         ",
  "                                                                                 ",
  "   --center center            Specify the center location of grid. The valid     ",
  "                              entry will be 'none', 't_cell' or 'c_cell' with    ",
  "                              default value 'none'. The grid refinement is       ",
  "                              assumed to be 2 in x and y-direction when center   ",
  "                              is not 'none'. 'c_cell' should be used for the grid",
  "                              used in MOM.                                       ",
  "                                                                                 ",
  "   --shift_fac                shift west by 180/shift_fac. Default value is 18.  ",
  "                                                                                 ",
  "   --do_schmidt               Set to do Schmidt transformation to create         ",
  "                              stretched grid. When do_schmidt is set, the        ",
  "                              following must be set: --stretch_factor,           ",
  "                              --target_lon and --target_lat.                     ",
  "                                                                                 ",
  "   --stretch_factor #         Stretching factor for the grid                     ",
  "                                                                                 ",
  "   --target_lon #             center longitude of the highest resolution tile    ",
  "                                                                                 ",
  "   --target_lat #             center latitude of the highest resolution tile     ",
  "                                                                                 ",
  "   --nest_grid                set to create nest grid as well as the global grid.",
  "                              This option could only be set when grid_type is    ",
  "                              'gnomonic_ed'. When it is set, besides 6 tile grid ",
  "                              files created, there is one more nest grid with    ",
  "                              file name = $grid_name.tile${parent_tile}.nest.nc  ",
  "                                                                                 ",
  "   --parent_tile #            Specify the parent tile number of nest grid.       ",
  "                                                                                 ",
  "   --refine_ratio #           Specify the refinement ratio for nest grid.        ",
  "                                                                                 ",
  "   --istart_nest #            Specify the starting i-direction index of nest     ",
  "                              grid in parent tile supergrid(Fortran index).      ",
  "                                                                                 ",
  "   --iend_nest #              Specify the ending i-direction index of nest       ",
  "                              grid in parent tile supergrid(Fortran index).      ",  
  "                                                                                 ",
  "   --jstart_nest #            Specify the starting j-direction index of nest     ",
  "                              grid in parent tile supergrid(Fortran index).      ",
  "                                                                                 ",
  "   --jend_nest #              Specify the ending j-direction index of nest       ",
  "                              grid in parent tile supergrid(Fortran index).      ",
  "                                                                                 ",
  "   --halo #                   halo size to used in the atmosphere cubic sphere   ",
  "                              model. It only needs to be specified when          ",
  "                              --nest_grid is set.                                ",
  "                                                                                 ",
  "   --great_circle_algorithm   When specified, great_circle_algorithm will be     ",
  "                              used to compute grid cell area.                    ",
  "                                                                                 ",
  "   --out_halo #               extra halo size data to be written out. This is    ",
  "                              only works for gnomonic_ed.                        ",
  "                                                                                 ",
  "   --verbose                  Will print out running time message when this      ",
  "                              option is set. Otherwise the run will be silent    ",
  "                              when there is no error.                            ",
  "                                                                                 ",
  "   Example                                                                       ",      
  "                                                                                 ",
  "                                                                                 ",
  "   1. generating regular lon-lat grid (supergrid size 60x20)                     ",
  "      > make_hgrid --grid_type regular_lonlat_grid --nxbnd 2 --nybnd 2           ",
  "        --xbnd 0,30 --ybnd 50,60  --nlon 60 --nlat 20                            ",
  "                                                                                 ",
  "   2. generating tripolar grid with various grid resolution and C-cell centered  ",
  "      using monotonic bi-cub spline interpolation.                               ",
  "      > make_hgrid --grid_type tripolar_grid --nxbnd 2 --nybnd 7 --xbnd -280,80  ",
  "                   --ybnd -82,-30,-10,0,10,30,90 --nlon 720                      ",
  "                   --nlat 104,48,40,40,48,120 --grid_name om3_grid               ",
  "                   --center c_cell                                               ",
  "                                                                                 ",
  "   3. generating tripolar grid with various grid resolution and C-cell centered  ",
  "      using legacy algorithm (create GFDL CM2/ocean-like grid)                   ",
  "      > make_hgrid --grid_type tripolar_grid --nxbnd 2 --nybnd 7 --xbnd -280,80  ",
  "                   --ybnd -82,-30,-10,0,10,30,90 --dlon 1.0,1.0                  ",
  "                   --dlat 1.0,1.0,0.6666667,0.3333333,0.6666667,1.0,1.0          ",
  "                   --grid_name om3_grid --center c_cell                          ",
  "                                                                                 ",  
  "   4. generating simple cartesian grid(supergrid size 20x20)                     ",
  "      > make_hgrid --grid_type simple_cartesian_grid --xbnd 0,30 --ybnd 50,60    ",
  "                   --nlon 20 --nlat 20  --simple_dx 1000 --simple_dy 1000        ",
  "                                                                                 ",
  "   5. generating conformal cubic grid. (supergrid size 60x60 for each tile)      ",
  "      > make_hgrid --grid_type conformal_cubic_grid --nlon 60 --nratio 2         ",
  "                                                                                 ",
  "   6. generating gnomonic cubic grid with equal_dist_face_edge(C48 grid)         ",
  "      > make_hgrid --grid_type gnomonic_ed --nlon 96                             ",
  "                                                                                 ",
  "   7. generating gnomonic cubic streched grid.                                   ",
  "      > make_hgrid --grid_type gnomonic_ed --nlon 180 --do_schmidt               ",
  "                   --stretch_factor 3 --target_lat 40. --target_lon 20.          ",
  "                                                                                 ",
  "   8. generating spectral grid. (supergrid size 128x64)                          ",
  "      > make_hgrid --grid_type spectral_grid --nlon 128 --nlat 64                ",
  "                                                                                 ",
  "   9. Through    user-defined grids                                              ",
  "      > make_hgrid --grid_type from_file --my_grid_file my_grid_file             ",
  "                   --nlon 4 --nlat 4                                             ",
  "                                                                                 ",
  "       contents of sample my_grid_file                                           ",
  "         The first line of my_grid_file will be text ( will be ignored)          ",     
  "         followed by nlon+1 lines of real value of x-direction supergrid bound   ",
  "         location. Then another line of text ( will be ignored), followed by     ",
  "         nlat+1 lines of real value of y-direction supergrid bound location.     ",
  "                                                                                 ",     
  "         For example:                                                            ",
  "                                                                                 ",
  "            x-grid                                                               ",
  "            0.0                                                                  ",
  "            5.0                                                                  ",
  "            10.0                                                                 ",
  "            15.0                                                                 ",
  "            20.0                                                                 ",
  "            y-grid                                                               ",
  "            -10                                                                  ",
  "            10                                                                   ",
  "            20                                                                   ",
  "            30                                                                   ",
  "            40                                                                   ", 
  "                                                                                 ",
  "   10. generating f_plane_grids                                                   ",
  "      > make_hgrid --grid_type f_plane_grid --f_plane_latitude 55 --nxbnd 2      ",
  "                   --nybnd 2 --xbnd 0,30 --ybnd 50,60  --nlon 60 --nlat 20       ",
  "                                                                                 ",
  " A note on generating cyclic regular lon-lat grids when center = 'c_cell':-      ",
  " It is possible to have an improperly centered boundary unless care is taken to  ",
  " ensure local symmetry at the join.                                              ",
  " A correctly formed grid is only guaranteed if the first three values of the     ",
  " --xbnd argument mirror the last 3 and the first two 'nlon' arguments mirror the ",
  " last 2.                                                                         ",
  "                                                                                 ",
  " For example for a global grid make_hgrid should be invoked as                   ",
  "       > make_hgrid --grid_type regular_lonlat_grid ...                          ",
  "                    --xbnd 0,X1,X2,...,360-X2,360-X1,360                         ",
  "                    --nlon N1,N2,...,N2,N1 --center c_cell                       ",
  "                                                                                 ",
  " As an example                                                                   ",
  "                                                                                 ",
  "       > make_hgrid --grid_type regular_lonlat_grid --nxbnd 7 --nybnd 2          ",
  "                    --xbnd 0,15,30,300,330,345,360 --ybnd 50,60                  ",
  "                    --nlon 4,2,6,4,2,4 --nlat 2 --center c_cell                  ",
  "                                                                                 ",
  "                                                                                 ",
  " results in a valid cyclic grid whereas (note the second last value of nlon)     ",
  "                                                                                 ",
  "       > make_hgrid --grid_type regular_lonlat_grid --nxbnd 7 --nybnd 2          ",
  "                    --xbnd 0,15,30,300,330,345,360 --ybnd 50,60                  ",
  "                    --nlon 4,2,6,4,4,4 --nlat 2 --center c_cell                  ",
  "                                                                                 ",
  "                                                                                 ",
  " is not properly centered across 0,360                                           ",
  "                                                                                 ",
  " An informational message is issued if the leftmost and rightmost  resolutions   ",
  " differ  by more than 1 part in 10E6                                             ",
  "                                                                                 ",
  "",
  NULL };

char grid_version[] = "0.2";
char tagname[] = "$Name: fre-nctools-bronx-10 $";


void fill_cubic_grid_halo(int nx, int ny, int halo, double *data, double *data1_all,
			  double *data2_all, int tile, int ioff, int joff)
{
  int lw, le, ls, ln;
  int ntiles,nxp,nyp,nxph,nyph,i,j;

  
  nxp = nx+ioff;
  nyp = ny+joff;
  nxph = nx+ioff+2*halo;
  nyph = ny+joff+2*halo;
  
  for(i=0; i<nxph*nyph; i++) data[i] = MISSING_VALUE;
  
  /* first copy computing domain data */
  for(j=1; j<=nyp; j++)
    for(i=1; i<=nxp; i++)
      data[j*nxph+i] = data1_all[tile*nxp*nyp+(j-1)*nxp+(i-1)];
  
  ntiles=6;
  if(tile%2 == 1) { /* tile 2, 4, 6 */
    lw = (tile+ntiles-1)%ntiles;
    le = (tile+ntiles+2)%ntiles;
    ls = (tile+ntiles-2)%ntiles;
    ln = (tile+ntiles+1)%ntiles;
    for(j=1; j<=nyp; j++) {
      data[j*nxph] = data1_all[lw*nxp*nyp+(j-1)*nxp+nx-1]; /* west halo */
      data[j*nxph+nxp+1] = data2_all[le*nxp*nyp+ioff*nxp+nyp-j]; /*east halo */
    }

    for(i=1; i<=nxp; i++) {
      data[i] = data2_all[ls*nxp*nyp+(nxp-i)*nxp+(nx-1)]; /*south */
      data[(nyp+1)*nxph+i] = data1_all[ln*nxp*nyp+joff*nxp+i-1]; /*north */
    }
  }
  else { /* tile 1, 3, 5 */
    lw = (tile+ntiles-2)%ntiles;
    le = (tile+ntiles+1)%ntiles;
    ls = (tile+ntiles-1)%ntiles;
    ln = (tile+ntiles+2)%ntiles;
    for(j=1; j<=nyp; j++) {
      data[j*nxph] = data2_all[lw*nxp*nyp+(ny-1)*nxp+nyp-j]; /* west halo */
      data[j*nxph+nxp+1] = data1_all[le*nxp*nyp+(j-1)*nxp+joff]; /*east halo */
    }

    for(i=1; i<=nxp; i++) {
      data[i] = data1_all[ls*nxp*nyp+(ny-1)*nxp+i-1]; /*south */
      data[(nyp+1)*nxph+i] = data2_all[ln*nxp*nyp+(nxp-i)*nxp+joff]; /*north */
    }    

  }
}

int main(int argc, char* argv[])
{
  int  nratio = 1;
  char method[32] = "conformal";
  char orientation[32] = "center_pole";
  int  nxbnds=2, nybnds=2, nxbnds0=0, nybnds0=0, nxbnds1=0, nybnds1=0, nxbnds2=0, nybnds2=0;
  int  nxbnds3=0, nybnds3=0;
  double xbnds[MAXBOUNDS], ybnds[MAXBOUNDS];
  double dx_bnds[MAXBOUNDS], dy_bnds[MAXBOUNDS];
  int nlon[MAXBOUNDS-1], nlat[MAXBOUNDS-1];
  char grid_type[128]="regular_lonlat_grid";
  char my_grid_file[MAXBOUNDS][STRINGLEN];
  double f_plane_latitude  = 100;
  double lat_join=65.;
  double shift_fac = 18.0;
  int do_schmidt = 0;
  double stretch_factor = 0.0;
  double target_lon   = 0.0;
  double target_lat   = 0.0;
  int    nest_grid    = 0;
  int    parent_tile  = 0;
  int    refine_ratio = 0;
  int    istart_nest  = 0;
  int    iend_nest    = 0;
  int    jstart_nest  = 0;
  int    jend_nest    = 0;
  int    halo = 0;
  int    out_halo=0;
  int    present_stretch_factor = 0;
  int    present_target_lon = 0;
  int    present_target_lat = 0;
  int    use_great_circle_algorithm = 0;
  unsigned int verbose = 0;
  double simple_dx=0, simple_dy=0;
  int nx, ny, nxp, nyp, ntiles=1, ntiles_global=1;
  int *nxl=NULL, *nyl=NULL;
  unsigned int ntiles_file;
  double *x=NULL, *y=NULL, *dx=NULL, *dy=NULL, *angle_dx=NULL, *angle_dy=NULL, *area=NULL;

  int isc, iec, jsc, jec;
  int  use_legacy;
  char history[2560];
  char gridname[128] = "horizontal_grid";
  char center[32] = "none";
  char geometry[32] = "spherical";
  char projection[32] = "none";
  char arcx[32] = "small_circle";
  char north_pole_tile[32] = "0.0 90.0";
  char north_pole_arcx[32] = "0.0 90.0";
  char discretization[32]  = "logically_rectangular";
  char conformal[32]       = "true";
  char mesg[256], str[128];
  char entry[MAXBOUNDS*STRINGLEN];
  int n, errflg, c, i;  
  int option_index;

  static struct option long_options[] = {
    {"grid_type",       required_argument, NULL, 'a'},
    {"my_grid_file",    required_argument, NULL, 'b'},
    {"nxbnds",          required_argument, NULL, 'c'},
    {"nybnds",          required_argument, NULL, 'd'},
    {"xbnds",           required_argument, NULL, 'e'},
    {"ybnds",           required_argument, NULL, 'f'},
    {"nlon",            required_argument, NULL, 'g'},
    {"nlat",            required_argument, NULL, 'i'},
    {"lat_join",        required_argument, NULL, 'j'},
    {"nratio",          required_argument, NULL, 'k'},
    {"simple_dx",       required_argument, NULL, 'l'},
    {"simple_dy",       required_argument, NULL, 'm'},
    {"grid_name",       required_argument, NULL, 'q'},
    {"center",          required_argument, NULL, 'r'},
    {"dlon",            required_argument, NULL, 's'},
    {"dlat",            required_argument, NULL, 't'},
    {"f_plane_latitude",required_argument, NULL, 'u'},
    {"do_schmidt",      no_argument,       NULL, 'w'},
    {"stretch_factor",  required_argument, NULL, 'x'},
    {"target_lon",      required_argument, NULL, 'y'},
    {"target_lat",      required_argument, NULL, 'z'},
    {"nest_grid",       no_argument,       NULL, 'A'},
    {"refine_ratio",    required_argument, NULL, 'B'},
    {"parent_tile",     required_argument, NULL, 'C'},
    {"istart_nest",     required_argument, NULL, 'D'},
    {"iend_nest",       required_argument, NULL, 'E'},
    {"jstart_nest",     required_argument, NULL, 'F'},
    {"jend_nest",       required_argument, NULL, 'G'},
    {"halo",            required_argument, NULL, 'H'},
    {"shift_fac",       required_argument, NULL, 'I'},
    {"great_circle_algorithm", no_argument, NULL, 'J'},
    {"out_halo",        required_argument, NULL, 'K'},
    {"help",            no_argument,       NULL, 'h'},
    {"verbose",         no_argument,       NULL, 'v'},
    {0, 0, 0, 0},
  };

  /* start parallel */
  mpp_init(&argc, &argv);
  mpp_domain_init();  

  /* There is no need to run this tool in parallel, so we limit this tool
     to be run on single processor*/
  if(mpp_npes() > 1) mpp_error( "make_hgrid: make_hgrid must be run one processor, contact developer");
  
  /*
   * process command line
   */
  errflg = argc <3;
  
  while ((c = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
    switch (c) {
    case 'a':
      strcpy(grid_type, optarg);
      break;
    case 'b':
      strcpy(entry, optarg);
      tokenize(entry, ",", STRINGLEN, MAXBOUNDS, (char *)my_grid_file, &ntiles_file);
      break;
    case 'c':
      nxbnds0 = atoi(optarg);
      break;
    case 'd':
      nybnds0 = atoi(optarg);
      break;        
    case 'e':
      strcpy(entry, optarg);
      nxbnds1 = get_double_entry(entry, xbnds);
      break;
    case 'f':
      strcpy(entry, optarg);
      nybnds1 = get_double_entry(entry, ybnds);
      break;
    case 'g':
      strcpy(entry, optarg);
      nxbnds2 = get_int_entry(entry, nlon);
      break;
    case 'i':
      strcpy(entry, optarg);
      nybnds2 = get_int_entry(entry, nlat);
      break;
    case 'j':
      lat_join = atof(optarg);
      break; 
    case 'k':
      nratio = atoi(optarg);
      break;      
    case 'l':
      simple_dx = atof(optarg);
      break;
    case 'm':
      simple_dy = atof(optarg);
      break;
    case 'q':
      strcpy(gridname, optarg);
      break;
    case 'r':
      strcpy(center, optarg);
      break;
    case 's':
      strcpy(entry, optarg);
      nxbnds3 = get_double_entry(entry, dx_bnds);
      break;
    case 't':
      strcpy(entry, optarg);
      nybnds3 = get_double_entry(entry, dy_bnds);
      break;
    case 'u':
      f_plane_latitude = atof(optarg);
      break;
    case 'w':
      do_schmidt = 1;
      break;
    case 'x':
      present_stretch_factor = 1;
      stretch_factor = atof(optarg);
      break;
    case 'y':
      present_target_lon = 1;
      target_lon = atof(optarg);
      break;
    case 'z':
      present_target_lat = 1;
      target_lat = atof(optarg);
      break;
    case 'A':
      nest_grid = 1;
      break;
    case 'B':
      refine_ratio = atoi(optarg);
      break;
    case 'C':
      parent_tile = atoi(optarg);
      break;
    case 'D':
      istart_nest = atoi(optarg);
      break;
    case 'E':
      iend_nest = atoi(optarg);
      break;
    case 'F':
      jstart_nest = atoi(optarg);
      break;
    case 'G':
      jend_nest = atoi(optarg);
      break;
    case 'H':
      halo = atoi(optarg);
      break;
    case 'I':
      shift_fac = atof(optarg);
      break;
    case 'J':
      use_great_circle_algorithm = 1;
      break;
    case 'K':
      out_halo = atoi(optarg);
      break;      
    case 'v':
      verbose = 1;
      break;
    case 'h':
      errflg++;
      break;
    case '?':
      errflg++;      
    }      
  }
  
  if (errflg ) {
    char **u = usage;
    while (*u) { fprintf(stderr, "%s\n", *u); u++; }
    exit(2);
  }  

  /* define history to be the history in the grid file */
  strcpy(history,argv[0]);

  for(i=1;i<argc;i++) {
    strcat(history, " ");
    strcat(history, argv[i]);
  }
  
  if(mpp_pe() == mpp_root_pe() && verbose) printf("==>NOTE: the grid type is %s\n",grid_type);

  if(strcmp(grid_type,"regular_lonlat_grid") ==0 )
    my_grid_type = REGULAR_LONLAT_GRID;
  else if(strcmp(grid_type,"tripolar_grid") ==0 )
    my_grid_type = TRIPOLAR_GRID;
  else if(strcmp(grid_type,"from_file")==0 )
    my_grid_type = FROM_FILE;
  else if(strcmp(grid_type, "simple_cartesian_grid")==0)
    my_grid_type = SIMPLE_CARTESIAN_GRID;
  else if(strcmp(grid_type, "spectral_grid") ==0 )
    my_grid_type = SPECTRAL_GRID;
  else if(strcmp(grid_type, "conformal_cubic_grid") ==0 )
    my_grid_type = CONFORMAL_CUBIC_GRID;
  else if(strcmp(grid_type, "gnomonic_ed") ==0 )
    my_grid_type = GNOMONIC_ED;
  else if(strcmp(grid_type, "f_plane_grid") == 0 )
    my_grid_type = F_PLANE_GRID;
  else if(strcmp(grid_type, "beta_plane_grid") == 0 )
    my_grid_type = BETA_PLANE_GRID;
  else
    mpp_error("make_hgrid: only grid_type = 'regular_lonlat_grid', 'tripolar_grid', 'from_file', "
	      "'gnomonic_ed', 'conformal_cubic_grid', 'simple_cartesian_grid', "
              "'spectral_grid', 'f_plane_grid' and 'beta_plane_grid' is implemented");  

  if(my_grid_type != GNOMONIC_ED && out_halo  != 0)
    mpp_error("make_hgrid: out_halo should not be set when grid_type = gnomonic_ed");
  if(out_halo !=0 && out_halo != 1)
    mpp_error("make_hgrid: out_halo should be 0 or 1");
  
 if( my_grid_type != GNOMONIC_ED && do_schmidt )
    mpp_error("make_hgrid: --do_schmidt should not be set when grid_type is not 'gnomonic_ed'");
  
  use_legacy = 0;
  /* check the command-line arguments to make sure the value are suitable */
  if( my_grid_type == REGULAR_LONLAT_GRID || my_grid_type == TRIPOLAR_GRID ||
      my_grid_type == F_PLANE_GRID || my_grid_type == BETA_PLANE_GRID ) {
    int num_specify;
    nxbnds = nxbnds0; nybnds = nybnds0;
    if( nxbnds <2 || nybnds < 2) mpp_error("make_hgrid: grid type is 'regular_lonlat_grid', 'tripolar_grid', 'f_plane_grid' or 'beta_plane_grid', "
					   "both nxbnds and nybnds should be no less than 2");
    if( nxbnds != nxbnds1 ) mpp_error("make_hgrid: grid type is 'regular_lonlat_grid, 'tripolar_grid', 'f_plane_grid' or 'beta_plane_grid', "
				      "nxbnds does not match number of entry in xbnds");
    if( nybnds != nybnds1 ) mpp_error("make_hgrid: grid type is 'regular_lonlat_grid, 'tripolar_grid', 'f_plane_grid' or 'beta_plane_grid', "
				      "nybnds does not match number of entry in ybnds");
    num_specify = 0;
    if( nxbnds2 > 0 && nybnds2 > 0 ) num_specify ++;
    if( nxbnds3 > 0 && nybnds3 > 0 ) {
      num_specify ++;
      use_legacy = 1;
    }

    if( num_specify == 0 ) mpp_error("make_hgrid: grid type is 'regular_lonlat_grid', 'tripolar_grid', 'f_plane_grid' or 'beta_plane_grid', "
				     "need to specify one of the pair --nlon --nlat or --dlon --dlat");
    if( num_specify == 2 ) mpp_error("make_hgrid: grid type is 'regular_lonlat_grid', 'tripolar_grid', 'f_plane_grid' or 'beta_plane_grid', "
				     "can not specify both --nlon --nlat and --dlon --dlat");
    if( use_legacy ) {
      if( nxbnds != nxbnds3 ) mpp_error("make_hgrid: grid type is 'tripolar_grid', 'tripolar_grid', 'f_plane_grid' or 'beta_plane_grid', "
					"nxbnds does not match number of entry in dlon");
      if( nybnds != nybnds3 ) mpp_error("make_hgrid: grid type is 'tripolar_grid', 'tripolar_grid', 'f_plane_grid' or 'beta_plane_grid', "
					"nybnds does not match number of entry in dlat");
    }
    else {
      if( nxbnds != nxbnds2+1 ) mpp_error("make_hgrid: grid type is 'tripolar_grid', 'tripolar_grid', 'f_plane_grid' or 'beta_plane_grid', "
					"nxbnds does not match number of entry in nlon");
      if( nybnds != nybnds2+1 ) mpp_error("make_hgrid: grid type is 'tripolar_grid', 'tripolar_grid', 'f_plane_grid' or 'beta_plane_grid', "
					"nybnds does not match number of entry in nlat");
    }
  }  


  if( my_grid_type == CONFORMAL_CUBIC_GRID || my_grid_type == GNOMONIC_ED ) {
    ntiles = 6;
    ntiles_global = 6;
  }
  
  if(  my_grid_type != GNOMONIC_ED && nest_grid )
    mpp_error("make_hgrid: --nest_grid can be set only when grid_type = 'gnomonic_ed'");
  
  if( my_grid_type == TRIPOLAR_GRID ) {
    strcpy(projection, "tripolar");
    if( nxbnds != 2) mpp_error("make_hgrid: grid type is 'tripolar_grid', nxbnds should be 2");
  }
  else if( my_grid_type == FROM_FILE ) {
    /* For ascii file, nlon and nlat should be specified through --nlon, --nlat
       For netcdf file, grid resolution will be read from grid file
    */
    
    if(ntiles_file == 0) mpp_error("make_hgrid: grid_type is 'from_file', but my_grid_file is not specified");
    ntiles = ntiles_file;
    for(n=0; n<ntiles; n++) {
      if(strstr(my_grid_file[n],".nc") ) {
	/* get the grid size for each tile, the grid is on model grid, should need to multiply by 2 */
	int fid;
	fid = mpp_open(my_grid_file[n], MPP_READ);
        if(mpp_dim_exist(fid, "grid_xt") ) {
	  if( mpp_dim_exist(fid, "grid_yt") == 0)
	    mpp_error("make_hgrid: grid_yt should be a dimension when grid_xt is a dimension");
	  nlon[n] = mpp_get_dimlen(fid, "grid_xt")*2;
	  nlat[n] = mpp_get_dimlen(fid, "grid_yt")*2;
	}
	else if(mpp_dim_exist(fid, "rlon") ) {
	  if( mpp_dim_exist(fid, "rlat") == 0)
	    mpp_error("make_hgrid: rlat should be a dimension when rlon is a dimension");
          nlon[n] = mpp_get_dimlen(fid, "rlon")*2;
	  nlat[n] = mpp_get_dimlen(fid, "rlat")*2;
	}
	else if(mpp_dim_exist(fid, "lon") ) {
	  if( mpp_dim_exist(fid, "lat") == 0)
	    mpp_error("make_hgrid: lat should be a dimension when lon is a dimension");
          nlon[n] = mpp_get_dimlen(fid, "lon")*2;
	  nlat[n] = mpp_get_dimlen(fid, "lat")*2;
	}
	else if(mpp_dim_exist(fid, "i") ) {
	  if( mpp_dim_exist(fid, "j") == 0)
	    mpp_error("make_hgrid: j should be a dimension when i is a dimension");
          nlon[n] = mpp_get_dimlen(fid, "i")*2;
	  nlat[n] = mpp_get_dimlen(fid, "j")*2;
	}
	else if(mpp_dim_exist(fid, "x") ) {
	  if( mpp_dim_exist(fid, "y") == 0)
	    mpp_error("make_hgrid: y should be a dimension when x is a dimension");
          nlon[n] = mpp_get_dimlen(fid, "x")*2;
	  nlat[n] = mpp_get_dimlen(fid, "y")*2;
	}
	
	else {
	  mpp_error("make_hgrid: none of grid_xt, rlon, lon, x, and i is a dimension in input file");
	}	
	mpp_close(fid);
      }
      else {
	if(nxbnds2 != ntiles || nybnds2 != ntiles ) mpp_error("make_hgrid: grid type is 'from_file', number entry entered "
						"through --nlon and --nlat should be equal to number of files "
							  "specified through --my_grid_file");
      }
    }
    	/* for simplify purpose, currently we assume all the tile have the same grid size */
    for(n=1; n<ntiles; n++) {
      if( nlon[n] != nlon[0] || nlat[n] != nlat[0])  mpp_error("make_hgrid: grid_type is from_file, all the tiles should "
							       "have same grid size, contact developer");
    }
  }
  else if( my_grid_type == SIMPLE_CARTESIAN_GRID ) {
    strcpy(geometry, "planar");
    strcpy(north_pole_tile, "none");
    if(nxbnds1 != 2 || nybnds1 != 2 ) mpp_error("make_hgrid: grid type is 'simple_cartesian_grid', number entry entered "
						"through --xbnds and --ybnds should be 2");
    if(nxbnds2 != 1 || nybnds2 != 1 ) mpp_error("make_hgrid: grid type is 'simple_cartesian_grid', number entry entered "
						"through --nlon and --nlat should be 1");
    if(simple_dx == 0 || simple_dy == 0) mpp_error("make_hgrid: grid_type is 'simple_cartesian_grid', "
						   "both simple_dx and simple_dy both should be specified");
  }
  else if( my_grid_type == SPECTRAL_GRID ) {
    if(nxbnds2 != 1 || nybnds2 != 1 ) mpp_error("make_hgrid: grid type is 'spectral_grid', number entry entered "
						"through --nlon and --nlat should be 1");    
  }
  else if( my_grid_type == CONFORMAL_CUBIC_GRID ){
    strcpy(projection, "cube_gnomonic");
    strcpy(conformal, "FALSE");
    if(nxbnds2 != 1 ) mpp_error("make_hgrid: grid type is 'conformal_cubic_grid', number entry entered "
				"through --nlon should be 1");
    if(nratio < 1) mpp_error("make_hgrid: grid type is 'conformal_cubic_grid', nratio should be a positive integer");
  }
  else if( my_grid_type == GNOMONIC_ED ) {
    strcpy(projection, "cube_gnomonic");
    strcpy(conformal, "FALSE");
    if( do_schmidt ) {
      if( present_stretch_factor == 0 || present_target_lon == 0 || present_target_lat == 0 )
	mpp_error("make_hgrid: grid type is 'gnomonic_ed, --stretch_factor, --target_lon "
		  "and --target_lat must be set when --do_schmidt is set");
    }

    if(nest_grid) {
      if(refine_ratio == 0) mpp_error("make_hgrid: --refine_ratio must be set when --nest_grid is set");
      if(parent_tile == 0 && mpp_pe()==mpp_root_pe()) {
	 printf("NOTE from make_hgrid: parent_tile is 0, the output grid will have resolution refine_ration*nlon");
      }
      else {
	if(istart_nest == 0) mpp_error("make_hgrid: --istart_nest must be set when --nest_grid is set");
	if(iend_nest == 0) mpp_error("make_hgrid: --iend_nest must be set when --nest_grid is set");
	if(jstart_nest == 0) mpp_error("make_hgrid: --jstart_nest must be set when --nest_grid is set");
	if(jend_nest == 0) mpp_error("make_hgrid: --jend_nest must be set when --nest_grid is set");
	if(halo == 0 ) mpp_error("make_hgrid: --halo must be set when --nest_grid is set");
	ntiles++;   /* one more tile for the nest region */
      }
    }
    if(nxbnds2 != 1 ) mpp_error("make_hgrid: grid type is 'gnomonic_cubic_grid', number entry entered "
				"through --nlon should be 1");
  }
  else if( my_grid_type == F_PLANE_GRID ||  my_grid_type == BETA_PLANE_GRID) {
    if(f_plane_latitude > 90 || f_plane_latitude < -90.)
      mpp_error("make_hgrid: f_plane_latitude should be between -90 and 90.");
    if(f_plane_latitude > ybnds[nybnds-1] || f_plane_latitude < ybnds[0] ) {
      if(mpp_pe() == mpp_root_pe())
	printf("Warning from make_hgrid: f_plane_latitude is not inside the latitude range of the grid\n");
    }
    if(mpp_pe() == mpp_root_pe())
      printf("make_hgrid: setting geometric factor according to f-plane with f_plane_latitude = %g\n", f_plane_latitude );
  }

  nxl = (int *)malloc(ntiles*sizeof(int));
  nyl = (int *)malloc(ntiles*sizeof(int));
  
  /* get super grid size */
  if(use_legacy) {
    nxl[0] = get_legacy_grid_size(nxbnds, xbnds, dx_bnds);
    nyl[0] = get_legacy_grid_size(nybnds, ybnds, dy_bnds);
  }
  else {
    if( my_grid_type == GNOMONIC_ED || my_grid_type == CONFORMAL_CUBIC_GRID ) {
      for(n=0; n<ntiles_global; n++) {	
	nxl[n] = nlon[0];
	nyl[n] = nxl[n];
	if(nest_grid && parent_tile == 0) {
	  nxl[n] *= refine_ratio;
	  nyl[n] *= refine_ratio;
	}
      }
      if(ntiles > ntiles_global) {
        nxl[ntiles_global] = (iend_nest-istart_nest+1)*refine_ratio;
        nyl[ntiles_global] = (jend_nest-jstart_nest+1)*refine_ratio;
      }
    }
    else {
      nxl[0] = 0;
      nyl[0] = 0;
      for(n=0; n<nxbnds-1; n++) nxl[0] += nlon[n];
      for(n=0; n<nybnds-1; n++) nyl[0] += nlat[n];  
    }
  }
  nx = nxl[0];
  ny = nyl[0];
  nxp = nx + 1;
  nyp = ny + 1;

  if(strcmp(center,"none") && strcmp(center,"c_cell") && strcmp(center,"t_cell") )
    mpp_error("make_hgrid: center should be 'none', 'c_cell' or 't_cell' ");
  
 /* create grid information */
  {
    int size1, size2, size3, size4;
    
    size1 = nxp*nyp*ntiles;
    size2 = nxp*(nyp+1)*ntiles;
    size3 = (nxp+1)*nyp*ntiles;
    size4 = nxp    *nyp*ntiles;
  
    if( nest_grid) { /* nest grid is the last tile */
      size1 += nxl[ntiles-1]*nyl[ntiles-1];
      size2 += nxl[ntiles-1]*nyl[ntiles-1];
      size3 += nxl[ntiles-1]*nyl[ntiles-1];
      size4 += nxl[ntiles-1]*nyl[ntiles-1];
    }

    x        = (double *) malloc(size1*sizeof(double));
    y        = (double *) malloc(size1*sizeof(double));
    dx       = (double *) malloc(size2*sizeof(double));
    dy       = (double *) malloc(size3*sizeof(double));
    area     = (double *) malloc(size4*sizeof(double));
    angle_dx = (double *) malloc(size1*sizeof(double));
    if( strcmp(conformal,"true") !=0 )
      angle_dy = (double *) malloc(size1*sizeof(double));
  }

  isc = 0;
  iec = nx-1;
  jsc = 0;
  jec = ny-1;

  if(my_grid_type==REGULAR_LONLAT_GRID)
    create_regular_lonlat_grid(&nxbnds, &nybnds, xbnds, ybnds, nlon, nlat, dx_bnds, dy_bnds,
			       use_legacy, &isc, &iec, &jsc, &jec, x, y, dx, dy, area,
			       angle_dx, center, use_great_circle_algorithm);
  else if(my_grid_type==TRIPOLAR_GRID)
    create_tripolar_grid(&nxbnds, &nybnds, xbnds, ybnds, nlon, nlat, dx_bnds, dy_bnds,
			 use_legacy, &lat_join, &isc, &iec, &jsc, &jec, x, y, dx, dy,
			 area, angle_dx, center, verbose, use_great_circle_algorithm);
  else if(my_grid_type==FROM_FILE) {
    for(n=0; n<ntiles; n++) {
      int n1, n2, n3, n4;
      n1 = n * nxp * nyp;
      n2 = n * nx  * nyp;
      n3 = n * nxp * ny;
      n4 = n * nx  * ny;
      create_grid_from_file(my_grid_file[n], &nx, &ny, x+n1, y+n1, dx+n2, dy+n3, area+n4, angle_dx+n1, use_great_circle_algorithm);
    }
  }
  else if(my_grid_type==SIMPLE_CARTESIAN_GRID)
    create_simple_cartesian_grid(xbnds, ybnds, &nx, &ny, &simple_dx, &simple_dy, &isc, &iec, &jsc, &jec,
				 x, y, dx, dy, area, angle_dx );
  else if(my_grid_type==SPECTRAL_GRID)
    create_spectral_grid(&nx, &ny, &isc, &iec, &jsc, &jec, x, y, dx, dy, area, angle_dx, use_great_circle_algorithm );
  else if(my_grid_type==CONFORMAL_CUBIC_GRID)
    create_conformal_cubic_grid(&nx, &nratio, method, orientation, x, y, dx, dy, area, angle_dx, angle_dy );
  else if(my_grid_type==GNOMONIC_ED)
    create_gnomonic_cubic_grid(grid_type, nxl, nyl, x, y, dx, dy, area, angle_dx, angle_dy,
			       shift_fac, do_schmidt, stretch_factor, target_lon, target_lat,
			       nest_grid, parent_tile, refine_ratio,
			       istart_nest, iend_nest, jstart_nest, jend_nest, halo );
  else if((my_grid_type==F_PLANE_GRID) || (my_grid_type==BETA_PLANE_GRID))
    create_f_plane_grid(&nxbnds, &nybnds, xbnds, ybnds, nlon, nlat, dx_bnds, dy_bnds,
			use_legacy, f_plane_latitude, &isc, &iec, &jsc, &jec, x, y, dx, dy, area, angle_dx, center);
 
  /* write out data */
  {
    int fid, id_tile, id_x, id_y, id_dx, id_dy, id_area, id_angle_dx, id_angle_dy, id_arcx;
    int dimlist[5], dims[2], i, j, l, m;
    size_t start[4], nwrite[4];
    char tilename[128] = "";
    char outfile[128] = "";
    int pos_c, pos_e, pos_n, pos_t;

    pos_c = 0;
    pos_e = 0;
    pos_t = 0;
    pos_n = 0;    
    for(n=0 ; n< ntiles; n++) {

      sprintf(tilename, "tile%d", n+1);
      if(ntiles>1)
	sprintf(outfile, "%s.tile%d.nc", gridname, n+1);
      else
	sprintf(outfile, "%s.nc", gridname);
      fid = mpp_open(outfile, MPP_WRITE);
      /* define dimenison */
      nx = nxl[n];
      ny = nyl[n];
      nxp = nx+1;
      nyp = ny+1;
      dimlist[0] = mpp_def_dim(fid, "string", STRINGLEN);
      dimlist[1] = mpp_def_dim(fid, "nx", nx+2*out_halo);
      dimlist[2] = mpp_def_dim(fid, "ny", ny+2*out_halo);
      dimlist[3] = mpp_def_dim(fid, "nxp", nxp+2*out_halo);
      dimlist[4] = mpp_def_dim(fid, "nyp", nyp+2*out_halo);
      /* define variable */
      if( strcmp(north_pole_tile, "none") == 0) /* no north pole, then no projection */
	id_tile = mpp_def_var(fid, "tile", MPP_CHAR, 1, dimlist, 4, "standard_name", "grid_tile_spec",
			      "geometry", geometry, "discretization", discretization, "conformal", conformal );
      else if( strcmp(projection, "none") == 0) 
	id_tile = mpp_def_var(fid, "tile", MPP_CHAR, 1, dimlist, 5, "standard_name", "grid_tile_spec",
			      "geometry", geometry, "north_pole", north_pole_tile, "discretization",
			      discretization, "conformal", conformal );
      else
	id_tile = mpp_def_var(fid, "tile", MPP_CHAR, 1, dimlist, 6, "standard_name", "grid_tile_spec",
			      "geometry", geometry, "north_pole", north_pole_tile, "projection", projection,
			      "discretization", discretization, "conformal", conformal );
      
      dims[0] = dimlist[4]; dims[1] = dimlist[3];
      id_x = mpp_def_var(fid, "x", MPP_DOUBLE, 2, dims, 2, "standard_name", "geographic_longitude",
        		 "units", "degree_east");
      if(out_halo>0) mpp_def_var_att_double(fid, id_x, "_FillValue", MISSING_VALUE);
      id_y = mpp_def_var(fid, "y", MPP_DOUBLE, 2, dims, 2, "standard_name", "geographic_latitude",
			 "units", "degree_north");
      if(out_halo>0) mpp_def_var_att_double(fid, id_y, "_FillValue", MISSING_VALUE);
      dims[0] = dimlist[4]; dims[1] = dimlist[1];
      id_dx = mpp_def_var(fid, "dx", MPP_DOUBLE, 2, dims, 2, "standard_name", "grid_edge_x_distance",
			  "units", "meters");
      if(out_halo>0) mpp_def_var_att_double(fid, id_dx, "_FillValue", MISSING_VALUE);
      dims[0] = dimlist[2]; dims[1] = dimlist[3];
      id_dy = mpp_def_var(fid, "dy", MPP_DOUBLE, 2, dims, 2, "standard_name", "grid_edge_y_distance",
			  "units", "meters");
      if(out_halo>0) mpp_def_var_att_double(fid, id_dy, "_FillValue", MISSING_VALUE);
      dims[0] = dimlist[2]; dims[1] = dimlist[1];
      id_area = mpp_def_var(fid, "area", MPP_DOUBLE, 2, dims, 2, "standard_name", "grid_cell_area",
			    "units", "m2" );
      if(out_halo>0) mpp_def_var_att_double(fid, id_area, "_FillValue", MISSING_VALUE);
      dims[0] = dimlist[4]; dims[1] = dimlist[3];
      id_angle_dx = mpp_def_var(fid, "angle_dx", MPP_DOUBLE, 2, dims, 2, "standard_name",
				    "grid_vertex_x_angle_WRT_geographic_east", "units", "degrees_east");
      if(out_halo>0) mpp_def_var_att_double(fid, id_angle_dx, "_FillValue", MISSING_VALUE);
      if(strcmp(conformal, "true") != 0) {
	id_angle_dy = mpp_def_var(fid, "angle_dy", MPP_DOUBLE, 2, dims, 2, "standard_name",
				  "grid_vertex_y_angle_WRT_geographic_north", "units", "degrees_north");
	if(out_halo>0) mpp_def_var_att_double(fid, id_angle_dy, "_FillValue", MISSING_VALUE);
      }
      if( strcmp(north_pole_arcx, "none") == 0)
	id_arcx = mpp_def_var(fid, "arcx", MPP_CHAR, 1, dimlist, 1, "standard_name", "grid_edge_x_arc_type" );
      else
	id_arcx = mpp_def_var(fid, "arcx", MPP_CHAR, 1, dimlist, 2, "standard_name", "grid_edge_x_arc_type",
			      "north_pole", north_pole_arcx );
      mpp_def_global_att(fid, "grid_version", grid_version);
      mpp_def_global_att(fid, "code_version", tagname);
      if(use_great_circle_algorithm) mpp_def_global_att(fid, "great_circle_algorithm", "TRUE");
      if(n>=ntiles_global) mpp_def_global_att(fid, "nest_grid", "TRUE");
      mpp_def_global_att(fid, "history", history);
      
      mpp_end_def(fid);
      for(m=0; m<4; m++) { start[m] = 0; nwrite[m] = 0; }
      nwrite[0] = strlen(tilename);
      mpp_put_var_value_block(fid, id_tile, start, nwrite, tilename );

      if(out_halo ==0) {
        mpp_put_var_value(fid, id_x, x+pos_c);
	mpp_put_var_value(fid, id_y, y+pos_c);
	mpp_put_var_value(fid, id_dx, dx+pos_n);
	mpp_put_var_value(fid, id_dy, dy+pos_e);
	mpp_put_var_value(fid, id_area, area+pos_t);
	mpp_put_var_value(fid, id_angle_dx, angle_dx+pos_c);
	if(strcmp(conformal, "true") != 0) mpp_put_var_value(fid, id_angle_dy, angle_dy+pos_c);
      }
      else {
	double *tmp;

	tmp = (double *)malloc((nxp+2*out_halo)*(nyp+2*out_halo)*sizeof(double));
	fill_cubic_grid_halo(nx,ny,out_halo,tmp,x,x,n,1,1);
	mpp_put_var_value(fid, id_x, tmp);
	fill_cubic_grid_halo(nx,ny,out_halo,tmp,y,y,n,1,1);
	mpp_put_var_value(fid, id_y, tmp);
        fill_cubic_grid_halo(nx,ny,out_halo,tmp,angle_dx,angle_dx,n,1,1);
	mpp_put_var_value(fid, id_angle_dx, tmp);
	if(strcmp(conformal, "true") != 0) {
	  fill_cubic_grid_halo(nx,ny,out_halo,tmp,angle_dy,angle_dy,n,1,1);
	mpp_put_var_value(fid, id_angle_dy, tmp);
	}
	
        fill_cubic_grid_halo(nx,ny,out_halo,tmp,dx,dy,n,0,1);
	mpp_put_var_value(fid, id_dx, tmp);
        fill_cubic_grid_halo(nx,ny,out_halo,tmp,dy,dx,n,1,0);
	mpp_put_var_value(fid, id_dy, tmp);
        fill_cubic_grid_halo(nx,ny,out_halo,tmp,area,area,n,0,1);
	mpp_put_var_value(fid, id_area, tmp);
	free(tmp);
      }
	
      nwrite[0] = strlen(arcx);
      mpp_put_var_value_block(fid, id_arcx, start, nwrite, arcx );
      mpp_close(fid);
      pos_c += nxp*nyp;
      pos_e += nxp*ny;
      pos_n += nx*nyp;
      pos_t += nx*ny;
    }
  }

  free(x);
  free(y);
  free(dx);
  free(dy);
  free(area);
  free(angle_dx);
  if(strcmp(conformal, "true") != 0) free(angle_dy);
  if(mpp_pe() == mpp_root_pe() && verbose) printf("generate_grid is run successfully. \n");

  mpp_end();

  return 0;
  
};  /* end of main */


