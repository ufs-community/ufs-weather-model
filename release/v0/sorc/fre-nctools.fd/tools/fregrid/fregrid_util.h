#ifndef FREGRID_UTIL_H_
#define FREGRID_UTIL_H_
#include "globals.h"
void set_mosaic_data_file(int ntiles, const char *mosaic_file, const char *dir, File_config *file,
			  const char *filename);
void set_field_struct(int ntiles, Field_config *field, int nvar, char * varname, File_config *file);
void get_input_grid(int ntiles, Grid_config *grid, Bound_config *bound, const char *mosaic_file, unsigned int opcode,
                    int *great_circl_algorithm, int save_weight_only);
void get_output_grid_from_mosaic(int ntiles, Grid_config *grid, const char *mosaic_file, unsigned int opcode,
                                 int *great_circl_algorithm);
void get_output_grid_by_size(int ntiles, Grid_config *grid, double lonbegin, double lonend, double latbegin, double latend, 
                             int nlon, int nlat, int finer_steps, int center_y, unsigned int opcode);
void get_input_metadata(int ntiles, int nfiles, File_config *file1, File_config *file2,
		        Field_config *scalar, Field_config *u_comp, Field_config *v_comp,
			const Grid_config *grid, int kbegin, int kend, int lbegin, int lend, unsigned int opcode,
                        char *associated_file_dir);
void set_output_metadata (int ntiles_in, int nfiles, const File_config *file1_in, const File_config *file2_in,
			  const Field_config *scalar_in, const Field_config *u_in, const Field_config *v_in,
			  int ntiles_out, File_config *file1_out, File_config *file2_out, Field_config *scalar_out,
			  Field_config *u_out, Field_config *v_out, const Grid_config *grid_out, const VGrid_config *vgrid_out, 
                          const char *history, const char *tagname, unsigned int opcode);
void get_field_attribute( int ntiles, Field_config *field);
void copy_field_attribute( int ntiles_out, Field_config *field_in, Field_config *field_out);
void set_remap_file( int ntiles, const char *mosaic_file, const char *remap_file, Interp_config *interp, unsigned int *opcode, int save_weight_only);
void write_output_time(int ntiles, File_config *output, int level);
void get_input_data(int ntiles, Field_config *field, Grid_config *grid, Bound_config *bound,
		    int varid, int level_z, int level_n, int level_t, int extrapolate, double stop_crit);
void get_test_input_data(char *test_case, double test_param, int ntiles, Field_config *field,
			 Grid_config *grid, Bound_config *bound, unsigned int opcode);
void get_input_output_cell_area(int ntiles_in, Grid_config *grid_in, int ntiles_out, Grid_config *grid_out, unsigned int opcode);
void allocate_field_data(int ntiles, Field_config *field, Grid_config *grid, int nz);
void write_field_data(int ntiles, Field_config *field, Grid_config *grid, int varid, int level_z, int level_n, int level_t);
void set_weight_inf(int ntiles, Grid_config *grid, const char *weight_file, const char *weight_field, int has_cell_measure_att);
void get_output_vgrid( VGrid_config *vgrid, const char *vgrid_file );
void get_input_vgrid( VGrid_config *vgrid, const char *vgrid_file, const char *field );
void setup_vertical_interp(VGrid_config *vgrid_in, VGrid_config *vgrid_out);
void do_vertical_interp(VGrid_config *vgrid_in, VGrid_config *vgrid_out, Grid_config *grid_out, Field_config *field, int varid);
#endif
