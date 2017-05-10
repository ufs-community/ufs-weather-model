#ifndef READ_MOSAIC_H_
#define READ_MOSAIC_H_

int read_mosaic_xgrid_size( const char *xgrid_file );
#ifdef OVERLOAD_R4
void read_mosaic_xgrid_order1(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, float *area );
void read_mosaic_xgrid_order1_region(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, float *area, int *isc, int *iec );
void read_mosaic_xgrid_order2(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, 
                              float *area, float *di, float *dj );
float get_global_area(void);
#else
void read_mosaic_xgrid_order1(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, double *area );
void read_mosaic_xgrid_order1_region(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, double *area, int *isc, int *iec );
void read_mosaic_xgrid_order2(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, 
                              double *area, double *di, double *dj );
double get_global_area(void);
#endif
int read_mosaic_ntiles(const char *mosaic_file);
int read_mosaic_ncontacts(const char *mosaic_file);
void read_mosaic_grid_sizes(const char *mosaic_file, int *nx, int *ny);
void read_mosaic_contact(const char *mosaic_file, int *tile1, int *tile2, int *istart1, int *iend1,
			 int *jstart1, int *jend1, int *istart2, int *iend2, int *jstart2, int *jend2);
void read_mosaic_grid_data(const char *mosaic_file, const char *name, int nx, int ny,
                           double *data, int level, int ioff, int joff);

#endif
