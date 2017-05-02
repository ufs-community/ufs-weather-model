/***********************************************************************
                      mosaic_util.h
    This header file provide some utilities routine that will be used in many tools.
    
    contact: Zhi.Liang@noaa.gov
***********************************************************************/
#ifndef MOSAIC_UTIL_H_
#define MOSAIC_UTIL_H_

#ifndef RANGE_CHECK_CRITERIA
#define RANGE_CHECK_CRITERIA 0.05
#endif

#define min(a,b) (a<b ? a:b)
#define max(a,b) (a>b ? a:b)
#define SMALL_VALUE ( 1.e-10 )
struct Node{
  double x, y, z, u, u_clip;
  int intersect; /* indicate if this point is an intersection, 0 = no, 1= yes, 2=both intersect and vertices */ 
  int inbound;      /* -1 uninitialized, 0 coincident, 1 outbound, 2 inbound */
  int initialized; /* = 0 means empty list */
  int isInside;   /* = 1 means one point is inside the other polygon, 0 is not, -1 undecided. */
  int subj_index; /* the index of subject point that an intersection follow. */
  int clip_index; /* the index of clip point that an intersection follow */
  struct Node *Next;
};


void error_handler(const char *msg);
int nearest_index(double value, const double *array, int ia);
int lon_fix(double *x, double *y, int n_in, double tlon);
double minval_double(int size, const double *data);
double maxval_double(int size, const double *data);
double avgval_double(int size, const double *data);
void latlon2xyz(int size, const double *lon, const double *lat, double *x, double *y, double *z); 
void xyz2latlon(int size, const double *x, const double *y, const double *z, double *lon, double *lat);
double box_area(double ll_lon, double ll_lat, double ur_lon, double ur_lat);
double poly_area(const double lon[], const double lat[], int n);
double poly_area_dimensionless(const double lon[], const double lat[], int n);
double poly_area_no_adjust(const double x[], const double y[], int n);
int fix_lon(double lon[], double lat[], int n, double tlon);
void tokenize(const char * const string, const char *tokens, unsigned int varlen,
	      unsigned int maxvar, char * pstring, unsigned int * const nstr);
double great_circle_distance(double *p1, double *p2);
double spherical_excess_area(const double* p_ll, const double* p_ul,
			     const double* p_lr, const double* p_ur, double radius);
void vect_cross(const double *p1, const double *p2, double *e );
double spherical_angle(const double *v1, const double *v2, const double *v3);
void normalize_vect(double *e);
void unit_vect_latlon(int size, const double *lon, const double *lat, double *vlon, double *vlat);
double great_circle_area(int n, const double *x, const double *y, const double *z);
double * cross(const double *p1, const double *p2);
double dot(const double *p1, const double *p2);
int intersect_tri_with_line(const double *plane, const double *l1, const double *l2, double *p,
			     double *t);
int invert_matrix_3x3(long double m[], long double m_inv[]);
void mult(long double m[], long double v[], long double out_v[]);
double metric(const double *p);
int insidePolygon(struct Node *node, struct Node *list );
int inside_a_polygon( double *lon1, double *lat1, int *npts, double *lon2, double *lat2);

void rewindList(void);
struct Node *getNext();
void initNode(struct Node *node);
void addEnd(struct Node *list, double x, double y, double z, int intersect, double u, int inbound, int inside);
int addIntersect(struct Node *list, double x, double y, double z, int intersect, double u1, double u2, 
                int inbound, int is1, int ie1, int is2, int ie2);
int length(struct Node *list);
int samePoint(double x1, double y1, double z1, double x2, double y2, double z2);
int sameNode(struct Node node1, struct Node node2);
void addNode(struct Node *list, struct Node nodeIn);
struct Node *getNode(struct Node *list, struct Node inNode);
struct Node *getNextNode(struct Node *list);
void copyNode(struct Node *node_out, struct Node node_in);
void printNode(struct Node *list, char *str);
int intersectInList(struct Node *list, double x, double y, double z);
void insertAfter(struct Node *list, double x, double y, double z, int intersect, double u, int inbound,
		 double x2, double y2, double z2);
double gridArea(struct Node *grid);
int isIntersect(struct Node node);
int getInbound( struct Node node );
struct Node *getLast(struct Node *list);
int getFirstInbound( struct Node *list, struct Node *nodeOut);
void getCoordinate(struct Node node, double *x, double *y, double *z);
void getCoordinates(struct Node *node, double *p);
void setCoordinate(struct Node *node, double x, double y, double z);
void setInbound(struct Node *interList, struct Node *list);
int isInside(struct Node *node);
void set_reproduce_siena_true(void);


#endif
