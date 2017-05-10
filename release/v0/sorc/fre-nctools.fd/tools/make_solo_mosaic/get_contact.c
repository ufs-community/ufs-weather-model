#include <stdlib.h>
#include <math.h>
#include "mosaic_util.h"
#include "get_contact.h"


/************************************************************************************************
         int get_align_contact(
  This routine will return number of algined contacts bewteen two tiles (line contact).
  This routine assume the starting and ending points of the  contact line are coincidence with
  the grid points of both tiles. lrg_rectangle tiles are assumed. 
  It will return the contact information, which includes 

*************************************************************************************************/

double* east_bound(const double *data, int nx, int ny);
double* west_bound(const double *data, int nx, int ny);
double* south_bound(const double *data, int nx, int ny);
double* north_bound(const double *data, int nx, int ny);
#define EPSLN (1.0e-10)

int get_contact_index( int size1, int size2, double *x1, double *y1, double *x2, double *y2, double periodx,
		       double periody, int *start1, int *end1, int *start2, int *end2);
int get_overlap_index( double x1, double y1, int nx2, int ny2, const double *x2, const double *y2,
		       int *i2, int *j2 );

int get_align_contact(int tile1, int tile2, int nx1, int ny1, int nx2, int ny2, 
                      const double *x1, const double *y1, const double *x2, 
                      const double *y2, double periodx, double periody,
		      int *istart1, int *iend1, int *jstart1, int *jend1, 
                      int *istart2, int *iend2, int *jstart2, int *jend2)
{
  double *xb1, *yb1, *xb2, *yb2;
  int ncontact, start1, end1, start2, end2;
  
  ncontact = 0;
  /* East bound of tile1 and west bound of tile2 */
  xb1 = east_bound(x1, nx1, ny1);
  yb1 = east_bound(y1, nx1, ny1);
  xb2 = west_bound(x2, nx2, ny2);
  yb2 = west_bound(y2, nx2, ny2);
  if( get_contact_index( ny1, ny2, xb1, yb1, xb2, yb2, periodx, 0.0, &start1, &end1, &start2, &end2) ) {
    istart1[ncontact] = nx1-1;
    iend1[ncontact]   = nx1-1;
    istart2[ncontact] = 1;
    iend2[ncontact]   = 1;    
    jstart1[ncontact] = start1;
    jend1[ncontact]   = end1;
    jstart2[ncontact] = start2;
    jend2[ncontact]   = end2;
    ncontact++;   
  }

  /* East bound of tile1 and SOUTH bound of tile2, tile1 and tile must be different tile */
  if(tile1 != tile2) {
    free(xb2);
    free(yb2);
    xb2 = south_bound(x2, nx2, ny2);
    yb2 = south_bound(y2, nx2, ny2);
    if( get_contact_index( ny1, nx2, xb1, yb1, xb2, yb2, 0.0, 0.0, &start1, &end1, &start2, &end2) ) {
      istart1[ncontact] = nx1-1;
      iend1[ncontact]   = nx1-1;
      istart2[ncontact] = start2;
      iend2[ncontact]   = end2;    
      jstart1[ncontact] = start1;
      jend1[ncontact]   = end1;
      jstart2[ncontact] = 1;
      jend2[ncontact]   = 1;
      ncontact++;
    }
  }
  free(xb1);
  free(yb1);   
  free(xb2);
  free(yb2);  
  
  /* South bound of tile1 and NORTH bound of tile2 */
  xb1 = south_bound(x1, nx1, ny1);
  yb1 = south_bound(y1, nx1, ny1);
  xb2 = north_bound(x2, nx2, ny2);
  yb2 = north_bound(y2, nx2, ny2);
  if( get_contact_index( nx1, nx2, xb1, yb1, xb2, yb2, 0.0, periody, &start1, &end1, &start2, &end2) ) {
    istart1[ncontact] = start1;
    iend1[ncontact]   = end1;
    istart2[ncontact] = start2;
    iend2[ncontact]   = end2;    
    jstart1[ncontact] = 1;
    jend1[ncontact]   = 1;
    jstart2[ncontact] = ny2-1;
    jend2[ncontact]   = ny2-1;
    ncontact++;
  }

  /* South bound of tile1 and East bound of tile2, tile1 and tile must be different tile*/  
  if(tile1 != tile2 ) {
    free(xb2);
    free(yb2);
    xb2 = east_bound(x2, nx2, ny2);
    yb2 = east_bound(y2, nx2, ny2);
    if( get_contact_index( nx1, ny2, xb1, yb1, xb2, yb2, 0.0, 0.0, &start1, &end1, &start2, &end2) ) {
      istart1[ncontact] = start1;
      iend1[ncontact]   = end1;
      istart2[ncontact] = nx2-1;
      iend2[ncontact]   = nx2-1;    
      jstart1[ncontact] = 1;
      jend1[ncontact]   = 1;
      jstart2[ncontact] = start2;
      jend2[ncontact]   = end2;
      ncontact++;
    }
  }
  free(xb1);
  free(yb1);   
  free(xb2);
  free(yb2);  

  /* to avoid duplicate, the following will be done only when tile1 not equal to tile2 */
  if(tile1 != tile2) {
    /* West bound of tile1 and east bound of tile2*/
    xb1 = west_bound(x1, nx1, ny1);
    yb1 = west_bound(y1, nx1, ny1);
    xb2 = east_bound(x2, nx2, ny2);
    yb2 = east_bound(y2, nx2, ny2);
    if( get_contact_index( ny1, ny2, xb1, yb1, xb2, yb2, periodx, 0.0, &start1, &end1, &start2, &end2) ) {
      istart1[ncontact] = 1;
      iend1[ncontact]   = 1;
      istart2[ncontact] = nx2-1;
      iend2[ncontact]   = nx2-1;    
      jstart1[ncontact] = start1;
      jend1[ncontact]   = end1;
      jstart2[ncontact] = start2;
      jend2[ncontact]   = end2;
      ncontact++;
    }
    free(xb2);
    free(yb2);

    /* West bound of tile1 and North bound of tile2 */  
    xb2 = north_bound(x2, nx2, ny2);
    yb2 = north_bound(y2, nx2, ny2);
    if( get_contact_index( ny1, nx2, xb1, yb1, xb2, yb2, 0.0, 0.0, &start1, &end1, &start2, &end2) ) {
      istart1[ncontact] = 1;
      iend1[ncontact]   = 1;
      istart2[ncontact] = start2;
      iend2[ncontact]   = end2;    
      jstart1[ncontact] = start1;
      jend1[ncontact]   = end1;
      jstart2[ncontact] = ny2-1;
      jend2[ncontact]   = ny2-1;
      ncontact++;
    }
    free(xb1);
    free(yb1);   
    free(xb2);
    free(yb2);  

  
    /* North bound of tile1 and South bound of tile2 */
    xb1 = north_bound(x1, nx1, ny1);
    yb1 = north_bound(y1, nx1, ny1);
    xb2 = south_bound(x2, nx2, ny2);
    yb2 = south_bound(y2, nx2, ny2);
    if( get_contact_index( nx1, nx2, xb1, yb1, xb2, yb2, 0.0, periody, &start1, &end1, &start2, &end2) ) {
      istart1[ncontact] = start1;
      iend1[ncontact]   = end1;
      istart2[ncontact] = start2;
      iend2[ncontact]   = end2;    
      jstart1[ncontact] = ny1-1;
      jend1[ncontact]   = ny1-1;
      jstart2[ncontact] = 1;
      jend2[ncontact]   = 1;
      ncontact++;
    }
    free(xb2);
    free(yb2);

    /* North bound of tile1 and West bound of tile2 */  
    xb2 = west_bound(x2, nx2, ny2);
    yb2 = west_bound(y2, nx2, ny2);
    if( get_contact_index( nx1, ny2, xb1, yb1, xb2, yb2, 0.0, 0.0, &start1, &end1, &start2, &end2) ) {
      istart1[ncontact] = start1;
      iend1[ncontact]   = end1;
      istart2[ncontact] = 1;
      iend2[ncontact]   = 1;    
      jstart1[ncontact] = ny1-1;
      jend1[ncontact]   = ny1-1;
      jstart2[ncontact] = start2;
      jend2[ncontact]   = end2;
      ncontact++;
    }
    free(xb1);
    free(yb1);   
    free(xb2);
    free(yb2);  
  }

  /* when tile1 = tile2, we need to consider about folded. Only foled north is considered here */
  if(tile1 == tile2) {
    int i, folded = 1;
    double dx;
    int num;
 
    num = 0;
    xb1 = north_bound(x1, nx1, ny1);
    yb1 = north_bound(y1, nx1, ny1);
    for(i=0; i<nx1/2; i++) {
      if( yb1[i] != yb1[nx1-i-1] ) {
	folded = 0;
	break;
      }
      dx = fabs(xb1[i] - xb1[nx1-i-1]);
      if( dx !=0 && dx != 360 ) {
        if(dx == 180) 
           num++;
        else {
	   folded = 0;
	   break;
        }
      }
    }
    /* if number of pair with lontitude difference == 180 is more than 1, it is not tripolar grid */
    if(num>1) folded = 0;

    if(folded) {
      istart1[ncontact] = 1;
      iend1[ncontact]   = nx1/2;
      istart2[ncontact] = nx1-1;
      iend2[ncontact]   = nx1/2+1;    
      jstart1[ncontact] = ny1-1;
      jend1[ncontact]   = ny1-1;
      jstart2[ncontact] = ny1-1;
      jend2[ncontact]   = ny1-1;
      ncontact++;
    }
    free(xb1);
    free(yb1);
  }
  
  return ncontact;
  
};

/* For simplying reason, we are only deal with nested overlap.
 */

int get_overlap_contact( int tile1, int tile2, int nx1, int ny1, int nx2, int ny2, 
                      const double *x1, const double *y1, const double *x2, 
                      const double *y2, int *istart1, int *iend1, int *jstart1, int *jend1, 
                      int *istart2, int *iend2, int *jstart2, int *jend2)
{
  int l1, l2, count;
  int p1x[4], p1y[4], p2x[4], p2y[4];
  double dx, dy;


  /* check the four corner of tile 2 */
  
  /* southwest corner */
  l2 = 0;
  p2x[0] = 0; p2y[0]=0;
  count = get_overlap_index(x2[l2], y2[l2], nx1, ny1, x1, y1, p1x, p1y);
  /* southeast corner */
  p2x[1] = nx2-1; p2y[1] = 0;
  if( count > 0 ) {
    l2 = nx2-1;
    count = get_overlap_index(x2[l2], y2[l2], nx1, ny1, x1, y1, p1x+1, p1y+1);
  }
  /* northwest corner */
  p2x[2] = 0; p2y[2] = ny2-1;
  if( count > 0 ) {
    l2 = (ny2-1)*nx2;
    count = get_overlap_index(x2[l2], y2[l2], nx1, ny1, x1, y1, p1x+2, p1y+2);
  }      
  /* northeast corner */
  p2x[3] = nx2-1; p2y[3] = ny2-1;
  if( count > 0 ) {
    l2 = (ny2-1)*nx2+ nx2-1;
    count = get_overlap_index(x2[l2], y2[l2], nx1, ny1, x1, y1, p1x+3, p1y+3);
  }  
  /* check the four corner of tile 1 if count==0 */
  if(count == 0) {
  
    /* southwest corner */
    l1 = 0;
    p1x[0] = 0; p1y[0]=0;
    count = get_overlap_index(x1[l1], y1[l1], nx2, ny2, x2, y2, p2x, p2y);
    if(count == 0) return 0;

    /* southeast corner */
    p1x[1] = nx1-1; p1y[1] = 0;
    l1 = nx1-1;
    count = get_overlap_index(x1[l1], y1[l1], nx2, ny2, x2, y2, p2x+1, p2y+1);
    if(count == 0) return 0;

    /* northwest corner */
    p1x[2] = 0; p1y[2] = ny1-1;
    l1 = (ny1-1)*nx1;
    count = get_overlap_index(x1[l1], y1[l1], nx2, ny2, x2, y2, p2x+2, p2y+2);
    if(count == 0) return 0;
    
    /* northeast corner */
    p1x[3] = nx1-1; p1y[3] = ny1-1;
    l1 = (ny1-1)*nx1+ nx1-1;
    count = get_overlap_index(x1[l1], y1[l1], nx2, ny2, x2, y2, p2x+3, p2y+3);
    if(count == 0) return 0;
    
  }    
  
  
/*   for(j1=0; j1<ny1; j1++) for(i1=0;i1<nx1; i1++) { */
/*     l1 = j1*nx1+i1; */
/*     for(j2=0;j2<ny2; j2++) for(i2=0; i2<nx2; i2++) { */
/*       l2 = j2*nx2+i2; */
/*       dx = fabs(x1[l1]- x2[l2]); */
/*       dy = fabs(y1[l1]- y2[l2]); */
/*       if( dx < EPSLN && dy <EPSLN ) { */
/* 	p1x[0] = i1; */
/* 	p1y[0] = j1; */
/* 	p2x[0] = i2; */
/* 	p2y[0] = j2; */
/* 	goto found_sw; */
/*       } */
/*     } */
/*   } */
/*   return 0;    */
/*   found_sw: */
/*   for(j1=0; j1<ny1; j1++) for(i1=nx1-1;i1>=0; i1--) { */
/*     l1 = j1*nx1+i1; */
/*     for(j2=0;j2<ny2; j2++) for(i2=nx2-1; i2>=0; i2--) { */
/*       l2 = j2*nx2+i2; */
/*       dx = fabs(x1[l1]- x2[l2]); */
/*       dy = fabs(y1[l1]- y2[l2]); */
/*       if( dx < EPSLN && dy <EPSLN ) { */
/* 	p1x[1] = i1; */
/* 	p1y[1] = j1; */
/* 	p2x[1] = i2; */
/* 	p2y[1] = j2; */
/* 	goto found_se; */
/*       } */
/*     } */
/*   } */
/*   return 0;    */
/*   found_se: */
  
/*   for(j1=ny1-1; j1>=0; j1--) for(i1=0;i1<nx1; i1++) { */
/*     l1 = j1*nx1+i1; */
/*     for(j2=ny2-1;j2>=0; j2--) for(i2=0; i2<nx2; i2++) { */
/*       l2 = j2*nx2+i2; */
/*       dx = fabs(x1[l1]- x2[l2]); */
/*       dy = fabs(y1[l1]- y2[l2]); */
/*       if( dx < EPSLN && dy <EPSLN ) { */
/* 	p1x[2] = i1; */
/* 	p1y[2] = j1; */
/* 	p2x[2] = i2; */
/* 	p2y[2] = j2; */
/* 	goto found_nw; */
/*       } */
/*     } */
/*   } */
/*   return 0;  */
/*   found_nw: */
  
/*   for(j1=ny1-1; j1>=0; j1--) for(i1=nx1-1;i1>=0; i1--) { */
/*     l1 = j1*nx1+i1; */
/*     for(j2=ny2-1;j2>=0; j2--) for(i2=nx2-1; i2>=0; i2--) { */
/*       l2 = j2*nx2+i2; */
/*       dx = fabs(x1[l1]- x2[l2]); */
/*       dy = fabs(y1[l1]- y2[l2]); */
/*       if( dx < EPSLN && dy <EPSLN ) { */
/* 	p1x[3] = i1; */
/* 	p1y[3] = j1; */
/* 	p2x[3] = i2; */
/* 	p2y[3] = j2; */
/* 	goto found_ne; */
/*       } */
/*     } */
/*   } */
/*   return 0; */
  
/*   found_ne: */
  /* currently we assume there is no rotation for the overlapped grid */
  /* The following must be true
     p1y[0] == p1y[1]   p2y[0] == p2y[1]
     p1y[2] == p1y[3]   p2y[2] == p2y[3]
     p1x[0] == p1x[2]   p2x[0] == p2x[2]
     p1x[1] == p1x[3]   p2x[1] == p2x[3]
  */
  if( p1y[0] != p1y[1] ) return 0;
  if( p2y[0] != p2y[1] ) return 0;
  if( p1x[0] != p1x[2] ) return 0;
  if( p2x[0] != p2x[2] ) return 0;
  if( p1y[2] != p1y[3] ) return 0;
  if( p2y[2] != p2y[3] ) return 0;
  if( p1x[1] != p1x[3] ) return 0;
  if( p2x[1] != p2x[3] ) return 0;

  /* exclude align contact p1y[0] == p1y[2] .or. p1x[0] == p1x[1]  */
  if( p1y[0] == p1y[2] ||  p1x[0] == p1x[1] ) return 0;
  if( p2y[0] == p2y[2] ||  p2x[0] == p2x[1] ) return 0;

  istart1[0] = p1x[0]+1;
  iend1[0]   = p1x[1];
  jstart1[0] = p1y[0]+1;
  jend1[0]   = p1y[3];
  istart2[0] = p2x[0]+1;
  iend2[0]   = p2x[1];
  jstart2[0] = p2y[0]+1;
  jend2[0]   = p2y[3];    
    
  return 1;

}
  
int get_overlap_index( double x1, double y1, int nx2, int ny2, const double *x2, const double *y2,
		       int *i2, int *j2 ) {
  int i, j, l2;
  double dx, dy;
  
  for(j=0; j<ny2; j++) for(i=0;i<nx2; i++) {
    l2 = j*nx2+i;
    dx = fabs(x1- x2[l2]);
    dy = fabs(y1- y2[l2]);
    if( dx < EPSLN && dy <EPSLN ) {
      i2[0] = i;
      j2[0] = j;
      return 1;
    }
  }

  return 0;
}


int get_contact_index( int size1, int size2, double *x1, double *y1, double *x2, double *y2, double periodx,
		       double periody, int *start1, int *end1, int *start2, int *end2)
{
  int i1, i2;
  double dx, dy;

  /* Find the first point in tile 1 cocindent with a point in tile2  */
  *start1 = -1;
  *start2 = -1;
  for(i1=0; i1<size1; i1++) {
    for(i2=0; i2<size2; i2++) {
      dx = fabs(x1[i1]- x2[i2]);
      dx = min(dx, fabs(dx-periodx));
      dy = fabs(y1[i1]- y2[i2]);
      dy = min(dy, fabs(dy-periody));
      if( dx < EPSLN && dy <EPSLN ) {
	*start1 = i1+1;
	*start2 = i2+1;
	goto foundstart;
      }
    }
  }

  return 0;
    
  foundstart:

  /* Find the last point in tile 1 cocindent with a point in tile2 */
  *end1 = -1;
  *end2 = -1;
  for(i1=size1-1; i1>=0; i1--) {
    for(i2=size2-1; i2>=0; i2--) {
      dx = fabs(x1[i1]- x2[i2]);
      dx = min(dx, fabs(dx-periodx));
      dy = fabs(y1[i1]- y2[i2]);
      dy = min(dy, fabs(dy-periody));
      if( dx < EPSLN && dy <EPSLN ) {
	*end1 = i1+1;
	*end2 = i2+1;
	goto foundend;
      }
    }
  }

  return 0;
    
 foundend: if( *start1 == *end1 || *start2 == *end2 ) return 0;

  if(*start1 > *end1 )
    (*start1)--;
  else
    (*end1)--;

  if(*start2 > *end2 )
    (*start2)--;
  else
    (*end2)--;  

  return 1;
    
};


double* west_bound(const double *data, int nx, int ny)
{
  int i;
  double *bnd;
  bnd = (double *)malloc(ny*sizeof(double));
  for(i=0; i<ny; i++) bnd[i] = data[i*nx];
  return bnd;
}

double* east_bound(const double *data, int nx, int ny)
{
  int i;
  double *bnd;
  bnd = (double *)malloc(ny*sizeof(double));
  for(i=0; i<ny; i++) bnd[i] = data[i*nx+nx-1];
  return bnd;
}

double* south_bound(const double *data, int nx, int ny)
{
  int i;
  double *bnd;
  bnd = (double *)malloc(nx*sizeof(double));
  for(i=0; i<nx; i++) bnd[i] = data[i];
  return bnd;
}

double* north_bound(const double *data, int nx, int ny)
{
  int i;
  double *bnd;
  bnd = (double *)malloc(nx*sizeof(double));
  for(i=0; i<nx; i++) bnd[i] = data[(ny-1)*nx+i];
  return bnd;
}
