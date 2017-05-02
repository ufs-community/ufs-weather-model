/*
      **** MppDomain.cpp ****
      MppDomain package
      NOTE: only mpi is implemented here. if needed, shmem version
      will be added on in the future.
      Contact: Zhi.Liang@noaa.gov
*/
#include <stdlib.h>
#include <math.h>
#include "mpp.h"
#include "mpp_domain.h"

/***********************************************************
   global variables
***********************************************************/
int pe, npes, root_pe;
#define MAX_BUFFER_SIZE 10000000
double rBuffer[MAX_BUFFER_SIZE];
double sBuffer[MAX_BUFFER_SIZE];

/************************************************************
         void mpp_domain_init()
   initialization routine. get the processor information.
   memory allocation.
************************************************************/
void mpp_domain_init( )
{
  pe      = mpp_pe();
  npes    = mpp_npes();
  root_pe = mpp_root_pe();  

}; /* mpp_domain_init */

/*********************************************************** 
         void mpp_domain_end()
         release memory.
***********************************************************/
void mpp_domain_end ()
{
  /* will add something here if needed */
};

/************************************************************
    void mpp_define_layout()
    define domain layout based on given grid resolution and
    number of processors
***********************************************************/

void mpp_define_layout(int ni, int nj, int ndivs, int layout[])
{
  int idiv, jdiv;
  float fdiv;

  /*first try to divide ndivs in the domain aspect ratio:
    if imperfect aspect, reduce idiv till it divides ndivs */
  fdiv = sqrt((1.0*ndivs*ni)/nj); 
  idiv = floor(fdiv);
  if(fdiv-idiv > 0.5) idiv = idiv + 1;
  idiv = (idiv>1) ? idiv : 1;  /*for isz=1 line above can give 0*/
  while( ndivs%idiv != 0 ) {
    idiv = idiv - 1;
  }                 /*will terminate at idiv=1 if not before*/
  jdiv = ndivs/idiv;
  layout[0] = idiv;
  layout[1] = jdiv;
};

/***************************************************************
 void mpp_compute_extent(int npts, int ndivs, int *ibegin, int *iend)

 Compute extent of 1-D decomposition

 problem of dividing nx points into n domains maintaining symmetry
 i.e nx=18 n=4 4554 and 5445 are solutions but 4455 is not.
 this will always work for nx even n even or odd
 this will always work for nx odd, n odd
 this will never  work for nx odd, n even: for this case we supersede the mirror calculation
 symmetrize = .NOT. ( mod(ndivs,2).EQ.0 .AND. mod(ieg-isg+1,2).EQ.1 )
 nx even n odd fails if n>nx/2

***************************************************************/
void mpp_compute_extent(int npts, int ndivs, int *ibegin, int *iend)
{
  int ndivs_is_odd, npts_is_odd, symmetrize;
  int isg, ieg, is, ie;
  int imax, ndmax, ndmirror;
  int ndiv;
  
  if(ndivs > npts ) {
     mpp_error("mpp_compute_extent: more divisions requested than rows available. " );
  }

  ndivs_is_odd = ndivs%2;
  npts_is_odd = npts%2;
  symmetrize = 0;
  if( ndivs_is_odd && npts_is_odd ) symmetrize = 1; 
  if( ndivs_is_odd == 0 && npts_is_odd == 0 ) symmetrize = 1;
  if( ndivs_is_odd && npts_is_odd == 0 && ndivs < npts/2 ) symmetrize = 1;

  isg = 0;
  ieg = npts-1;
  is = isg;
  for(ndiv=0; ndiv<ndivs; ndiv++){
    /*mirror domains are stored in the list and retrieved if required. */
    if( ndiv == 0 ) { /* initialize max points and max domains */
      imax = ieg;
      ndmax = ndivs;
    }
    /* do bottom half of decomposition, going over the midpoint for odd ndivs */
    if( ndiv < (ndivs-1)/2+1 ) {
      /*domain is sized by dividing remaining points by remaining domains */
      ie = is + ceil((imax-is+1.0)/(ndmax-ndiv) ) - 1;
      ndmirror = (ndivs-1) - ndiv; /* mirror domain */
      if( ndmirror > ndiv && symmetrize ) { /* only for domains over the midpoint */
	/*mirror extents, the max(,) is to eliminate overlaps */
	ibegin[ndmirror] = max( isg+ieg-ie, ie+1 );
	iend[ndmirror]   = max( isg+ieg-is, ie+1 );
	imax = ibegin[ndmirror] - 1;
	ndmax --;
      }
    }
    else {
      if( symmetrize ){
	/*do top half of decomposition by retrieving saved values */
	is = ibegin[ndiv];
	ie = iend[ndiv];
      }
      else {
	ie = is + ceil((imax-is+1.0)/(ndmax-ndiv)) - 1;
      }
    }
 
    ibegin[ndiv] = is;
    iend[ndiv] = ie;
    if( ie < is )mpp_error("mpp_domains(mpp_compute_extent): domain extents must be positive definite." );
    if( ndiv == ndivs-1 && iend[ndiv] != ieg )
      mpp_error( "mpp_domains(mpp_compute_extent): domain extents do not span space completely." );
    is = ie + 1;
  }

}; /* mpp_compute_extent */


/***********************************************************
    void mpp_define_domain_1d(int size, domain1D *domain )
    define 1-D domain decomposition.
**********************************************************/
void mpp_define_domain_1d(int npts, int ndivs, domain1D *domain )
{
  int n, npts_left, pos, size;

  domain->beglist = (int *)malloc(ndivs*sizeof(int));
  domain->endlist = (int *)malloc(ndivs*sizeof(int));

  mpp_compute_extent(npts, ndivs, domain->beglist, domain->endlist);
  
  if(npes == ndivs) {
    domain->start = domain->beglist[pe];
    domain->end   = domain->endlist[pe];
    domain->size  = domain->end - domain->start + 1;
    domain->sizeg = npts;
  }
  
}; /* mpp_define_domain_1d */


/************************************************************
  void define_domain(int ni, int nj, int layout[], int xhalo, int yhalo, domain2D *domain  )
   define 2D domain decomposition
************************************************************/

void mpp_define_domain2d(int ni, int nj, int layout[], int xhalo, int yhalo, domain2D *domain )
{
  domain1D domx, domy;
  int i, j, posx, posy, n; 

  domain->isclist = (int *)malloc(layout[0]*layout[1]*sizeof(int));
  domain->ieclist = (int *)malloc(layout[0]*layout[1]*sizeof(int));
  domain->jsclist = (int *)malloc(layout[0]*layout[1]*sizeof(int));
  domain->jeclist = (int *)malloc(layout[0]*layout[1]*sizeof(int));  

  mpp_define_domain_1d(ni, layout[0], &domx);
  mpp_define_domain_1d(nj, layout[1], &domy); 

  n = 0;
  for(j=0; j<layout[1]; j++) {
    for(i=0; i<layout[0]; i++) {
      domain->isclist[n] = domx.beglist[i]+xhalo;
      domain->ieclist[n] = domx.endlist[i]+xhalo;
      domain->jsclist[n] = domy.beglist[j]+yhalo;
      domain->jeclist[n] = domy.endlist[j]+yhalo;
      n++;
    }
  }

  domain->xhalo = xhalo;
  domain->yhalo = yhalo;
  domain->isc   = domain->isclist[pe];
  domain->iec   = domain->ieclist[pe];
  domain->jsc   = domain->jsclist[pe];
  domain->jec   = domain->jeclist[pe];
  domain->isd   = domain->isc - xhalo;
  domain->ied   = domain->iec + xhalo;
  domain->jsd   = domain->jsc - yhalo;
  domain->jed   = domain->jec + yhalo;  
  domain->nxc   = domain->iec - domain->isc + 1;
  domain->nyc   = domain->jec - domain->jsc + 1;
  domain->nxd   = domain->ied - domain->isd + 1;
  domain->nyd   = domain->jed - domain->jsd + 1;
  domain->nxg   = ni;
  domain->nyg   = nj;

  mpp_delete_domain1d(&domx);
  mpp_delete_domain1d(&domy);
  
}; /* mpp_define_domain2d */

/****************************************************************************
  void mpp_delete_domain1d(domain1D *domain);
  release the memory assigned to 1-D domain
*****************************************************************************/
void mpp_delete_domain1d(domain1D *domain)
{

  free(domain->beglist);
  free(domain->endlist);
  
}; /* mpp_delete_domain1d */

/****************************************************************************
  void mpp_delete_domain2d(domain2D *domain);
  release the memory assigned to 2-D domain
*****************************************************************************/
void mpp_delete_domain2d(domain2D *domain)
{

  free(domain->isclist);
  free(domain->ieclist);
  free(domain->jsclist);
  free(domain->jeclist);
  
}; /* mpp_delete_domain2d */



/***********************************************************
     void get_get_compute_domain(omain2D domain, int *is, int *ie, int *js, int *je)
    get the compute domain decomposition
***********************************************************/

void mpp_get_compute_domain2d(domain2D domain, int *is, int *ie, int *js, int *je)
{
  *is = domain.isc;
  *ie = domain.iec;
  *js = domain.jsc;
  *je = domain.jec;
}; /* mpp_get_compute_domain */

/***********************************************************
     void get_get_compute_domain( int *is, int *ie, int *js, int *je)
    get the compute domain decomposition of current pe list.
***********************************************************/

void mpp_get_compute_domains2d(domain2D domain, int *is, int *ie, int *js, int *je)
{
  int n;
  
  for(n=0; n<npes;n++) {
    is[n] = domain.isclist[n];
    ie[n] = domain.ieclist[n];
    js[n] = domain.jsclist[n];
    je[n] = domain.jeclist[n];
  }

}/* mpp_get_compute_domains */

/**********************************************************
    void mpp_get_data_domain(int *is, int *ie, int *js, int *je)
   get the data domain decomposition
***********************************************************/

void mpp_get_data_domain2d(domain2D domain, int *is, int *ie, int *js, int *je)
{
  *is = domain.isd;
  *ie = domain.ied;
  *js = domain.jsd;
  *je = domain.jed;
}; /* mpp_get_data_domain */ 

/**********************************************************
    void mpp_get_global_domain(int *nx, int *ny )
   get the global domain size
***********************************************************/

void mpp_get_global_domain2d(domain2D domain, int *nx, int *ny )
{
  *nx = domain.nxg;
  *ny = domain.nyg;
}; /* mpp_get_global_domain */ 

/****************************************************************
   void mpp_get_shift(int sizex, int sizey, int *ishift, int *jshift)
   Return the shift value. For non-symmetry domain, ishift will be always 0.
   For symmetric domain, shift value can be 0 or 1 depending on position.

****************************************************************/
void mpp_get_shift(domain2D domain, int sizex, int sizey, int *ishift, int *jshift)
{
  int nxc, nyc, nxd, nyd;
  
  nxc = domain.nxc; nyc = domain.nyc;
  nxd = domain.nxd; nyd = domain.nyd;

  if( (sizex == nxc && sizey == nyc) || (sizex == nxd && sizey == nyd) ) { 
    *ishift = 0;
    *jshift = 0;
  }
  else if ( (sizex == nxc+1 && sizey == nyc) || (sizex == nxd+1 && sizey == nyd)) { 
    *ishift = 1;
    *jshift = 0;
  }
  else if ( (sizex == nxc && sizey == nyc+1) || (sizex == nxd && sizey == nyd+1) ) { 
    *ishift = 0;
    *jshift = 1;
  }
  else if ( (sizex == nxc+1 && sizey == nyc+1) || (sizex == nxd+1 && sizey == nyd+1) ) { 
    *ishift = 1;
    *jshift = 1;
  }
  else
    mpp_error("mpp_domain: data should be on either compute or data domain.");

}; /* mpp_get_shift */

/*************************************************************
    mpp_global_field_all_double(domain2D domain, int sizex, int sizey, const double *ldata, double *gdata)
    get the global data on all the pe.
    ldata is on compute domain and gdata is on global domain
************************************************************/
void mpp_global_field_all_double(domain2D domain, int sizex, int sizey, const double* ldata, double* gdata)
{
  double *send_buffer=NULL, *recv_buffer=NULL;
  int i, j, n, ni, nj, ii, jj, l, p, recv_size;
  int ishift, jshift, nxc, nyc, nxd, nyd, nxg;
  int is, ie, js, je, isd, jsd;
  int send_buffer_is_allocated;

  mpp_get_shift( domain, sizex, sizey, &ishift, &jshift);
  is = domain.isc;
  ie = domain.iec + ishift;
  js = domain.jsc;
  je = domain.jec + jshift;
  
  nxc = ie-is+1;
  nyc = je-js+1;
  isd = domain.isd;
  jsd = domain.jsd;
  nxd = domain.nxd + ishift;
  nyd = domain.nyd + jshift;
  nxg = domain.nxg + ishift;

  /* first fill the send buffer */
  send_buffer = (double *)malloc(nxc*nyc*sizeof(double));
  if( sizex == nxc && sizey == nyc ){ /* data is on compute domain */
    /* for one pe case, just simply copy ldata to gdata */
    if(npes == 1) {
      for(i=0; i<nxc*nyc; i++) gdata[i]=ldata[i];
      return;
    }
    else {
      for(i=0; i<nxc*nyc; i++) send_buffer[i] = ldata[i];
    }
  }
  else if( sizex == nxd && sizey == nyd  ){ /* data is on data domain */
    n = 0;
    /* for one pe case, just simply copy ldata to gdata */
    if(npes == 1) {
      for(j=js;j<=je;j++) for(i=is;i<=ie;i++) {
	gdata[n++] = ldata[(j-jsd)*nxd+(i-isd)];
      }
      return;
    }
    else {
      for(j=js;j<=je;j++) for(i=is;i<=ie;i++) {
	send_buffer[n++] = ldata[(j-jsd)*nxd+(i-isd)];
      }
    }
  }
  else
    mpp_error("mpp_domain(mpp_global_field_all_double: data should be on compute/data domain");
  
  /* send the data */
  for(p=0;p<npes;p++) {
    mpp_send_double(send_buffer, nxc*nyc, p);
  }

  /* received the data */
  for(p=0;p<npes;p++) {
    recv_size = (domain.ieclist[p]-domain.isclist[p]+1+ishift)*(domain.jeclist[p]-domain.jsclist[p]+1+jshift);
    recv_buffer = ( double *) malloc(recv_size*sizeof(double));
    mpp_recv_double(recv_buffer, recv_size, p );
    n = 0;
    for(j=domain.jsclist[p]; j<=domain.jeclist[p]+jshift; j++){
      for(i=domain.isclist[p]; i<=domain.ieclist[p]+ishift; i++){
	gdata[j*nxg+i] = recv_buffer[n++];
      }
    }
    free(recv_buffer);
  }

  mpp_sync_self();
  
  free(send_buffer);
}; /* mpp_global_field_all_double */


/*************************************************************
    mpp_global_field_double(domain2D domain, int sizex, int sizey, const double *ldata, double *gdata)
    get the global data on root pe.
    ldata is on compute domain and gdata is on global domain
************************************************************/
void mpp_global_field_double(domain2D domain, int sizex, int sizey, const double* ldata, double* gdata)
{
  double *send_buffer=NULL, *recv_buffer=NULL;
  int i, j, n, ni, nj, ii, jj, l, p, recv_size;
  int ishift, jshift, nxc, nyc, nxd, nyd, nxg;
  int is, ie, js, je, isd, jsd;
  
  mpp_get_shift( domain, sizex, sizey, &ishift, &jshift);
  is = domain.isc;
  ie = domain.iec + ishift;
  js = domain.jsc;
  je = domain.jec + jshift;
  
  nxc = ie-is+1;
  nyc = je-js+1;
  isd = domain.isd;
  jsd = domain.jsd;
  nxd = domain.nxd + ishift;
  nyd = domain.nyd + jshift;
  nxg = domain.nxg + ishift;
  
  /* all other pe except root pe will send data to root pe */
  if( pe != root_pe) {
    if( sizex == nxc && sizey == nyc ){ /* data is on compute domain */
      mpp_send_double(ldata, sizex*sizey, root_pe);
    }
    else if( sizex == nxd && sizey == nyd  ){ /* data is on data domain */
      send_buffer = (double *)malloc(nxc*nyc*sizeof(double));
      n = 0;
      for(j=js;j<=je;j++) {
	for(i=is;i<=ie;i++) send_buffer[n++] = ldata[(j-jsd)*nxd+(i-isd)];
      }
      mpp_send_double(send_buffer, nxc*nyc, root_pe);
    }
    else
      mpp_error("mpp_domain: data should be on compute/data domain");
  }

  /* receive from other pe on root pe   */
  if( pe == root_pe ) {
    for(p=0;p<npes;p++) {
      if( p == root_pe) {
	if( sizex == nxc && sizey == nyc  ){ /* data is on compute domain */
	  n = 0;
	  for(j=js;j<=je;j++) {
	    for(i=is;i<=ie;i++) gdata[j*nxg+i] = ldata[n++];
	  }
	}
	else {
	  for(j=js;j<=je;j++) {
	    for(i=is;i<=ie;i++) gdata[j*nxg+i] = ldata[(j-jsd)*nxd+(i-isd)];
	  }
	}
      }
      else {
	recv_size = (domain.ieclist[p]-domain.isclist[p]+1+ishift)*(domain.jeclist[p]-domain.jsclist[p]+1+jshift);
	recv_buffer = ( double *) malloc(recv_size*sizeof(double));
	mpp_recv_double(recv_buffer, recv_size, p );
	n = 0;
	for(j=domain.jsclist[p]; j<=domain.jeclist[p]+jshift; j++){
	  for(i=domain.isclist[p]; i<=domain.ieclist[p]+ishift; i++){
	    gdata[j*nxg+i] = recv_buffer[n++];
	  }
	}
	free(recv_buffer);
      }
    }
  }

  mpp_sync_self();
  
  if(send_buffer != NULL) free(send_buffer);
}; /* mpp_global_field_double */

/*************************************************************
    mpp_global_field_int(domain2D domain, int sizex, int sizey, const int *ldata, int *gdata)
    get the global data on root pe.
    ldata is on compute domain and gdata is on global domain
************************************************************/
void mpp_global_field_int(domain2D domain, int sizex, int sizey, const int* ldata, int* gdata)
{
  int *send_buffer=NULL, *recv_buffer=NULL;
  int i, j, n, ni, nj, ii, jj, l, p, recv_size;
  int ishift, jshift, nxc, nyc, nxd, nyd, nxg;
  int is, ie, js, je, isd, jsd;
  
  mpp_get_shift( domain, sizex, sizey, &ishift, &jshift);
  is = domain.isc;
  ie = domain.iec + ishift;
  js = domain.jsc;
  je = domain.jec + jshift;
  
  nxc = ie-is+1;
  nyc = je-js+1;
  isd = domain.isd;
  jsd = domain.jsd;
  nxd = domain.nxd + ishift;
  nyd = domain.nyd + jshift;
  nxg = domain.nxg + ishift;
  
  /* all other pe except root pe will send data to root pe */
  if( pe != root_pe) {
    if( sizex == nxc && sizey == nyc ){ /* data is on compute domain */
      mpp_send_int(ldata, sizex*sizey, root_pe);
    }
    else if( sizex == nxd && sizey == nyd  ){ /* data is on data domain */
      send_buffer = (int *)malloc(nxc*nyc*sizeof(int));
      n = 0;
      for(j=js;j<=je;j++) {
	for(i=is;i<=ie;i++) send_buffer[n++] = ldata[(j-jsd)*nxd+(i-isd)];
      }
      mpp_send_int(send_buffer, nxc*nyc, root_pe);
    }
    else
      mpp_error("mpp_domain: data should be on compute/data domain");
  }

  /* receive from other pe on root pe   */
  if( pe == root_pe ) {
    for(p=0;p<npes;p++) {
      if( p == root_pe) {
	if( sizex == nxc && sizey == nyc  ){ /* data is on compute domain */
	  n = 0;
	  for(j=js;j<=je;j++) {
	    for(i=is;i<=ie;i++) gdata[j*nxg+i] = ldata[n++];
	  }
	}
	else {
	  for(j=js;j<=je;j++) {
	    for(i=is;i<=ie;i++) gdata[j*nxg+i] = ldata[(j-jsd)*nxd+(i-isd)];
	  }
	}
      }
      else {
	recv_size = (domain.ieclist[p]-domain.isclist[p]+1+ishift)*(domain.jeclist[p]-domain.jsclist[p]+1+jshift);
	recv_buffer = ( int *) malloc(recv_size*sizeof(int));
	mpp_recv_int(recv_buffer, recv_size, p );
	n = 0;
	for(j=domain.jsclist[p]; j<=domain.jeclist[p]+jshift; j++){
	  for(i=domain.isclist[p]; i<=domain.ieclist[p]+ishift; i++){
	    gdata[j*nxg+i] = recv_buffer[n++];
	  }
	}
	free(recv_buffer);
      }
    }
  }

  mpp_sync_self();
  
  if(send_buffer != NULL) free(send_buffer);
}; /* mpp_global_field_int */

/*************************************************************
    mpp_global_field_double_3D(domain2D domain, int sizex, int sizey, int sizez,
                               const double *ldata, double *gdata)
    get the global data on root pe.
    ldata is on compute domain and gdata is on global domain
************************************************************/
void mpp_global_field_double_3D(domain2D domain, int sizex, int sizey, int sizez,
				const double* ldata, double* gdata)
{
  double *send_buffer=NULL, *recv_buffer=NULL;
  int i, j, k, n, ni, nj, ii, jj, l, p, recv_size;
  int ishift, jshift, nxc, nyc, nxd, nyd, nxg, nyg;
  int is, ie, js, je, isd, jsd;
  int send_size;  

  mpp_get_shift( domain, sizex, sizey, &ishift, &jshift);
  is = domain.isc;
  ie = domain.iec + ishift;
  js = domain.jsc;
  je = domain.jec + jshift;
  
  nxc = ie-is+1;
  nyc = je-js+1;
  isd = domain.isd;
  jsd = domain.jsd;
  nxd = domain.nxd + ishift;
  nyd = domain.nyd + jshift;
  nxg = domain.nxg + ishift;
  nyg = domain.nyg + ishift;

  /* all other pe except root pe will send data to root pe */
  send_size = nxc*nyc*sizez;
  if( pe != root_pe) {
    if( sizex == nxc && sizey == nyc ){ /* data is on compute domain */
      mpp_send_double(ldata, sizex*sizey*sizez, root_pe);
    }
    else if( sizex == nxd && sizey == nyd  ){ /* data is on data domain */
      if( send_size > MAX_BUFFER_SIZE) {
         send_buffer = (double *)malloc(send_size*sizeof(double));
      }
      else {
         send_buffer = sBuffer;
      }
      n = 0;
      for(k=0; k<sizez; k++) for(j=js;j<=je;j++) {
	for(i=is;i<=ie;i++) send_buffer[n++] = ldata[k*nxd*nyd+(j-jsd)*nxd+(i-isd)];
      }
      mpp_send_double(send_buffer, send_size, root_pe);
    }
    else
      mpp_error("mpp_domain: data should be on compute/data domain");
  }
  /* receive from other pe on root pe   */
  if( pe == root_pe ) {
    for(p=0;p<npes;p++) {
      if( p == root_pe) {
	if( sizex == nxc && sizey == nyc  ){ /* data is on compute domain */
	  n = 0;
	  for(k=0;k<sizez;k++) for(j=js;j<=je;j++) {
	    for(i=is;i<=ie;i++) gdata[k*nxg*nyg+j*nxg+i] = ldata[n++];
	  }
	}
	else {
	  for(k=0;k<sizez;k++) for(j=js;j<=je;j++) {
	    for(i=is;i<=ie;i++) gdata[k*nxg*nyg+j*nxg+i] = ldata[k*nxd*nyd+(j-jsd)*nxd+(i-isd)];
	  }
	}
      }
      else {
	recv_size = (domain.ieclist[p]-domain.isclist[p]+1+ishift)
	  *(domain.jeclist[p]-domain.jsclist[p]+1+jshift)*sizez;
        if(recv_size>MAX_BUFFER_SIZE){
	   recv_buffer = ( double *) malloc(recv_size*sizeof(double));
        }
        else {
           recv_buffer = rBuffer;
        }
	mpp_recv_double(recv_buffer, recv_size, p );
	n = 0;
	for(k=0;k<sizez;k++) {
	  for(j=domain.jsclist[p]; j<=domain.jeclist[p]+jshift; j++){
	    for(i=domain.isclist[p]; i<=domain.ieclist[p]+ishift; i++){
	      gdata[k*nxg*nyg+j*nxg+i] = recv_buffer[n++];
	    }
	  }
	}
        if(recv_size>MAX_BUFFER_SIZE)
	   free(recv_buffer);
        else
           recv_buffer=NULL;
      }
    }
  }
  mpp_sync_self();
/*z1l: mpp_sync is needed when running on multiple processor job. Otherwisde the memory
usage will increase. For example, remap_land will fail when running on 270 processors */
  mpp_sync();

  if(send_buffer != NULL) {
     if(send_size>MAX_BUFFER_SIZE)
        free(send_buffer);
     else
        send_buffer = NULL;
  }
}; /* mpp_global_field_double */

/*******************************************************************************
  void mpp_gather_field_int(int lsize, int *ldata, int *gdata)
  gather integer data onto every processor
*******************************************************************************/
void mpp_gather_field_int(int lsize, int *ldata, int *gdata)
{
  int n, p, i;
  int *rbuffer=NULL;
  int *rsize=NULL;

  rsize = (int *)malloc(npes*sizeof(int));

  
  for(p = 0; p<npes; p++) {
    if(pe != p) { /* send to other pe. */
      mpp_send_int(&lsize, 1, p);
    }
  }

  for(p = 0; p<npes; p++) {
    if(pe != p) { /* recv from other pe. */
      mpp_recv_int(rsize+p, 1, p);
    }
  }

  mpp_sync_self();

  for(p = 0; p<npes; p++) {
    if(pe != p) { /* send to other pe. */
      if(lsize>0) mpp_send_int(ldata, lsize, p);
    }
  }  
  n = 0;
  /* receive from other pe and fill the gdata */
  for(p = 0; p<npes; p++) {
    if(pe != p) { /* recv from other pe. */
      if(rsize[p]>0) {
	rbuffer = ( int *) malloc(rsize[p]*sizeof(int));
	mpp_recv_int(rbuffer, rsize[p], p );
	for(i=0; i<rsize[p]; i++) gdata[n++] = rbuffer[i];
	free(rbuffer);
      }
    }
    else {
      for(i=0; i<lsize; i++) gdata[n++] = ldata[i];
    }
  }

  mpp_sync_self();
  free(rsize);
  
}; /* mpp_gather_field_int */

/*******************************************************************************
  void mpp_gather_field_int_root(int lsize, int *ldata, int *gdata)
  gather integer data onto root processor
*******************************************************************************/
void mpp_gather_field_int_root(int lsize, int *ldata, int *gdata)
{
  int n, p, i;
  int *rbuffer=NULL;
  int *rsize=NULL;

  rsize = (int *)malloc(npes*sizeof(int));
  
  /* all other pe except root pe will send data to root pe */
  if( pe != root_pe) {
      mpp_send_int(&lsize, 1, root_pe);
  }

  else {
    for(p = 0; p<npes; p++) {
       if(root_pe != p) mpp_recv_int(rsize+p, 1, p);
    }
  }

  mpp_sync_self();

  if( pe != root_pe) {
      if(lsize>0) mpp_send_int(ldata, lsize, root_pe);
  }  
  else {
    int cur_size;
    n = 0;
    cur_size = 0;
    /* receive from other pe and fill the gdata */
    for(p = 0; p<npes; p++) {
      if(root_pe != p) { /* recv from other pe. */
	if(rsize[p]>0) {
	  if( rsize[p] > cur_size ) {
	    if( rbuffer ) free(rbuffer);
	    rbuffer = ( int *) malloc(rsize[p]*sizeof(int));
	    cur_size = rsize[p];
	  }
	  mpp_recv_int(rbuffer, rsize[p], p );
	  for(i=0; i<rsize[p]; i++) gdata[n++] = rbuffer[i];
	}
      }
      else {
	for(i=0; i<lsize; i++) gdata[n++] = ldata[i];
      }
    }
    if(rbuffer) free(rbuffer);
  }

  mpp_sync_self();
  free(rsize);
  
}; /* mpp_gather_field_int_root */


/*******************************************************************************
  void mpp_gather_field_double_root(int lsize, double *ldata, double *gdata)
  gather double data onto root processor
*******************************************************************************/
void mpp_gather_field_double_root(int lsize, double *ldata, double *gdata)
{
  int n, p, i;
  double *rbuffer=NULL;
  int *rsize=NULL;

  rsize = (int *)malloc(npes*sizeof(int));
  
  /* all other pe except root pe will send data to root pe */
  if( pe != root_pe) {
      mpp_send_int(&lsize, 1, root_pe);
  }

  else {
    for(p = 0; p<npes; p++) {
       if(root_pe != p) mpp_recv_int(rsize+p, 1, p);
    }
  }

  mpp_sync_self();

  if( pe != root_pe) {
      if(lsize>0) mpp_send_double(ldata, lsize, root_pe);
  }  
  else {
    int cur_size;
    n = 0;
    cur_size = 0;
    /* receive from other pe and fill the gdata */
    for(p = 0; p<npes; p++) {
      if(root_pe != p) { /* recv from other pe. */
	if(rsize[p]>0) {
	  if( rsize[p] > cur_size ) {
	    if( rbuffer ) free(rbuffer);
	    rbuffer = ( double *) malloc(rsize[p]*sizeof(double));
	    cur_size = rsize[p];
	  }
	  mpp_recv_double(rbuffer, rsize[p], p );
	  for(i=0; i<rsize[p]; i++) gdata[n++] = rbuffer[i];
	}
      }
      else {
	for(i=0; i<lsize; i++) gdata[n++] = ldata[i];
      }
    }
    if(rbuffer) free(rbuffer);
  }

  mpp_sync_self();
  free(rsize);
  
}; /* mpp_gather_field_double_root */


/*******************************************************************************
  void mpp_gather_field_double(int lsize, double *ldata, double *gdata)
  gather doubleeger data onto every processor
*******************************************************************************/
void mpp_gather_field_double(int lsize, double *ldata, double *gdata)
{
  int n, p, i;
  double*rbuffer=NULL;
  int *rsize=NULL;

  rsize = (int *)malloc(npes*sizeof(int));
  
  for(p = 0; p<npes; p++) {
    if(pe != p) { /* send to other pe. */
      mpp_send_int(&lsize, 1, p);
    }
  }

  for(p = 0; p<npes; p++) {
    if(pe != p) { /* recv from other pe. */
      mpp_recv_int(rsize+p, 1, p);
    }
  }

  mpp_sync_self();

  for(p = 0; p<npes; p++) {
    if(pe != p) { /* send to other pe. */
      if(lsize>0) mpp_send_double(ldata, lsize, p);
    }
  }  
  n = 0;
  /* receive from other pe and fill the gdata */
  for(p = 0; p<npes; p++) {
    if(pe != p) { /* recv from other pe. */
      if(rsize[p]>0) {
	rbuffer = (double *) malloc(rsize[p]*sizeof(double));
	mpp_recv_double(rbuffer, rsize[p], p );
	for(i=0; i<rsize[p]; i++) gdata[n++] = rbuffer[i];
	free(rbuffer);
      }
    }
    else {
      for(i=0; i<lsize; i++) gdata[n++] = ldata[i];
    }
  }

  mpp_sync_self();
  free(rsize);
  
}; /* mpp_gather_field_double*/

