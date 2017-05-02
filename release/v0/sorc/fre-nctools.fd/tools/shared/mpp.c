#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>
#ifdef use_libMPI
#include <mpi.h>
#endif
#include "mpp.h"


/****************************************************
         global variables
*****************************************************/
int npes, root_pe, pe;
int *pelist=NULL;
const int tag = 1;
#ifdef use_libMPI  
MPI_Request *request;
#endif

/**************************************************************
                     void mpp_init()
     this routine will create communicator.
***************************************************************/

void mpp_init(int *argc, char ***argv)
{
  int n;
  
#ifdef use_libMPI
  MPI_Init(argc, argv); 
  MPI_Comm_rank(MPI_COMM_WORLD,&pe);
  MPI_Comm_size(MPI_COMM_WORLD,&npes);
  request = (MPI_Request *)malloc(npes*sizeof(MPI_Request));
  for(n=0; n<npes; n++) request[n] = MPI_REQUEST_NULL;
#else
  pe = 0;
  npes = 1;
#endif
  pelist = (int *)malloc(npes*sizeof(int));
  for(n=0; n<npes; n++) pelist[n] = n;
  root_pe = 0;
}; /* mpp_init */

/***********************************************************
               void mpp_end()
     This routine will terminate the parallel.
************************************************************/

void mpp_end()
{
#ifdef use_libMPI   
  MPI_Finalize();
#endif  
} /* mpp_end */

/*****************************************************************
           int mpp_pe()
      Returns processor ID.
******************************************************************/

int mpp_pe()
{
  return pe;
}; /* mpp_pe */


/**************************************************************
                 int mpp_npes()
      Returns processor count for current pelist.
**************************************************************/

int mpp_npes()
{
  return npes;
}; /* mpp_npes */

/*************************************************************
               int mpp_root_pe()
    return root processor of current pelist
*************************************************************/

int mpp_root_pe()
{
  return root_pe;
}; /* mpp_root_pe */

/*************************************************************
               int* mpp_get_pelist()
    return current pelist
*************************************************************/

int* mpp_get_pelist()
{
  return pelist;
}; /* mpp_get_pelist */

/************************************************************
    void mpp_sync_self()
 this is to check if current PE's outstanding puts are complete
*************************************************************/
void mpp_sync_self() {
  int n;
#ifdef use_libMPI     
  MPI_Status status;
  
  for(n=0; n<npes; n++) {
    if(request[n] != MPI_REQUEST_NULL) MPI_Wait( request+n, &status );
  }
#endif
  
}

/************************************************************
  void mpp_sync()
************************************************************/
void mpp_sync()
{
#ifdef use_libMPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

}

/*************************************************************
    void mpp_send_double(const double* data, int size, int to_pe)
      send data to "to_pe"
*************************************************************/

void mpp_send_double(const double* data, int size, int to_pe)
{
#ifdef use_libMPI    
  MPI_Status status;
  /* make sure only one message from pe->to_pe in queue */
  if(request[to_pe] != MPI_REQUEST_NULL) {
    MPI_Wait( request+to_pe, &status );
  }
    
  MPI_Isend(data, size, MPI_DOUBLE, to_pe, tag, MPI_COMM_WORLD, request+to_pe);
#endif
  
}; /* mpp_send_double */


/*************************************************************
    void mpp_send_int(const int* data, int size, int to_pe)
      send data to "to_pe"
*************************************************************/

void mpp_send_int(const int* data, int size, int to_pe)
{
#ifdef use_libMPI    
  MPI_Status status;
  if(request[to_pe] != MPI_REQUEST_NULL) {
    MPI_Wait( request+to_pe, &status );
  }  

  MPI_Isend(data, size, MPI_INT, to_pe, tag, MPI_COMM_WORLD, request+to_pe);
#endif
  
}; /* mpp_send_int */

/***********************************************************
    void mpp_recv_double(double* data, int size, int from_pe)
     receive data from "from_pe"
***********************************************************/

void mpp_recv_double(double* data, int size, int from_pe)
{
#ifdef use_libMPI      
  MPI_Status status;
  MPI_Recv(data, size, MPI_DOUBLE, from_pe, MPI_ANY_TAG, MPI_COMM_WORLD,&status);
#endif  
}; /* mpp_recv_double */

/***********************************************************
    void mpp_recv_int(int* data, int size, int from_pe)
     receive data from "from_pe"
***********************************************************/

void mpp_recv_int(int* data, int size, int from_pe)
{
#ifdef use_libMPI      
  MPI_Status status;
  MPI_Recv(data, size, MPI_INT, from_pe, MPI_ANY_TAG, MPI_COMM_WORLD,&status);
#endif  
}; /* mpp_recv_int */


/*******************************************************************************
  int mpp_sum_int(int count, int *data)
  sum integer over all the pes.
*******************************************************************************/
void mpp_sum_int(int count, int *data)
{

#ifdef use_libMPI
  int i;
  int *sum;
  sum = (int *)malloc(count*sizeof(int));
  MPI_Allreduce(data, sum, count, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  for(i=0; i<count; i++)data[i] = sum[i];
  free(sum);
#endif


}; /* mpp_sum_int */

/*******************************************************************************
  void mpp_sum_double(int count, double *data)
  sum double over all the pes.
*******************************************************************************/
void mpp_sum_double(int count, double *data)
{

#ifdef use_libMPI
  int i;
  double *sum;
  sum = (double *)malloc(count*sizeof(double));
  MPI_Allreduce(data, sum, count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for(i=0; i<count; i++)data[i] = sum[i];
  free(sum);  
#endif


}; /* mpp_sum_double */

void mpp_min_double(int count, double *data)
{

#ifdef use_libMPI
  int i;
  double *minval;
  minval = (double *)malloc(count*sizeof(double));
  MPI_Allreduce(data, minval, count, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  for(i=0; i<count; i++) data[i] = minval[i];
  free(minval);
#endif

}; /* mpp_min_double */


void mpp_max_double(int count, double *data)
{

#ifdef use_libMPI
  int i;
  double *maxval;
  maxval = (double *)malloc(count*sizeof(double));
  MPI_Allreduce(data, maxval, count, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  for(i=0; i<count; i++) data[i] = maxval[i];
  free(maxval);
#endif

}; /* mpp_max_double */

/***********************************************************
    void mpp_error(char *str)
    error handler: will print out error message and then abort
***********************************************************/

void mpp_error(char *str)
{
  fprintf(stderr, "Error from pe %d: %s\n", pe, str );
#ifdef use_libMPI      
  MPI_Abort(MPI_COMM_WORLD, -1);
#else
  exit(1);
#endif  
}; /* mpp_error */

double get_mem_usage(void)

{
  double mem;
  
#if defined(__sgi) || defined(__aix) || defined(__SX)
#define RUSAGE_SELF      0         /* calling process */
#define RUSAGE_CHILDREN  -1        /* terminated child processes */
 struct rusage my_rusage;
 int iret;

 my_rusage.ru_maxrss = 0;
 iret = getrusage(RUSAGE_SELF,&my_rusage);
 mem = my_rusage.ru_maxrss;
 mem /= 1000;

#else
 char filename[]="/proc/self/status";
 char mesg[256];
 FILE *fp;
 char line[80], units[32];

 fp = fopen(filename, "r");
 if(!fp) {
   strcpy(mesg, "tool_util.c: Can not open ascii file ");
   strcat(mesg,filename);
   mpp_error(mesg);
 }

 mem = 0;
 while(fgets(line, 80, fp) != NULL)  /* get a line, up to 80 chars from fp.  done if NULL */
   {
     if( strncmp(line, "VmHWM:", 6) == 0) {
       sscanf(line+6, "%lf %s", &mem, units);
       if( strcmp(units, "kB") == 0) mem = mem/1024;
       break;
     }
   }
 fclose(fp); 
#endif

 return mem;
 
}

void print_time(const char* text, double t)
{
  double tmin, tmax, tavg;
  
  tmin=t;
  tmax=t;
  tavg=t;
  mpp_min_double(1, &tmin);
  mpp_max_double(1, &tmax);
  mpp_sum_double(1, &tavg);
  tavg /= mpp_npes();
  if( mpp_pe() == mpp_root_pe() ) {
    printf("Running time for %s, min=%g, max=%g, avg=%g\n", text, tmin, tmax, tavg); 
  }

}


void print_mem_usage(const char* text)
{
  double m, mmin, mmax, mavg;


    m = get_mem_usage();
    mmin = m;
    mmax = m;
    mavg = m;
    mpp_min_double(1, &mmin);
    mpp_max_double(1, &mmax);
    mpp_sum_double(1, &mavg);
    mavg /= mpp_npes();
    if( mpp_pe() == mpp_root_pe() ) {
      printf("Memuse(MB) at %s, min=%g, max=%g, avg=%g\n", text, mmin, mmax, mavg); 
    }
}
      
