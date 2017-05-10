#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sched.h>
#include <errno.h>
#include <sys/resource.h>
#include <sys/syscall.h>

static pid_t gettid(void)
{
  return syscall(__NR_gettid);
}

/*
 * Returns this thread's CPU affinity, if bound to a single core,
 * or else -1.
 */
int get_cpu_affinity(void)
{
#ifdef use_libMPI
  cpu_set_t coremask;		/* core affinity mask */

  CPU_ZERO(&coremask);
  if (sched_getaffinity(gettid(),sizeof(cpu_set_t),&coremask) != 0) {
    fprintf(stderr,"Unable to get thread %d affinity. %s\n",gettid(),strerror(errno));
  }

  int cpu;
  int first_cpu = -1;	/* first CPU in range */
  int last_cpu = -1;	/* last CPU in range */
  for (cpu=0;cpu < CPU_SETSIZE;cpu++) {
    if (CPU_ISSET(cpu,&coremask)) {
      if (first_cpu == -1) {
         first_cpu = cpu;
      } else {
        last_cpu = cpu;
      }
    }
  }

  return (last_cpu == -1) ? first_cpu : -1;
#else
  return 0;
#endif

}

int get_cpu_affinity_(void) { return get_cpu_affinity(); }	/* Fortran interface */


/*
 * Set CPU affinity to one core.
 */
void set_cpu_affinity( int cpu )
{
  cpu_set_t coremask;		/* core affinity mask */

  CPU_ZERO(&coremask);
  CPU_SET(cpu,&coremask);
  if (sched_setaffinity(gettid(),sizeof(cpu_set_t),&coremask) != 0) {
    fprintf(stderr,"Unable to set thread %d affinity. %s\n",gettid(),strerror(errno));
  }
}

void set_cpu_affinity_(int *cpu) { set_cpu_affinity(*cpu); }	/* Fortran interface */
