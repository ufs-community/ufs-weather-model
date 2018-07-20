/*
 * This work (Common Community Physics Package), identified by NOAA, NCAR,
 * CU/CIRES, is free of known copyright restrictions and is placed in the
 * public domain.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @file ccpp_memory.c
 *
 * Routines for memory profiling for different architectures..
 *
 * @ingroup CCPP
 * @{
 **/

#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "unistd.h"
#include "sys/types.h"
#if __unix__
#include "sys/sysinfo.h"
#elif __APPLE__
#include <sys/sysctl.h>
#include <sys/mount.h>
#include <mach/mach.h>
#endif
#ifdef MPI
#include <mpi.h>
#endif

#include "ccpp_memory.h"

#if __unix__

/**
 * Parses an output line of /proc/self/status and returns an integer
 * value. Assumes that a digit will be found and the line ends in " Kb".
 *
 * @param[in]  line  The line to parse (from /proc/self/status)
 * @retval     i     Integer value of digits parsed in line
 **/
static int parseLine(char* line){
	int i = strlen(line);
	const char* p = line;
	while (*p <'0' || *p > '9') p++;
	line[i-3] = '\0';
	i = atoi(p);
	return i;
}

/**
 * Parses /proc/self/status and returns virtual memory per process.
 * The value returned is in kB (i.e. same units as /proc/self/status).
 *
 * @retval  result  Virtual memory used by this process in kB
 **/
static int getVirtMemPerProcess(){
	FILE* file = fopen("/proc/self/status", "r");
	int result = -1;
	char line[128];

	while (fgets(line, 128, file) != NULL){
		if (strncmp(line, "VmSize:", 7) == 0){
			result = parseLine(line);
			break;
		}
	}
	fclose(file);
	return result;
}

/**
 * Parses /proc/self/status and returns physical memory per process.
 * The value returned is in kB (i.e. same units as /proc/self/status).
 *
 * @retval  result  Physical memory used by this process in kB
 **/
static int getPhysMemPerProcess(){
	FILE* file = fopen("/proc/self/status", "r");
	int result = -1;
	char line[128];

	while (fgets(line, 128, file) != NULL){
		if (strncmp(line, "VmRSS:", 6) == 0){
			result = parseLine(line);
			break;
		}
	}
	fclose(file);
	return result;
}

/**
 * Parses /proc/self/status and returns maximum resident size (physical memory)
 * per process. The value returned is in kB (same units as /proc/self/status).
 *
 * @retval  result  Maximum resident size used by this process in kB
 **/
static int getPhysMemPerProcessMax(){
	FILE* file = fopen("/proc/self/status", "r");
	int result = -1;
	char line[128];

	while (fgets(line, 128, file) != NULL){
		if (strncmp(line, "VmHWM:", 6) == 0){
			result = parseLine(line);
			break;
		}
	}
	fclose(file);
	return result;
}

/**
 * Extended memory statistics on Linux systems (multiple nodes).
 * Diagnostics are returned as formatted output (assumed to be used
 * for model development only, since seriously impacting runtime).
 *
 * @param[in]  mpicomm  MPI communicator
 * @param[in]  str      Formatted output string
 * @param[in]  lstr     Length of str
 * @retval     0        Successful
 * @retval     1        Not successful
 **/
int ccpp_memory_usage_c(const int mpicomm, char* str, int lstr){

	// Error handling
	int ierr;
	ierr = 0;

	// Local string with correct length
	char strtmp[lstr];

	// Get MPI rank and size
	int mpirank;
	int mpisize;
	int mpiroot;
#ifdef MPI
	ierr = MPI_Comm_rank(mpicomm, &mpirank);
	if (ierr!=MPI_SUCCESS) return 1;
	ierr = MPI_Comm_size(mpicomm, &mpisize);
	if (ierr!=MPI_SUCCESS) return 1;
#else
	mpirank = 0;
	mpisize = 1;
#endif

	// Retrieve hostname
	char hostname[128];
	gethostname(hostname, sizeof(hostname));

	// Retrieve PID
	int pid = (int)getpid();

	// Retrieve memory statistics

	struct sysinfo memInfo;
	sysinfo (&memInfo);

	// Available virtual memory this node - in Bytes
	long long totalVirtMem = memInfo.totalram;
	totalVirtMem += memInfo.totalswap;
	totalVirtMem *= memInfo.mem_unit;
	totalVirtMem = (double)totalVirtMem/1024.0/1024.0; // convert to MB

	// Virtual memory used this node - in Bytes
	long long virtMemUsed = memInfo.totalram - memInfo.freeram;
	virtMemUsed += memInfo.totalswap - memInfo.freeswap;
	virtMemUsed *= memInfo.mem_unit;
	virtMemUsed = (double)virtMemUsed/1024.0/1024.0; // convert to MB

	// Available physical memory this node - in Bytes
	long long totalPhysMem = memInfo.totalram;
	totalPhysMem *= memInfo.mem_unit;
	totalPhysMem = (double)totalPhysMem/1024.0/1024.0; // convert to MB

	// Physical memory used this node - in Bytes
	long long physMemUsed = memInfo.totalram - memInfo.freeram;
	physMemUsed *= memInfo.mem_unit;
	physMemUsed = (double)physMemUsed/1024.0/1024.0; // convert to MB

	// Memory statistics per process - in kB
	int virtMemPerProcess = (double)getVirtMemPerProcess()/1024.0; // convert to MB
	int physMemPerProcess = (double)getPhysMemPerProcess()/1024.0; // convert to MB
	int physMemPerProcessMax = (double)getPhysMemPerProcessMax()/1024.0; // convert to MB

#ifdef MPI
	// Average across MPI ranks
	int virtMemPerProcessAvg;
	int physMemPerProcessAvg;
	int physMemPerProcessMaxAvg;
	ierr = MPI_Allreduce(&virtMemPerProcess, &virtMemPerProcessAvg, 1, MPI_INT, MPI_SUM, mpicomm);
	if (ierr!=MPI_SUCCESS) return 1;
	virtMemPerProcessAvg = (int)((float)virtMemPerProcessAvg/(float)mpisize);
	ierr = MPI_Allreduce(&physMemPerProcess, &physMemPerProcessAvg, 1, MPI_INT, MPI_SUM, mpicomm);
	if (ierr!=MPI_SUCCESS) return 1;
	physMemPerProcessAvg = (int)((float)physMemPerProcessAvg/(float)mpisize);
	ierr = MPI_Allreduce(&physMemPerProcessMax, &physMemPerProcessMaxAvg, 1, MPI_INT, MPI_SUM, mpicomm);
	if (ierr!=MPI_SUCCESS) return 1;
	physMemPerProcessMaxAvg = (int)((float)physMemPerProcessMaxAvg/(float)mpisize);
#endif

	// Create formatted output
#ifndef MPI
	snprintf(strtmp, sizeof(strtmp), 
"--------------------------------------------------------------------------------\n\
Memory usage - MPI rank %d/%d\n\
\n\
Statistics for node %s:\n\
Total virtual  memory:      %lld MB\n\
Total virtual  memory used: %lld MB\n\
Total physical memory:      %lld MB\n\
Total physical memory used: %lld MB\n\
\n\
Statistics for process %d:\n\
Virtual    memory this process: %d MB\n\
Physical   memory this process: %d MB\n\
Max. phys. memory this process: %d MB\n\
--------------------------------------------------------------------------------\n\
", mpirank, mpisize, hostname, totalVirtMem, virtMemUsed, totalPhysMem, physMemUsed,
pid, virtMemPerProcess, physMemPerProcess, physMemPerProcessMax);
#else
	snprintf(strtmp, sizeof(strtmp), 
"--------------------------------------------------------------------------------\n\
Memory usage - MPI rank %d/%d\n\
\n\
Statistics for node %s:\n\
Total virtual  memory:      %lld MB\n\
Total virtual  memory used: %lld MB\n\
Total physical memory:      %lld MB\n\
Total physical memory used: %lld MB\n\
\n\
Statistics for process %d:\n\
Virtual    memory this process: %d MB\n\
Physical   memory this process: %d MB\n\
Max. phys. memory this process: %d MB\n\
\n\
Statistics across all processes:\n\
Virtual    memory average: %d MB\n\
Physical   memory average: %d MB\n\
Max. phys. memory average: %d MB\n\
--------------------------------------------------------------------------------\n\
", mpirank, mpisize, hostname, totalVirtMem, virtMemUsed, totalPhysMem, physMemUsed,
pid, virtMemPerProcess, physMemPerProcess, physMemPerProcessMax, virtMemPerProcessAvg,
physMemPerProcessAvg, physMemPerProcessMaxAvg);
#endif

	// Copy local string to output string (same size!)
	strcpy(str,strtmp);

	return ierr;
}

#elif __APPLE__

/**
 * Basic memory statistics on MacOSX systems (assumes single node).
 * Diagnostics are returned as formatted output (assumed to be used
 * for model development only, since seriously impacting runtime).
 *
 * @param[in]  mpicomm  MPI communicator
 * @param[in]  str      Formatted output string
 * @param[in]  lstr     Length of str
 * @retval     0        Successful
 * @retval     1        Not successful
 **/
int ccpp_memory_usage_c(const int mpicomm, char* str, int lstr){

	// Error handling
	int ierr;
	ierr = 0;

	// Local string with correct length
	char strtmp[lstr];

	// Get MPI rank and size
	int mpirank;
	int mpisize;
	int mpiroot;
#ifdef MPI
	ierr = MPI_Comm_rank(mpicomm, &mpirank);
	if (ierr!=MPI_SUCCESS) return 1;
	ierr = MPI_Comm_size(mpicomm, &mpisize);
	if (ierr!=MPI_SUCCESS) return 1;
#else
	mpirank = 0;
	mpisize = 1;
#endif

	// Retrieve memory statistics

	// Available physical memory - in Byte
	unsigned long long totalPhysMem;
	size_t len = sizeof(totalPhysMem);
	sysctlbyname("hw.memsize", &totalPhysMem, &len, NULL, 0);
	totalPhysMem = (double)totalPhysMem/1024.0/1024.0; // convert to MB

	// Maximum physical memory (maximum resident size) per process - in Byte
	struct rusage my_rusage;
	int physMemPerProcessMax;
	my_rusage.ru_maxrss = 0;
	getrusage(RUSAGE_SELF,&my_rusage);
	physMemPerProcessMax = (double)my_rusage.ru_maxrss/1024.0/1024.0; // convert to MB

	// Physical and virtual memory per process - in Byte
	struct task_basic_info t_info;
	mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
	int virtMemPerProcess;
	int physMemPerProcess;
	task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count);
	virtMemPerProcess = (double)t_info.virtual_size/1024.0/1024.0; // convert to MB
	physMemPerProcess = (double)t_info.resident_size/1024.0/1024.0; // convert to MB

#ifdef MPI
	// Average across MPI ranks
	int virtMemPerProcessAvg;
	int physMemPerProcessAvg;
	int physMemPerProcessMaxAvg;
	ierr = MPI_Allreduce(&virtMemPerProcess, &virtMemPerProcessAvg, 1, MPI_INT, MPI_SUM, mpicomm);
	if (ierr!=MPI_SUCCESS) return 1;
	virtMemPerProcessAvg = (int)((float)virtMemPerProcessAvg/(float)mpisize);
	ierr = MPI_Allreduce(&physMemPerProcess, &physMemPerProcessAvg, 1, MPI_INT, MPI_SUM, mpicomm);
	if (ierr!=MPI_SUCCESS) return 1;
	physMemPerProcessAvg = (int)((float)physMemPerProcessAvg/(float)mpisize);
	ierr = MPI_Allreduce(&physMemPerProcessMax, &physMemPerProcessMaxAvg, 1, MPI_INT, MPI_SUM, mpicomm);
	if (ierr!=MPI_SUCCESS) return 1;
	physMemPerProcessMaxAvg = (int)((float)physMemPerProcessMaxAvg/(float)mpisize);
#endif

	// Create formatted output
#ifndef MPI
	snprintf(strtmp, sizeof(strtmp), 
"--------------------------------------------------------------------------------\n\
Memory usage - MPI rank %d/%d\n\
\n\
Size of physical memory: %llu MB\n\
\n\
Virtual    memory this process: %d MB\n\
Physical   memory this process: %d MB\n\
Max. phys. memory this process: %d MB\n\
--------------------------------------------------------------------------------\n\
", mpirank, mpisize, totalPhysMem, virtMemPerProcess, physMemPerProcess, physMemPerProcessMax);
#else
	snprintf(strtmp, sizeof(strtmp), 
"--------------------------------------------------------------------------------\n\
Memory usage - MPI rank %d/%d\n\
\n\
Size of physical memory: %llu MB\n\
\n\
Virtual    memory this process: %d MB\n\
Physical   memory this process: %d MB\n\
Max. phys. memory this process: %d MB\n\
\n\
Statistics across all processes:\n\
Virtual    memory average: %d MB\n\
Physical   memory average: %d MB\n\
Max. phys. memory average: %d MB\n\
--------------------------------------------------------------------------------\n\
", mpirank, mpisize, totalPhysMem, virtMemPerProcess, physMemPerProcess, physMemPerProcessMax,
virtMemPerProcessAvg, physMemPerProcessAvg, physMemPerProcessMaxAvg);
#endif

	// Copy local string to output string (same size!)
	strcpy(str,strtmp);

	return ierr;
}

#endif

/**
 * @}
 **/
