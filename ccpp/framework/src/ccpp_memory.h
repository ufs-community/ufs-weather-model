/**
 * @file ccpp_memory.h
 *
 * Memory usage routines.
 *
 * @ingroup CCPP
 * @{
 **/
#ifndef CCPP_MEMORY_H
#define CCPP_MEMORY_H

#ifdef __cplusplus
extern "C"
{
#endif

#if __unix__

/** Parses lines in /proc/self/status **/
static int parseLine(char*);

/** Virtual memory used per process **/
static int getVirtMemPerProcess(void);

/** Physical memory used per process **/
static int getPhysMemPerProcess(void);

/** Maximum physical memory used per process **/
static int getPhysMemPerProcessMax(void);

/** Memory usage statistics **/
int ccpp_memory_usage_c(const int, char*, int);

#elif __APPLE__

/** Memory usage statistics **/
int ccpp_memory_usage_c(const int, char*, int);

#endif

#ifdef __cplusplus
}                               /* extern "C" */
#endif

#endif                          /* CCPP_MEMORY_H */

/**
 * @}
 **/
