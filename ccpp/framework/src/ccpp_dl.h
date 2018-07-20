/**
 * @file ccpp_dl.h
 *
 * The function pointer routines using dynamic loaded shared objects.
 *
 * @ingroup CCPP
 * @{
 **/
#ifndef CCPP_DL_H
#define CCPP_DL_H

#ifdef __cplusplus
extern "C"
{
#endif

/** Function libaray and cap function initialization routine. **/
int ccpp_dl_open(const char *, const char *, const char *, void **, void **);

/** Function library closing/unloading routine. **/
int ccpp_dl_close(void **);

/** Function pointer physics cap function call. **/
int ccpp_dl_call(void **, void **);

#ifdef __cplusplus
}                               /* extern "C" */
#endif

#endif                          /* CCPP_DL_H */

/**
 * @}
 **/
