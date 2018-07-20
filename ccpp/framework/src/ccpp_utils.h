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
 * \file ccpp_utils.h
 *
 * CCPP utility functions.
 *
 * \ingroup CCPP
 * \{
 **/
#ifndef CCPP_UTILS_H
#define CCPP_UTILS_H

#ifdef __cplusplus
extern "C"
{
#endif

/** Resolves the absolute path when given a relative path. **/
int ccpp_abs_path(const char *, char **);

#ifdef __cplusplus
}                               /* extern "C" */
#endif

#endif                          /* CCPP_UTILS_H */

/**
 * \}
 **/
