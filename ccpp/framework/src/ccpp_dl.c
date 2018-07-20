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
 * @file ccpp_dl.c
 *
 * Routines for the function/subroutine calls using dynamic loaded shared
 * objects.
 *
 * @ingroup CCPP
 * @{
 **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <dlfcn.h>
#include <err.h>
#include <sysexits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "ccpp_dl.h"

/** Shared library prefix and suffix for different platforms **/
static const char prefix[] = "lib";
#if __APPLE__
static const char suffix[] = ".dylib";
#elif __unix__
static const char suffix[] = ".so";
#endif

/**
 * Function call initialization routine.
 *
 * This dlopen()'s the library specified and tries to
 * obtain a handle to the function/scheme cap.
 *
 * @param[in]  scheme    The scheme name to call.
 * @param[in]  lib       The library continaing the physics scheme.
 * @param[in]  ver       The library version number.
 * @param[out] shdl      The scheme function pointer handle.
 * @param[out] lhdl      The library handle.
 * @retval     0         If it was sucessful
 * @retval     1         If there was an error
 **/
int
ccpp_dl_open(const char *scheme, const char *lib, const char *ver,
	     void **shdl, void **lhdl)
{
	int i = 0;
	int n = 0;
	const char cap[] = "_cap";
	const char *l = NULL;
	char *library = NULL;
	char *scheme_cap = NULL;
	char *error = NULL;
	struct stat sbuf = {0};

	/* Did we get an actual library file? */
	if (stat(lib, &sbuf) == 0) {
		l = lib;
	} else {
		/* Generate the library name with the platform suffix */
		n = (strlen(prefix) + strlen(lib) + strlen(suffix)
		     + strlen(ver) +2) *sizeof(char);
		library = malloc(n);
		memset(library, 0, n);
		if (strcmp(ver, "") != 0) {
#ifdef __APPLE__
			snprintf(library, n, "%s%s.%s%s", prefix, lib,
				 ver, suffix);
#elif defined(__linux__) || defined(__unix__)
			snprintf(library, n, "%s%s%s.%s", prefix, lib,
				 suffix, ver);
#else
		 	warnx("CCPP library name not configured for this operating system");
		 	return(EXIT_FAILURE);
#endif
		} else {
			snprintf(library, n, "%s%s%s", prefix, lib, suffix);
		}
		l = library;
	}

	/* Generate the scheme cap function name */
	n = (strlen(scheme) +strlen(cap) +1)*sizeof(char);
	scheme_cap = malloc(n);
	memset(scheme_cap, 0, n);

	n = strlen(scheme);
	for (i=0; i < n; ++i) {
		scheme_cap[i] = tolower(scheme[i]);
	}

	strncat(scheme_cap, cap, n);

	/* Open a handle to the library */
	*lhdl = dlopen(l, RTLD_NOW);
	if (!*lhdl) {
		warnx("%s", dlerror());
		return(EXIT_FAILURE);
	}

	dlerror();
	*(void **)shdl = dlsym(*lhdl, scheme_cap);
	if ((error = dlerror()) != NULL)  {
		warnx("%s", error);
		return(EXIT_FAILURE);
	}

	/* Free the library filename */
	if (library) {
		free(library);
		library = NULL;
	}

	/* Free the scheme cap function name */
	if (scheme_cap) {
		free(scheme_cap);
		scheme_cap = NULL;
	}

	return(EXIT_SUCCESS);
}

/**
 * Function call library closing routine.
 *
 * @param[in] lhdl      The library handle.
 * @retval     0        If it was sucessful
 * @retval     1        If there was an error
 **/
int
ccpp_dl_close(void **lhdl)
{
	char *error = NULL;

	dlerror();
	dlclose(*lhdl);
	if ((error = dlerror()) != NULL)  {
		warnx("%s", error);
		return(EXIT_FAILURE);
	}

	return(EXIT_SUCCESS);
}

/**
 * The function cap calling routine.
 *
 * @param[in] f_ptr     The scheme function pointer to call.
 * @param[in] data      The opaque ccpp_t data type to pass.
 * @retval     0        If it was sucessful
 * @retval    !=0       If there was an error
 **/
int
ccpp_dl_call(void **f_ptr, void **data)
{
	int (*fun)(void **) = NULL;

	*(int **)(&fun) = *f_ptr;

	return(fun(data));
}

/**
 * @}
 **/
