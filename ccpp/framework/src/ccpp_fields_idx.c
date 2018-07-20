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
 * @file ccpp_fields_idx.c
 *
 * @brief Routines and functions to generate and lookup fields/variables
 *        needed for the physics routines.
 *
 * @details The fields are stored in an array of C pointers within the
 *          ccpp_t type. There is also an index array in this type.
 *          We poppulate this index array with the standard name of
 *          each variable in the fields array. We use a binary search
 *          on the sorted index array to retreive the array index for
 *          the field witin the fields array.
 *
 * TODO 
 * - Test the sort and lookup times for qsort() and bsearch().
 * - Implement this as a hash-map instead.
 *
 * @ingroup Physics
 * @{
 **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <err.h>
#include <sysexits.h>
#include <assert.h>

#include "ccpp_fields_idx.h"

/**
 * Comparison function.
 *
 * Compares the name of two index elements using strcmp().
 * It returns an integer less than, equal to, or greater than
 * zero if the name in f1 is  found, respectively, to be less
 * than, to match, or be greater than the name in f2.
 *
 * @param[in] f1   The first field.
 * @param[in] f2   The second field.
 **/
static int
cmp(const void *f1, const void *f2)
{
	struct ccpp_field *f_1;
	struct ccpp_field *f_2;
	f_1 = *(struct ccpp_field * const *) f1;
	f_2 = *(struct ccpp_field * const *) f2;
	return strcmp(f_1->name, f_2->name);
}

/**
 * Initialization routine.
 *
 * Allocates an array for the field indices.
 *
 * @param[in,out] index The index array.
 * @retval        0     If it was sucessful.
 * @retval        1     If there was an error.
 **/
int
ccpp_field_idx_init(void **index)
{
	struct ccpp_field_idx *f_index;
	f_index = NULL;

	*index = (struct ccpp_field_idx *)malloc(sizeof(struct ccpp_field_idx));
	if (*index == NULL) {
		warnx("Unable to allocate field index");
		return(EXIT_FAILURE);
	}

	f_index = (struct ccpp_field_idx *)(*index);

	f_index->sorted = 0;
	f_index->n      = 0;
	f_index->max    = CCPP_FIELD_IDX_MAX;
	f_index->fields = malloc(CCPP_FIELD_IDX_MAX * sizeof(struct ccpp_field *));

	return(EXIT_SUCCESS);
}

/**
 * Finalization routine.
 *
 * Deallocates the field indices array.
 *
 * @param[in] index The index array.
 * @retval    0     If it was sucessful.
 **/
int
ccpp_field_idx_finalize(void **index)
{
	int i;

	struct ccpp_field_idx *f_index;

	f_index = (struct ccpp_field_idx *)(*index);

	for (i = 0; i < f_index->n; ++i) {
		if (f_index->fields[i]->name) {
			free(f_index->fields[i]->name);
			f_index->fields[i]->name = NULL;
		}
		free(f_index->fields[i]);
		f_index->fields[i] = NULL;
	}
	free(f_index->fields);
	f_index->fields = NULL;

	free(f_index);
	f_index = NULL;

	return(EXIT_SUCCESS);
}

/**
 * Add/Insert a field into the index.
 *
 * @param[in]     name  The name to add to the index array.
 * @param[in,out] index The index array.
 * @retval        > 0   The index location.
 * @retval        -1    If there was an error.
 **/
int
ccpp_field_idx_add(const char *name, void **index)
{
	struct ccpp_field_idx *f_index;
	int n;
	size_t len;
	f_index = (struct ccpp_field_idx *)(*index);
	n = 0;
	len = 0;

	n = f_index->n;
	if (n == f_index->max) {
		if (ccpp_field_idx_grow(index)) {
			warnx("Unable to grow field index array");
			return(-1);
		}
	}
	f_index->fields[n] = malloc(sizeof(struct ccpp_field));

	len = strlen(name);

	f_index->fields[n]->name = malloc((len + 1) * sizeof(char));

	strncpy(f_index->fields[n]->name, name, len * sizeof(char));
	f_index->fields[n]->name[len] = '\0';
	f_index->fields[n]->n = n+1;
	f_index->sorted = 0;
	f_index->n++;

	return(n+1);
}


/**
 * Find the index number of a field.
 *
 * @param[in]     name     The field name to find the index array.
 * @param[in,out] index    The index array.
 * @retval        > 0      The position in the index array of the requested field.
 * @retval        -1       If there was an error.
 **/
int
ccpp_field_idx_find(const char *name, void **index)
{
	int n;
	struct ccpp_field  *key;
	struct ccpp_field **res;
	struct ccpp_field_idx *f_index;
	n = 0;
	key = NULL;
	res = NULL;
	f_index = (struct ccpp_field_idx *)(*index);

	if (f_index->sorted == 0) {
		ccpp_field_idx_sort(index);
	}

	key = malloc(sizeof(struct ccpp_field));
	n = strlen(name);
	key->name = malloc((n+1) * sizeof(char));
	strncpy(key->name, name, n);
	key->name[n] = '\0';

	res = bsearch(&key, f_index->fields, f_index->n,
		          sizeof(struct ccpp_field *), cmp);
	if (*res == NULL) {
		warnx("Unable to find in index: %s", name);
		return(-1);
	}

	free(key->name);
	free(key);

	return((*res)->n);
}

/**
 * Sort the index by calling qsort() and using cmp() as the
 * comparison function.
 *
 * @param[in,out] index    The index array.
 * @retval        0        If there was no error.
 **/
static int
ccpp_field_idx_sort(void **index)
{
	struct ccpp_field_idx *f_index;
	f_index = (struct ccpp_field_idx *)(*index);

	qsort(f_index->fields, f_index->n, sizeof(struct ccpp_field *), cmp);
	f_index->sorted = 1;

	return(EXIT_SUCCESS);
}

/**
 * Grow the index field array.
 *
 * @param[in,out] index    The index array.
 * @retval        0        If there was no error.
 **/
static int
ccpp_field_idx_grow(void **index)
{
	// Warn user that field index array needs to grow
	warnx("Growing field index array");

	struct ccpp_field_idx *f_index;
	struct ccpp_field **new;
	int new_max;
	f_index = (struct ccpp_field_idx *)(*index);
	new = NULL;
	new_max = 0;

	new_max = f_index->max * CCPP_FIELD_IDX_GROW;

	new = realloc(f_index->fields, new_max * sizeof(struct ccpp_field *));
	if (new == NULL) {
		warnx("Unable to expand the field index array");
		return(EXIT_FAILURE);
	}
	f_index->fields = new;
	f_index->max = new_max;

	return(EXIT_SUCCESS);
}

/**
 * Get the maximum number of fields the index array can hold.
 *
 * @param[in,out] index    The index array.
 * @retval        >= 0     The maximum number of fields.
 **/
int
ccpp_field_idx_max(void **index)
{
	struct ccpp_field_idx *f_index;
	f_index = (struct ccpp_field_idx *)(*index);

	assert(f_index->max > 0);

	return(f_index->max);

}

/**
 * @}
 **/
