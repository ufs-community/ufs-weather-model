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
 * A test to make sure the field array is growable.
 **/

#include <stdio.h>
#include <stdlib.h>

#include "ccpp_fields_idx.h"

int
main(int argc, char **argv)
{
	int i = 0;
	int n = 100;
	char f[10] = {0};
	void *cdata = NULL;

	if (ccpp_field_idx_init(&cdata)) {
		return(EXIT_FAILURE);
	}

	for (i = 0; i < n; ++i) {
		sprintf(f, "f_%d", i);
		if (ccpp_field_idx_add(f, &cdata) <= 0) {
			return(EXIT_FAILURE);
		}
	}

	i = ccpp_field_idx_find("f_90", &cdata);
	printf("%d\n", i);

	if (ccpp_field_idx_finalize(&cdata)) {
		return(EXIT_FAILURE);
	}

	return(EXIT_SUCCESS);
}
