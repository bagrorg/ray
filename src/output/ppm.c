#include "ppm.h"

#include <stdio.h>
#include <errno.h>
#include <string.h>

#include "common/assert.h"

void write_ppm(const char *path, size_t w, size_t h, const RGB *data) {
	FILE *f = fopen(path, "w");

	RAY_VERIFY(f != NULL, "Failed to open %s for output (%s)", path, strerror(errno));

	fprintf(f, "P6\n");
	fprintf(f, "%ld %ld\n", w, h);
	fprintf(f, "255\n");
	fwrite(data, sizeof(*data), w * h, f);

	fclose(f);
}


