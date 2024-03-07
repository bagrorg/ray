#ifndef PPM_H
#define PPM_H

#include <stdint.h>
#include <stddef.h>

#include "common/scene.h"

void write_ppm(const char *path, size_t w, size_t h, const RGB *data);

#endif
