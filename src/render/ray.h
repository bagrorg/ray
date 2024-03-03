#ifndef RAY_H
#define RAY_H

#include <stdbool.h>

#include "common/math.h"
#include "common/scene.h"

typedef struct {
	vec3 o;			// origin
	vec3 d;			// direction
} ray;

bool intersect_box(const ray *r, const vec3 *szs, const vec3 *pos, const vec3 *rot, const RGB *col, RGB *dest_col);
bool intersect_plane(const ray *r, const vec3 *szs, RGB *dest_col);
bool intersect_ellipsoid(const ray *r, const vec3 *szs, RGB *dest_col);


#endif
