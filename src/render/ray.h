#ifndef RAY_H
#define RAY_H

#include <stdbool.h>
#include "common/scene.h"

typedef struct {
	vec3 o;			// origin
	vec3 d;			// direction
} ray;

bool intersect_box(const ray *r, const vec3 *szs, const vec3 *pos, const vec4 *rot, const RGB *col, RGB *dest_col, float *dist);
bool intersect_plane(const ray *r, const vec3 *nrm, const vec3 *pos, const vec4 *rot, const RGB *col, RGB *dest_col, float *dist);
bool intersect_ellipsoid(const ray *r, const vec3 *rads, const vec3 *pos, const vec4 *rot, const RGB *col, RGB *dest_col, float *dist);

RGB *render(const scene *s);

#endif
