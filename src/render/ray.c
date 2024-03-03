#include "ray.h"

#include "cglm/cglm.h"

bool intersect_box(const ray *r, const vec3 *szs, const vec3 *pos, const vec3 *rot, const RGB *col, RGB *dest_col) {

	return false;
}

bool intersect_plane(const ray *r, const vec3 *szs, RGB *dest_col);
bool intersect_ellipsoid(const ray *r, const vec3 *szs, RGB *dest_col);
