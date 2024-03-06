#include "ray.h"

#include "cglm/call/vec3.h"
#include "cglm/cglm.h"
#include "cglm/types.h"

#include "cglm/vec3.h"
#include "cglm/quat.h"

#include <math.h>

void swap(float *t1, float *t2) {
	float tmp = *t1;
	*t1 = *t2;
	*t2 = tmp;
}

bool intersect_box(const ray *r, const vec3 *szs, const vec3 *pos, const vec4 *rot, const RGB *col, RGB *dest_col, float *dist) {
	vec4 q_;
	glm_quat_inv(*rot, q_);

	vec3 adj_r_o_tmp;
	vec3 adj_r_o;
	glm_vec3_sub(r->o, pos, adj_r_o_tmp);
	glm_quat_rotatev(q_, adj_r_o_tmp, adj_r_o);

	vec3 adj_r_d;
	glm_quat_rotatev(q_, r->d, adj_r_d);


	float t1x, t2x, t1y, t2y, t1z, t2z;	

	t1x = ( (*szs)[0] - adj_r_o[0]) / adj_r_d[0];
	t2x = (-(*szs)[0] - adj_r_o[0]) / adj_r_d[0];

	t1y = ( (*szs)[1] - adj_r_o[1]) / adj_r_d[1];
	t2y = (-(*szs)[1] - adj_r_o[1]) / adj_r_d[1];

	t1z = ( (*szs)[2] - adj_r_o[2]) / adj_r_d[2];
	t2z = (-(*szs)[2] - adj_r_o[2]) / adj_r_d[2];

	if (t1x > t2x) {
		swap(&t1x, &t2x);
	}

	if (t1y > t2y) {
		swap(&t1y, &t2y);
	}

	if (t1z > t2z) {
		swap(&t1z, &t2z);
	}

	float t1 = fmaxf(t1x, fmaxf(t1y, t1z));
	float t2 = fmaxf(t2x, fmaxf(t2y, t2z));

	if (t1 > t2) {
		return false;
	}

	if (t1 > 0) {
		// t1 -- nearest solution
		*dist = t1;
		*dest_col = *col;
		return true;
	}

	if (t1 < 0 && t2 > 0) {
		// inside sphere, t2 -- front solution
		*dist = t2;
		*dest_col = *col;
		return true;
	}

	if (t2 < 0) {
		// Sphere is behind, TODO delete this branch?
		return false;
	}

	return false;
}

bool intersect_plane(const ray *r, const vec3 *nrm, const vec3 *pos, const vec4 *rot, const RGB *col, RGB *dest_col, float *dist) {
	vec4 q_;
	glm_quat_inv(*rot, q_);

	vec3 adj_r_o_tmp;
	vec3 adj_r_o;
	glm_vec3_sub(r->o, pos, adj_r_o_tmp);
	glm_quat_rotatev(q_, adj_r_o_tmp, adj_r_o);

	vec3 adj_r_d;
	glm_quat_rotatev(q_, r->d, adj_r_d);

	float t = -glm_vec3_dot(adj_r_o, *nrm) / glm_vec3_dot(adj_r_d, *nrm);
	
	if (t >= 0) {
		*dist = t;
		*dest_col = *col;
		return true;
	}

	return false;
}

bool intersect_ellipsoid(const ray *r, const vec3 *rads, const vec3 *pos, const vec4 *rot, const RGB *col, RGB *dest_col, float *dist) {
	vec4 q_;
	glm_quat_inv(*rot, q_);

	vec3 adj_r_o_tmp;
	vec3 adj_r_o;
	glm_vec3_sub(r->o, *pos, adj_r_o_tmp);
	glm_quat_rotatev(q_, adj_r_o_tmp, adj_r_o);

	vec3 adj_r_d;
	glm_quat_rotatev(q_, r->d, adj_r_d);

	vec3 dr, or;
	glm_vec3_div(adj_r_o, *rads, or);
	glm_vec3_div(adj_r_d, *rads, dr);

	float a, b, c;
	a = glm_vec3_dot(dr, dr);
	b = 2 * glm_vec3_dot(or, dr);
	c = glm_vec3_dot(or, or) - 1;

	float D = b * b - 4 * a * c;
	if (D < 0) {
		return false;
	}
	float Dsqrt = sqrtf(D);

	float t1 = (Dsqrt + b) / (2 * a);
	float t2 = (Dsqrt - b) / (2 * a);
	if (t2 > t1) {
		swap(&t1, &t2);
	}
	
	if (t1 > 0) {
		// t1 -- nearest solution
		*dist = t1;
		*dest_col = *col;
		return true;
	}

	if (t1 < 0 && t2 > 0) {
		// inside sphere, t2 -- front solution
		*dist = t2;
		*dest_col = *col;
		return true;
	}

	if (t2 < 0) {
		// Sphere is behind, TODO delete this branch?
		return false;
	}

	return false;
}

RGB process_ray(const scene *s, const ray *r) {
	size_t b_cnt = s->bxs.n;
	size_t p_cnt = s->plns.n;
	size_t e_cnt = s->elps.n;

	bool intersected = false;
	float t;
	RGB col = s->bg;

	for (size_t i = 0; i < b_cnt; i++) {
		float tmp_t;
		RGB tmp_c;
		bool succ = intersect_box(
				r, 
				&s->bxs.szs[i], 
				&s->bxs.comm_data.pos[i],
				&s->bxs.comm_data.rot[i], 
				&s->bxs.comm_data.col[i], 
				&tmp_c, 
				&tmp_t);
		
		if (succ) {
			if ((intersected && tmp_t < t) || !intersected) {
				t = tmp_t;
				col = tmp_c;
			}

			intersected = true;
		}
	}

	for (size_t i = 0; i < p_cnt; i++) {
		float tmp_t;
		RGB tmp_c;
		bool succ = intersect_plane(
				r, 
				&s->plns.nrms[i], 
				&s->plns.comm_data.pos[i],
				&s->plns.comm_data.rot[i], 
				&s->plns.comm_data.col[i], 
				&tmp_c, 
				&tmp_t);
		
		if (succ) {
			if ((intersected && tmp_t < t) || !intersected) {
				t = tmp_t;
				col = tmp_c;
			}

			intersected = true;
		}
	}

	for (size_t i = 0; i < e_cnt; i++) {
		float tmp_t;
		RGB tmp_c;
		bool succ = intersect_ellipsoid(
				r, 
				&s->elps.rads[i], 
				&s->elps.comm_data.pos[i],
				&s->elps.comm_data.rot[i], 
				&s->elps.comm_data.col[i], 
				&tmp_c, 
				&tmp_t);
		
		if (succ) {
			if ((intersected && tmp_t < t) || !intersected) {
				t = tmp_t;
				col = tmp_c;
			}

			intersected = true;
		}
	}

	return col;
}

RGB *render(const scene *s) {
	RGB *render = malloc(sizeof(RGB) * s->cam.h * s->cam.w);

	float aspect_ratio = (float) s->cam.w / s->cam.h;
  float fov_y = 2 * atan(tanf(s->cam.fov_x / 2) / aspect_ratio);

	for (size_t x = 0; x < s->cam.w; x++) {
		for (size_t y = 0; y < s->cam.h; y++) {
			ray r;
			glm_vec3_copy(s->cam.pos, r.o);

			float px = (2 * (float) x / s->cam.w - 1) * tanf(s->cam.fov_x / 2);
			float py = (2 * (float) y / s->cam.h - 1) * tanf(       fov_y / 2);

      vec3 tmp;
      glm_vec3_copy(s->cam.frwrd, r.d);

			glm_vec3_scale(s->cam.right, px, tmp);
      glm_vec3_add(tmp, r.d, r.d);

			glm_vec3_scale(s->cam.up, py, tmp);
      glm_vec3_add(tmp, r.d, r.d);

      glm_vec3_normalize(r.d);

			render[x + y * s->cam.w] = process_ray(s, &r);
		}
	}

	return render;
}

