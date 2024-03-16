#include "ray.h"

#include "cglm/types.h"

#include "cglm/vec3.h"
#include "cglm/quat.h"
#include "common/assert.h"
#include "common/scene.h"

#include <math.h>

void swap(float *t1, float *t2) {
	float tmp = *t1;
	*t1 = *t2;
	*t2 = tmp;
}

void prepare_ray(const ray *r, const vec3 *pos, const vec4 *rot, ray *r_new) {
	vec4 q_;
	glm_quat_inv(*rot, q_);

	glm_vec3_sub(r->o, *pos, r_new->o);
	glm_quat_rotatev(q_, r_new->o, r_new->o);     // POTENTIALLY DANGEROUS

	glm_quat_rotatev(q_, r->d, r_new->d);
}

void postprocess_normal(vec3 *N, const vec4 *rot) {
  glm_quat_rotatev(*rot, *N, *N); 
}

void compute_int_point(const ray *r, float t, vec3 *P) {
    vec3 step;
    glm_vec3_copy(r->o, *P);
    glm_vec3_copy(r->d, step);
    glm_vec3_scale(step, t, step);
    glm_vec3_add(*P, step, *P);
}

intsec intersect_box(const ray *r, const primitive *p) {
  intsec i = {
    .succ = false
  };
  ray adj_r;
  prepare_ray(r, &p->position, &p->rotation, &adj_r);

  vec3 ts1, ts2;
  glm_vec3_copy(p->sizes, ts1);
  glm_vec3_copy(p->sizes, ts2);
  glm_vec3_scale(ts1, -1, ts1);

  glm_vec3_sub(ts1, adj_r.o, ts1);
  glm_vec3_sub(ts2, adj_r.o, ts2);

  glm_vec3_div(ts1, adj_r.d, ts1);
  glm_vec3_div(ts2, adj_r.d, ts2);

  for (size_t i = 0; i < 3; i++) {
    if (ts1[i] > ts2[i]) {
      float tmp = ts1[i];
      ts1[i] = ts2[i];
      ts2[i] = tmp;
    } 
  }

  float t1 = fmaxf(ts1[0], fmaxf(ts1[1], ts1[2]));
  float t2 = fminf(ts2[0], fminf(ts2[1], ts2[2]));

	if (t1 > t2) {
    return i;
	}

	if (t1 >= 0) {
		// t1 -- nearest solution
    i.t = t1;
    glm_vec3_copy(p->colour, i.col);
    i.succ = true;
	} else if (t2 >= 0) {
		// t1 -- nearest solution
    i.t = t2;
    glm_vec3_copy(p->colour, i.col);
    i.succ = true;
	}

  if (i.succ) {
    vec3 P;
    compute_int_point(&adj_r, i.t, &P);
    glm_vec3_div(P, p->sizes, P);
    for (size_t i = 0; i < 3; i++) {
      if (fabs(P[i] - 1.0) > 0.00001) {
        P[i] = 0.0;
      }
    }
    glm_vec3_copy(P, i.N);
    glm_vec3_normalize(i.N);
    
    postprocess_normal(&i.N, &p->rotation);
  }

  return i;
}

intsec intersect_plane(const ray *r, const primitive *p) {
  intsec i = {
    .succ = false,
  };
  ray adj_r;
  prepare_ray(r, &p->position, &p->rotation, &adj_r);

	float t = -glm_vec3_dot(adj_r.o, p->normal) / glm_vec3_dot(adj_r.d, p->normal);
	
	if (t >= 0) {
    i.succ = true;
    i.t = t;
    glm_vec3_copy(p->colour, i.col);
    glm_vec3_copy(p->normal, i.N);
    glm_vec3_normalize(i.N);
    postprocess_normal(&i.N, &p->rotation);
	}

  return i;
}

intsec intersect_ellipsoid(const ray *r, const primitive *p) {
  intsec i = {
    .succ = false,
  };
  ray adj_r;
  prepare_ray(r, &p->position, &p->rotation, &adj_r);

	vec3 dr, or;
	glm_vec3_div(adj_r.o, p->radius, or);
	glm_vec3_div(adj_r.d, p->radius, dr);

	float a, b, c;
	a = glm_vec3_dot(dr, dr);
	b = 2 * glm_vec3_dot(or, dr);
	c = glm_vec3_dot(or, or) - 1;

	float D = b * b - 4 * a * c;
	if (D < 0) {
		return i;
	}
	float Dsqrt = sqrtf(D);

	float t1 = (Dsqrt - b) / (2 * a);
	float t2 = (-Dsqrt - b) / (2 * a);
	if (t2 < t1) {
		swap(&t1, &t2);
	}
  
	if (t1 > 0) {
		// t1 -- nearest solution
    i.t = t1;
    i.succ = true;
    i.l = OUTER;
    glm_vec3_copy(p->colour, i.col);
	} else if (t1 < 0 && t2 > 0) {
		// inside sphere, t2 -- front solution
    i.t = t2;
    i.succ = true;
    i.l = INNER;
    glm_vec3_copy(p->colour, i.col);
	} else if (t2 < 0) {
		// Sphere is behind, TODO delete this branch?
    return i;
	}

  if (i.succ) {
    vec3 P, R2;
    compute_int_point(&adj_r, i.t, &P);
    glm_vec3_mul(p->radius, p->radius, R2); 
    glm_vec3_div(P, R2, i.N);
    glm_vec3_normalize(i.N);
    postprocess_normal(&i.N, &p->rotation);
  }

  return i;
}

void attenuation(vec3 I, vec3 C, float R) {
  RAY_VERIFY(C[0] > 0, "C[0] == 0");

  float w = C[0] + C[1] * R + C[2] * R * R;
  glm_vec3_divs(I, w, I);
}

void process_light(const scene *s, const ray *r, const intsec *i, RGB *dest) {
  RGB col;
  vec3 sum = {(float) s->ambient[0] / 255, 
              (float) s->ambient[1] / 255, 
              (float) s->ambient[2] / 255};
  for (size_t l = 0; l < s->lights_directed.size; l++) {
    light_directed *ld = &((light_directed*)s->lights_directed.data)[l];

    vec3 P;
    glm_vec3_scale(r->d, i->t, P);
    glm_vec3_add(r->o, P, P);
    
    ray light_ray;
    glm_vec3_copy(P, light_ray.o);
    glm_vec3_copy(ld->direction, light_ray.d);

    // shift to prevent shadow acne
    vec3 shift;
    glm_vec3_copy(i->N, shift);
    glm_vec3_scale(shift, 1e-5, shift);
    glm_vec3_add(shift, light_ray.o, light_ray.o);
    
    intsec i_temp;
    process_ray(&i_temp, s, &light_ray, INFINITY);    

    if (i_temp.succ) {
      continue;
    }
    
    float light_dot = glm_vec3_dot(i->N, ld->direction);
    if (light_dot < 0) {
      continue;
    }
    
    vec3 w;
    glm_vec3_scale(ld->intensity, light_dot, w);
    glm_vec3_add(sum, w, sum);
  }

  for (size_t l = 0; l < s->lights_point.size; l++) {
    light_point *lp = &((light_point*)s->lights_point.data)[l];

    vec3 P;
    glm_vec3_scale(r->d, i->t, P);
    glm_vec3_add(r->o, P, P);

    vec3 light_dir;
    float R;
    glm_vec3_sub(lp->position, P, light_dir);
    R = glm_vec3_norm(light_dir);
    glm_vec3_normalize(light_dir);

    ray light_ray;
    glm_vec3_copy(P, light_ray.o);
    glm_vec3_copy(light_dir, light_ray.d);

    // shift to prevent shadow acne
    vec3 shift;
    glm_vec3_copy(i->N, shift);
    glm_vec3_scale(shift, 1e-5, shift);
    glm_vec3_add(shift, light_ray.o, light_ray.o);

    intsec i_temp;
    process_ray(&i_temp, s, &light_ray, R);
    if (i_temp.succ) {
      continue;
    }

    float light_dot = glm_vec3_dot(i->N, light_dir);
    if (light_dot < 0) {
      continue;
    }
    
    vec3 w;
    glm_vec3_scale(lp->intensity, light_dot, w);
    attenuation(w, lp->attenuation, R);

    glm_vec3_add(sum, w, sum);
  }

  glm_vec3_mul(i->col, sum, *dest);
}

void process_ray(intsec *res_glob, const scene *s, const ray *r, float max_depth) {
  res_glob->succ = false;

	for (size_t i = 0; i < s->primitives.size; i++) {
    primitive *p = &((primitive*) s->primitives.data)[i];
    intsec res;

    switch (p->type) {
      case PRIMITIVE_BOX:
        res = intersect_box(r, p);  
        break;
      case PRIMITIVE_PLANE:
        res = intersect_plane(r, p);
        break;
      case PRIMITIVE_ELLIPSOID:
        res = intersect_ellipsoid(r, p);
        break;
    };

		if (res.succ && res.t <= max_depth) {
			if ((res_glob->succ && res.t < res_glob->t) || !res_glob->succ) {
        *res_glob = res;
			}
		}
	}  
}

void process_ray_regular(intsec *res_glob, const scene *s, const ray *r) {
  process_ray(res_glob, s, r, INFINITY);
}

RGB *render(const scene *s) {
	RGB *render = malloc(sizeof(RGB) * s->cam.h * s->cam.w);

	float aspect_ratio = (float) s->cam.w / s->cam.h;
  float fov_y = 2 * atan(tanf(s->cam.fov_x / 2) / aspect_ratio);

	for (size_t x = 0; x < s->cam.w; x++) {
		for (size_t y = 0; y < s->cam.h; y++) {
			ray r;
			glm_vec3_copy(s->cam.pos, r.o);

			float px =  (2 * ((float) x + 0.5) / s->cam.w - 1) * tanf(s->cam.fov_x / 2);
			float py = -(2 * ((float) y + 0.5) / s->cam.h - 1) * tanf(       fov_y / 2);

      vec3 tmp;
      glm_vec3_copy(s->cam.frwrd, r.d);

			glm_vec3_scale(s->cam.right, px, tmp);
      glm_vec3_add(tmp, r.d, r.d);

			glm_vec3_scale(s->cam.up, py, tmp);
      glm_vec3_add(tmp, r.d, r.d);

      glm_vec3_normalize(r.d);

      intsec res;
			process_ray_regular(&res, s, &r);

      if (res.succ) {
        RGB col;
        process_light(s, &r, &res, &col);
        glm_vec3_copy(col, render[x + y * s->cam.w]);
      } else {
        glm_vec3_copy(s->bg, render[x + y * s->cam.w]);
      }
		}
	}

	return render;
}

