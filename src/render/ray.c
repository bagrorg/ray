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
    .p = NULL
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
    i.p = p;
	} else if (t2 >= 0) {
		// t1 -- nearest solution
    i.t = t2;
    i.p = p;
	}

  if (i.p != NULL) {
    vec3 P;
    compute_int_point(&adj_r, i.t, &P);
    glm_vec3_div(P, p->sizes, P);
    glm_vec3_div(adj_r.o, p->sizes, adj_r.o);
    i.l = INNER;

    for (size_t c = 0; c < 3; c++) {
      if (fabs(fabs(P[c]) - 1.0) > 0.00001) {
        P[c] = 0.0;
      }

      if (fabs(adj_r.o[c]) > 1.f) {
        i.l = OUTER;
      }
    }
    glm_vec3_copy(P, i.N);
    glm_vec3_normalize(i.N);
    if (i.l == INNER) {
      glm_vec3_negate(i.N);
    }
    
    postprocess_normal(&i.N, &p->rotation);
  }

  return i;
}

intsec intersect_plane(const ray *r, const primitive *p) {
  intsec i = {
    .p = NULL
  };
  ray adj_r;
  prepare_ray(r, &p->position, &p->rotation, &adj_r);

	float t = -glm_vec3_dot(adj_r.o, p->normal) / glm_vec3_dot(adj_r.d, p->normal);
	
	if (t >= 0) {
    i.p = p;
    i.t = t;

    glm_vec3_copy(p->normal, i.N);
    glm_vec3_normalize(i.N);
    postprocess_normal(&i.N, &p->rotation);
	}

  return i;
}

intsec intersect_ellipsoid(const ray *r, const primitive *p) {
  intsec i = {
    .p = NULL
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
    i.p = p;
    i.l = OUTER;
	} else if (t1 < 0 && t2 > 0) {
		// inside sphere, t2 -- front solution
    i.t = t2;
    i.p = p;
    i.l = INNER;
	} else if (t2 < 0) {
		// Sphere is behind, TODO delete this branch?
    return i;
	}

  if (i.p != NULL) {
    vec3 P, R2;
    compute_int_point(&adj_r, i.t, &P);
    glm_vec3_mul(p->radius, p->radius, R2); 
    glm_vec3_div(P, R2, i.N);
    glm_vec3_normalize(i.N);
    postprocess_normal(&i.N, &p->rotation);

    if (i.l == INNER) {
      glm_vec3_negate(i.N);
    }
  }

  return i;
}

void attenuation(vec3 I, vec3 C, float R) {
  RAY_VERIFY(C[0] > 0, "C[0] == 0");

  float w = C[0] + C[1] * R + C[2] * R * R;
  glm_vec3_divs(I, w, I);
}

intsec intersect_ray(const scene *s, const ray *r, float max_depth) {
  intsec res_glob = {
    .p = NULL
  };

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

		if (res.p != NULL && res.t <= max_depth) {
			if ((res_glob.p != NULL && res.t < res_glob.t) || res_glob.p == NULL) {
        res_glob = res;
			}
		}
	}

  return res_glob;
}

void process_light(const scene *s, const ray *r, const intsec *i, RGB *dest) {
  RGB col;
  vec3 sum = {(float) s->ambient[0], 
              (float) s->ambient[1], 
              (float) s->ambient[2]};

  // Direction lights processing
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
    
    intsec i_temp = intersect_ray(s, &light_ray, INFINITY);    

    if (i_temp.p != NULL) {
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
  
  // Point lights processing
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

    intsec i_temp = intersect_ray(s, &light_ray, R);
    if (i_temp.p != NULL) {
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

  glm_vec3_mul(i->p->colour, sum, *dest);
}

void process_ray(RGB *dest, const scene *s, const ray *r, float max_depth, size_t recuesion_depth);

void process_reflection(RGB *dest, const scene *s, const ray *r, const intsec *i, float max_depth, size_t recursion_depth) {
  ray reflected_ray;
  
  // REFL.O = R->O + R->D * i->t + eps * i->N
  glm_vec3_copy(r->d, reflected_ray.o);
  glm_vec3_scale(reflected_ray.o, i->t, reflected_ray.o);
  glm_vec3_add(r->o, reflected_ray.o, reflected_ray.o);

  vec3 shift;
  glm_vec3_copy(i->N, shift);
  glm_vec3_scale(shift, 1e-4, shift);
  glm_vec3_add(shift, reflected_ray.o, reflected_ray.o);

  // REFL.D = R->D - 2 * i->N * <i->N, R->D>
  glm_vec3_copy(i->N, reflected_ray.d);
  glm_vec3_scale(reflected_ray.d, 
                 -2 * glm_vec3_dot(r->d, i->N), 
                 reflected_ray.d);
  glm_vec3_add(reflected_ray.d, r->d, reflected_ray.d);
  glm_vec3_normalize(reflected_ray.d);

  process_ray(dest, s, &reflected_ray, max_depth, recursion_depth); 
}

void process_refraction(
  RGB *dest, 
  const scene *s, 
  const ray *r, 
  const intsec *i, 
  float max_depth, 
  size_t recursion_depth, 
  float cos_th1, 
  float sin_th2, 
  float mu1, 
  float mu2) {

  float cos_th2 = sqrtf(1 - sin_th2 * sin_th2);

  // TODO reduce amount of origin computations
  ray refracted_ray;
 
  // REFL.O = R->O + R->D * i->t + eps * i->N
  glm_vec3_copy(r->d, refracted_ray.o);
  glm_vec3_scale(refracted_ray.o, i->t, refracted_ray.o);
  glm_vec3_add(r->o, refracted_ray.o, refracted_ray.o);

  vec3 shift;
  glm_vec3_copy(i->N, shift);
  glm_vec3_scale(shift, -1e-4, shift);
  glm_vec3_add(shift, refracted_ray.o, refracted_ray.o);

  // REFL.D = left + right = mu1/mu2 * (R->D) + (mu1/mu2 cos_th1 - cos_th2) * i->N
  vec3 left, right;
  float mu1_mu2 = mu1 / mu2;
  glm_vec3_scale(r->d, mu1_mu2, left);
  glm_vec3_scale(i->N, mu1_mu2 * cos_th1 - cos_th2, right);
  glm_vec3_add(left, right, refracted_ray.d);
  glm_vec3_normalize(refracted_ray.d);
  
  process_ray(dest, s, &refracted_ray, max_depth, recursion_depth); 
}

void process_dielectric(RGB *dest, const scene *s, const ray *r, const intsec *i, float max_depth, size_t recursion_depth) {
  float mu1 = 1.f;
  float mu2 = i->p->ior;

  if (i->l == INNER) {
    swap(&mu1, &mu2);
  }

  float cos_th1 = -1 * glm_vec3_dot(i->N, r->d);
  float sin_th2 = mu1 / mu2 * sqrtf(1 - cos_th1 * cos_th1);

  RAY_VERIFY(cos_th1 > 0, "cos = %f < 0; %f %f %f; %d", cos_th1, i->N[0], i->N[1], i->N[2], i->l);
  
  RGB reflected_colour;
  process_reflection(&reflected_colour, s, r, i, max_depth, recursion_depth - 1);

  if (sin_th2 > 1) {
    glm_vec3_copy(reflected_colour, *dest);
    return;
  }

  RGB refracted_colour;
  process_refraction(&refracted_colour, s, r, i, max_depth, recursion_depth - 1, cos_th1, sin_th2, mu1, mu2); 
  float R0 = powf((mu1 - mu2) / (mu1 + mu2), 2);
  float R = R0 + (1 - R0) * powf(1 - cos_th1, 5);

  if (i->l == OUTER) {
    glm_vec3_mul(refracted_colour, i->p->colour, refracted_colour);
  }
  
  glm_vec3_scale(reflected_colour, R, reflected_colour);
  glm_vec3_scale(refracted_colour, 1 - R, refracted_colour);
  glm_vec3_add(reflected_colour, refracted_colour, *dest);
}

void process_diffuse(RGB *dest, const scene *s, const ray *r, const intsec *i) {
  process_light(s, r, i, dest);
}


void process_metal(RGB *dest, const scene *s, const ray *r, const intsec *i, float max_depth, size_t recursion_depth) {
  process_reflection(dest, s, r, i, max_depth, recursion_depth - 1);
  // Blend
  glm_vec3_mul(i->p->colour, *dest, *dest);
}


void process_ray(RGB *dest, const scene *s, const ray *r, float max_depth, size_t recursion_depth) {
  if (recursion_depth == 0) {
    (*dest)[0] = 0;
    (*dest)[1] = 0;
    (*dest)[2] = 0;
    return;
  }
  intsec res = intersect_ray(s, r, max_depth);

  if (res.p == NULL) {
    glm_vec3_copy(s->bg, *dest);
  } else { 
    switch (res.p->material) {
    case MATERIAL_DIFFUSE:
      process_diffuse(dest, s, r, &res);
      break;
    case MATERIAL_METAL:
      process_metal(dest, s, r, &res, max_depth, recursion_depth);
      break;
    case MATERIAL_DIELECTRIC:
      process_dielectric(dest, s, r, &res, max_depth, recursion_depth);
      break;
    }
  }
}

void process_ray_regular(RGB *dest, const scene *s, const ray *r, size_t recursion_depth) {
  process_ray(dest, s, r, INFINITY, recursion_depth);
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

			process_ray_regular(&render[x + y * s->cam.w], s, &r, s->rec_depth);
		}
	}

	return render;
}

