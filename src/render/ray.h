#ifndef RAY_H
#define RAY_H

#include <stdbool.h>
#include "common/scene.h"

typedef struct {
	vec3 o;			// origin
	vec3 d;			// direction
} ray;

typedef enum {
  INNER,
  OUTER
} loc;

typedef struct {
  RGB col;
  vec3 N;
  float t;
  loc l; 
  bool succ;
} intsec;

void process_ray(intsec *res_glob, const scene *s, const ray *r, float max_depth);
RGB *render(const scene *s);

#endif
