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
  const primitive *p;

  vec3 N;
  float t;
  loc l; 
} intsec;

RGB *render(const scene *s);

#endif
