#ifndef SCENE_H
#define SCENE_H

#include "vector.h"

#include <stdint.h>
#include <stddef.h>

#include "cglm/cglm.h"

typedef struct {
  uint8_t r, g, b;
} RGB;

///////////// GEOMETRY OBJECTS ///////////////

typedef enum {
  PRIMITIVE_BOX,
  PRIMITIVE_PLANE,
  PRIMITIVE_ELLIPSOID
} primitive_type;

typedef struct {
  primitive_type type;

	vec3 position;	
	vec4 rotation;	
	RGB colour;			
  
  union {
    vec3 normal;
    vec3 radius;
    vec3 sizes;
  };

} primitive;

typedef struct {
  size_t n;
  vec3 *dirs;
  vec3 *ints;
} lights_dir;

typedef struct {
  size_t n;
  vec3 *pos;
  vec3 *ints;
  vec3 *attens;
} lights_point;
///////////// //////////////// ///////////////

/////////////	    CAMERA	   ///////////////
typedef struct {
	vec3 pos;

	vec3 right;
	vec3 up;
	vec3 frwrd;

	float fov_x;

	size_t w;
	size_t h;
} camera;
///////////// //////////////// ///////////////

typedef struct {
	RGB bg;
  RGB ambient;
	camera cam;
  
  vector primitives;

  lights_point lp;
  lights_dir ld;
} scene;


void print_scene(const scene *s);

#endif
