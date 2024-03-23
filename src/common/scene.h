#ifndef SCENE_H
#define SCENE_H

#include "vector.h"

#include <stdint.h>
#include <stddef.h>

#include "cglm/cglm.h"

typedef vec3 RGB;

///////////// GEOMETRY OBJECTS ///////////////

typedef enum {
  PRIMITIVE_BOX,
  PRIMITIVE_PLANE,
  PRIMITIVE_ELLIPSOID
} primitive_type;

typedef enum {
  MATERIAL_DIFFUSE,
  MATERIAL_METAL,
  MATERIAL_DIELECTRIC,
} primitive_material;

typedef struct {
  primitive_type type;
  primitive_material material;

	vec3 position;	
	vec4 rotation;	
	RGB colour;			
  float ior;
  
  union {
    vec3 normal;
    vec3 radius;
    vec3 sizes;
  };

} primitive;

typedef struct {
  vec3 direction;
  vec3 intensity;
} light_directed;

typedef struct {
  vec3 position;
  vec3 intensity;
  vec3 attenuation;
} light_point;
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
  vector lights_directed;
  vector lights_point;
  
  size_t rec_depth;
} scene;


void print_scene(const scene *s);

#endif
