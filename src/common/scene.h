#ifndef SCENE_H
#define SCENE_H

#include <stdint.h>
#include <stddef.h>

#include "cglm/cglm.h"

typedef struct {
	uint8_t r;
	uint8_t g;
	uint8_t b;
} RGB;

///////////// GEOMETRY OBJECTS ///////////////
typedef struct {
	vec3 *pos;	
	vec4 *rot;	
	RGB *col;			
} comm_data;

typedef struct {
	size_t n;
	vec3 *nrms;
	comm_data comm_data;
} planes;

typedef struct {
	size_t n;
	vec3 *rads;
	comm_data comm_data;
} ellipsoids;

typedef struct {
	size_t n;
	vec3 *szs;
	comm_data comm_data;
} boxes;
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
	camera cam;

	planes plns;
	ellipsoids elps;
	boxes bxs;
} scene;

scene new_scene(size_t boxes, size_t ellipses, size_t planes);
void free_scene(scene *s);

void print_scene(const scene *s);

#endif
