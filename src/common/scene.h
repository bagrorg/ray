#ifndef SCENE_H
#define SCENE_H

#include <stdint.h>
#include <stddef.h>

typedef struct {
	uint8_t r;
	uint8_t g;
	uint8_t b;
} RGB;

typedef struct {
	float x;
	float y;
	float z;
} vec3;

typedef struct {
	float x;
	float y;
	float z;
	float w;
} vec4;


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

void print_scene(const scene *s);

#endif
