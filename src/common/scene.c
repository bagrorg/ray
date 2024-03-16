#include "scene.h"

#include <stdio.h>
#include <malloc.h>

void print_RGB(RGB rgb) {
	printf("(%f, %f, %f)", rgb[0], rgb[1], rgb[2]);
}

void print_vec3(const vec3 v) {
	printf("(%.4f, %.4f, %.4f)", v[0], v[1], v[2]);
}

void print_vec4(const vec4 v) {
	printf("(%.4f, %.4f, %.4f, %.4f)", v[0], v[1], v[2], v[3]);
}

void print_common_data(const primitive *p) {
	printf("\t\t\tposition = ");
	print_vec3(p->position);
	printf("\n");

	printf("\t\t\trotation = ");
	print_vec4(p->rotation);
	printf("\n");

	printf("\t\t\tcolor = ");
	print_RGB(p->colour);
	printf("\n");
}

void print_planes(const primitive *p) {
	printf("\tPLANE:\n");
  print_common_data(p);  
  printf("\t\t\tnormal = ");
  print_vec3(p->normal);
  printf("\n\n");
}

void print_boxes(const primitive *p) {
	printf("\tBOXES:\n");
  print_common_data(p);  
  printf("\t\t\tsizes = ");
  print_vec3(p->sizes);
  printf("\n\n");
}

void print_ellipsoids(const primitive *p) {
	printf("\tELLIPSOIDS:\n");
  print_common_data(p);  
  printf("\t\t\tradiuses = ");
  print_vec3(p->radius);
  printf("\n\n");
}

void print_cam(const camera *c) {
	printf("\tCAMERA:\n");

	printf("\t\tposition = ");
	print_vec3(c->pos);
	printf("\n");

	printf("\t\tright = ");
	print_vec3(c->right);
	printf("\n");

	printf("\t\tup = ");
	print_vec3(c->up);
	printf("\n");

	printf("\t\tforward = ");
	print_vec3(c->frwrd);
	printf("\n");

	printf("\t\tfov_x = %.4f\n", c->fov_x);

	printf("\t\twidth = %ld\n", c->w);
	printf("\t\theight = %ld\n", c->h);
	
	printf("\n");
	printf("\n");
}

void print_light_directed(const light_directed *ld) {
  printf("\tDIRECTED LIGHT:\n");

  printf("\t\tdirection = ");
  print_vec3(ld->direction);
  printf("\n");

  printf("\t\tintensity = ");
  print_vec3(ld->intensity);
  printf("\n");
  printf("\n");
}

void print_light_point(const light_point *lp) {
  printf("\tPOINT LIGHT:\n");

  printf("\t\tposition = ");
  print_vec3(lp->position);
  printf("\n");

  printf("\t\tintensity = ");
  print_vec3(lp->intensity);
  printf("\n");

  printf("\t\tattenuation = ");
  print_vec3(lp->attenuation);
  printf("\n");

  printf("\n");
}

void print_scene(const scene *s) {
	printf("SCENE:\n\tbackground = ");
	print_RGB(s->bg);	
	printf("\n");
	printf("\n");

	print_cam(&s->cam);
  
  for (size_t i = 0; i < s->primitives.size; i++) {
    primitive *p = &((primitive*)s->primitives.data)[i];

    switch (p->type) {
      case PRIMITIVE_BOX:
        print_boxes(p);
        break;
      case PRIMITIVE_PLANE:
        print_planes(p);
        break;
      case PRIMITIVE_ELLIPSOID:
        print_ellipsoids(p);
        break;
    }
  }

  for (size_t i = 0; i < s->lights_point.size; i++) {
    print_light_point(&((light_point*)s->lights_point.data)[i]);
  }

  for (size_t i = 0; i < s->lights_directed.size; i++) {
    print_light_directed(&((light_directed*)s->lights_directed.data)[i]);
  }
}

