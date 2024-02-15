#include "scene.h"

#include <stdio.h>

void print_RGB(RGB rgb) {
	printf("(%u, %u, %u)", rgb.r, rgb.g, rgb.b);
}

void print_vec3(vec3 v) {
	printf("(%.4f, %.4f, %.4f)", v.x, v.y, v.z);
}

void print_vec4(vec4 v) {
	printf("(%.4f, %.4f, %.4f, %.4f)", v.x, v.y, v.z, v.w);
}

void print_comm_data(const comm_data *data, size_t i) {
	printf("\t\t\tposition = ");
	print_vec3(data->pos[i]);
	printf("\n");

	printf("\t\t\trotation = ");
	print_vec4(data->rot[i]);
	printf("\n");

	printf("\t\t\tcolor = ");
	print_RGB(data->col[i]);
	printf("\n");
}

void print_planes(const planes *p) {
	printf("\tPLANES:\n");
	for (size_t i = 0; i < p->n; i++) {
		printf("\t\tPLANE #%ld:\n", i);

		printf("\t\t\tnormal = ");
		print_vec3(p->nrms[i]);
		printf("\n");

		print_comm_data(&p->comm_data, i);
	}
	printf("\n");
	printf("\n");
}

void print_boxes(const boxes *b) {
	printf("\tBOXES:\n");
	for (size_t i = 0; i < b->n; i++) {
		printf("\t\tBOX #%ld:\n", i);

		printf("\t\t\tnormal = ");
		print_vec3(b->szs[i]);
		printf("\n");

		print_comm_data(&b->comm_data, i);
	}
	printf("\n");
	printf("\n");
}

void print_ellipsoids(const ellipsoids *e) {
	printf("\tPLANES:\n");
	for (size_t i = 0; i < e->n; i++) {
		printf("\t\tPLANE #%ld:\n", i);

		printf("\t\t\tnormal = ");
		print_vec3(e->rads[i]);
		printf("\n");

		print_comm_data(&e->comm_data, i);
	}
	printf("\n");
	printf("\n");
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

void print_scene(const scene *s) {
	printf("SCENE:\n\t background = ");
	print_RGB(s->bg);	
	printf("\n");
	printf("\n");

	print_cam(&s->cam);
	print_planes(&s->plns);
	print_boxes(&s->bxs);
	print_ellipsoids(&s->elps);
}


