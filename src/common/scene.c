#include "scene.h"

#include <stdio.h>
#include <malloc.h>

void print_RGB(RGB rgb) {
	printf("(%u, %u, %u)", rgb.r, rgb.g, rgb.b);
}

void print_vec3(vec3 v) {
	printf("(%.4f, %.4f, %.4f)", v[0], v[1], v[2]);
}

void print_vec4(vec4 v) {
	printf("(%.4f, %.4f, %.4f, %.4f)", v[0], v[1], v[2], v[3]);
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

		printf("\t\t\tsizes = ");
		print_vec3(b->szs[i]);
		printf("\n");

		print_comm_data(&b->comm_data, i);
	}
	printf("\n");
	printf("\n");
}

void print_ellipsoids(const ellipsoids *e) {
	printf("\tELLIPSOIDS:\n");
	for (size_t i = 0; i < e->n; i++) {
		printf("\t\tELLIPSOID #%ld:\n", i);

		printf("\t\t\tradiuses = ");
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
	printf("SCENE:\n\tbackground = ");
	print_RGB(s->bg);	
	printf("\n");
	printf("\n");

	print_cam(&s->cam);
	print_planes(&s->plns);
	print_boxes(&s->bxs);
	print_ellipsoids(&s->elps);
}

comm_data new_comm_data(size_t cnt) { 
	comm_data data = {
		.col = malloc(sizeof(*data.col) * cnt),
		.pos = malloc(sizeof(*data.pos) * cnt),
		.rot = malloc(sizeof(*data.rot) * cnt)
	};

	return data;
}

planes new_planes(size_t planes_cnt) {
	planes p = {
		.n = planes_cnt,
		.nrms = malloc(sizeof(*p.nrms) * planes_cnt),
		.comm_data = new_comm_data(planes_cnt),
	};

	return p;
}

boxes new_boxes(size_t boxes_cnt) {
	boxes b = {
		.n = boxes_cnt,
		.szs = malloc(sizeof(*b.szs) * boxes_cnt),
		.comm_data = new_comm_data(boxes_cnt),
	};
	
	return b;
}

ellipsoids new_ellipsoids(size_t ellipsoids_cnt) {
	ellipsoids e = {
		.n = ellipsoids_cnt,
		.rads = malloc(sizeof(*e.rads) * ellipsoids_cnt),
		.comm_data = new_comm_data(ellipsoids_cnt),
	};
	
	return e;
}

scene new_scene(size_t boxes_cnt, size_t ellipsoids_cnt, size_t planes_cnt) {
	scene s = {
		.bxs = new_boxes(boxes_cnt),
		.elps = new_ellipsoids(ellipsoids_cnt),
		.plns = new_planes(planes_cnt),
	};
	
	return s;
}

void free_comm_data(comm_data *data) {
	free(data->rot);
	free(data->pos);
	free(data->col);
}

void free_boxes(boxes *b) {
	free(b->szs);
	free_comm_data(&b->comm_data);
}

void free_planes(planes *p) {
	free(p->nrms);
	free_comm_data(&p->comm_data);
}

void free_ellipsoids(ellipsoids *e) {
	free(e->rads);
	free_comm_data(&e->comm_data);
}

void free_scene(scene *s) {
	free_boxes(&s->bxs);
	free_planes(&s->plns);
	free_ellipsoids(&s->elps);
}

