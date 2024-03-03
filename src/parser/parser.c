#include "parser.h"

#include <errno.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#include "common/assert.h"

static const char *DIMENSIONS = "DIMENSIONS";
static const char *BG_COLOUR = "BG_COLOR";
static const char *CAMERA_POSITION = "CAMERA_POSITION";
static const char *CAMERA_RIGHT = "CAMERA_RIGHT";
static const char *CAMERA_UP = "CAMERA_UP";
static const char *CAMERA_FORWARD = "CAMERA_FORWARD";
static const char *CAMERA_FOV_X = "CAMERA_FOV_X";
static const char *NEW_PRIMITIVE = "NEW_PRIMITIVE";
static const char *PLANE = "PLANE";
static const char *ELLIPSOID = "ELLIPSOID";
static const char *BOX = "BOX";
static const char *POSITION = "POSITION";
static const char *ROTATION = "ROTATION";
static const char *COLOUR = "COLOR";

void process_dimensions(char *line, scene *s) {
	char *token = strtok(line, " ");
	
	RAY_VERIFY(token != NULL, "Dimensions line has less than 2 arguments!");

	token = strtok(NULL, " ");
	RAY_VERIFY(token != NULL, "Dimensions line has less than 2 arguments!");
	s->cam.w = atoll(token);

	token = strtok(NULL, " ");
	RAY_VERIFY(token != NULL, "Dimensions line has less than 2 arguments!");
	s->cam.h = atoll(token);
}

void process_float(char *line, float *dst) {
	char *token = strtok(line, " ");
	
	RAY_VERIFY(token != NULL, "Float line has less than 1 argument!");

	token = strtok(NULL, " ");
	RAY_VERIFY(token != NULL, "Float line has less than 1 argument!");
	*dst = atof(token);
}

void process_colour(char *line, RGB *dst) {
	float val;
	char *token = strtok(line, " ");

	RAY_VERIFY(token != NULL, "Colour line has less than 3 arguments!");

	token = strtok(NULL, " ");
	RAY_VERIFY(token != NULL, "Colour line has less than 3 arguments!");
	val = atof(token);
	dst->r = round(val * 255);

	token = strtok(NULL, " ");
	RAY_VERIFY(token != NULL, "Colour line has less than 3 arguments!");
	val = atof(token);
	dst->g = round(val * 255);

	token = strtok(NULL, " ");
	RAY_VERIFY(token != NULL, "Colour line has less than 3 arguments!");
	val = atof(token);
	dst->b = round(val * 255);
}

void process_vec3(char *line, vec3 *dst) {
	char *token = strtok(line, " ");

	RAY_VERIFY(token != NULL, "vec3 line has less than 3 arguments!");

	token = strtok(NULL, " ");
	RAY_VERIFY(token != NULL, "vec3 line has less than 3 arguments!");
	(*dst)[0] = atof(token);

	token = strtok(NULL, " ");
	RAY_VERIFY(token != NULL, "vec3 line has less than 3 arguments!");
	(*dst)[1] = atof(token);

	token = strtok(NULL, " ");
	RAY_VERIFY(token != NULL, "vec3 line has less than 3 arguments!");
	(*dst)[2] = atof(token);
}

void process_vec4(char *line, vec4 *dst) {
	char *token = strtok(line, " ");

	RAY_VERIFY(token != NULL, "vec4 line has less than 4 arguments!");

	token = strtok(NULL, " ");
	RAY_VERIFY(token != NULL, "vec4 line has less than 4 arguments!");
	(*dst)[0] = atof(token);

	token = strtok(NULL, " ");
	RAY_VERIFY(token != NULL, "vec4 line has less than 4 arguments!");
	(*dst)[1] = atof(token);

	token = strtok(NULL, " ");
	RAY_VERIFY(token != NULL, "vec4 line has less than 4 arguments!");
	(*dst)[2] = atof(token);

	token = strtok(NULL, " ");
	RAY_VERIFY(token != NULL, "vec4 line has less than 4 arguments!");
	(*dst)[3] = atof(token);
}

void prealloc_data(FILE *f, scene *s) {
	char *line = NULL;
	size_t len;
	ssize_t res;

	size_t boxes_cnt = 0;
	size_t ellipsoids_cnt = 0;
	size_t planes_cnt = 0;

	while ((res = getline(&line, &len, f)) != -1) {
		if (strncmp(BOX, line, strlen(BOX)) == 0) {
			boxes_cnt++;
		} else if (strncmp(PLANE, line, strlen(PLANE)) == 0) {
			ellipsoids_cnt++;
		} else if (strncmp(ELLIPSOID, line, strlen(ELLIPSOID)) == 0) {
			planes_cnt++;
		}
	}
	
	*s = new_scene(boxes_cnt, ellipsoids_cnt, planes_cnt);

	fseek(f, 0, SEEK_SET);
}

scene parse(const char *path) {
	ssize_t res;
	size_t len;
	char *line = NULL;
	FILE *f = fopen(path, "r");	

	RAY_VERIFY(f != NULL, "Problems with opening scene file (%s)", strerror(errno));

	scene s;	
	vec3 pos;
	RGB col;
	vec4 rot;

	size_t cur_plane = 0;
	size_t cur_box = 0;
	size_t cur_ellipsoid = 0;

	size_t *cur = NULL;
	comm_data *cur_data = NULL;

	prealloc_data(f, &s);

	while ((res = getline(&line, &len, f)) != -1) {
		if (strncmp(DIMENSIONS, line, strlen(DIMENSIONS)) == 0) {
			process_dimensions(line, &s);
		} else if (strncmp(BG_COLOUR, line, strlen(BG_COLOUR)) == 0) {
			process_colour(line, &s.bg);
		} else if (strncmp(CAMERA_POSITION, line, strlen(CAMERA_POSITION)) == 0) {
			process_vec3(line, &s.cam.pos);
		} else if (strncmp(CAMERA_RIGHT, line, strlen(CAMERA_RIGHT)) == 0) {
			process_vec3(line, &s.cam.right);
		} else if (strncmp(CAMERA_UP, line, strlen(CAMERA_UP)) == 0) {
			process_vec3(line, &s.cam.up);
		} else if (strncmp(CAMERA_FORWARD, line, strlen(CAMERA_FORWARD)) == 0) {
			process_vec3(line, &s.cam.frwrd);
		} else if (strncmp(CAMERA_FOV_X, line, strlen(CAMERA_FOV_X)) == 0) {
			process_float(line, &s.cam.fov_x);
		} else if (strncmp(NEW_PRIMITIVE, line, strlen(NEW_PRIMITIVE)) == 0) {
			if (cur != NULL) {
				memcpy(&cur_data->pos[*cur], &pos, sizeof(pos));
				memcpy(&cur_data->rot[*cur], &rot, sizeof(rot));

				cur_data->col[*cur] = col;
				*cur += 1;

				cur = NULL;
				cur_data = NULL;
			} 

			pos[0] = 0;
			pos[1] = 0;
			pos[2] = 0;

			col.r = 0;
			col.g = 0;
			col.b = 0;

			rot[0] = 0;
			rot[1] = 0;
			rot[2] = 0;
			rot[3] = 1;
		} else if (strncmp(PLANE, line, strlen(PLANE)) == 0) {
			cur = &cur_plane;
			cur_data = &s.plns.comm_data;
			process_vec3(line, &s.plns.nrms[cur_plane]);
		} else if (strncmp(ELLIPSOID, line, strlen(ELLIPSOID)) == 0) {
			cur = &cur_ellipsoid;
			cur_data = &s.elps.comm_data;
			process_vec3(line, &s.elps.rads[cur_ellipsoid]);
		} else if (strncmp(BOX, line, strlen(BOX)) == 0) {
			cur = &cur_box;
			cur_data = &s.bxs.comm_data;
			process_vec3(line, &s.bxs.szs[cur_box]);
		} else if (strncmp(POSITION, line, strlen(POSITION)) == 0) {
			process_vec3(line, &pos);
		} else if (strncmp(ROTATION, line, strlen(ROTATION)) == 0) {
			process_vec4(line, &rot);
		} else if (strncmp(COLOUR, line, strlen(COLOUR)) == 0) {
			process_colour(line, &col);
		} else if (strcmp(line, "\n") != 0) {
			RAY_WARNING("Unknown command '%s'", line);
		}
	}
	
	RAY_VERIFY(cur != NULL, "Corrupted scene file! No primitives provided!");

	memcpy(&cur_data->pos[*cur], &pos, sizeof(pos));
	memcpy(&cur_data->rot[*cur], &rot, sizeof(rot));
	*cur += 1;

	RAY_VERIFY(cur_box == s.bxs.n, "Parsed %ld boxes, but %ld expected!", cur_box, s.bxs.n);
	RAY_VERIFY(cur_ellipsoid == s.elps.n, "Parsed %ld ellipsoids, but %ld expected!", cur_ellipsoid, s.elps.n);
	RAY_VERIFY(cur_plane == s.plns.n, "Parsed %ld planes, but %ld expected!", cur_plane, s.plns.n);

	fclose(f);

	return s;
}

