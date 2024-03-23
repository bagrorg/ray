#include "parser.h"

#include <errno.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#include "common/assert.h"
#include "common/scene.h"

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
static const char *NEW_LIGHT = "NEW_LIGHT";
static const char *AMBIENT_LIGHT = "AMBIENT_LIGHT";
static const char *LIGHT_INTENSITY = "LIGHT_INTENSITY";
static const char *LIGHT_DIRECTION = "LIGHT_DIRECTION";
static const char *LIGHT_POSITION = "LIGHT_POSITION";
static const char *LIGHT_ATTENUATION = "LIGHT_ATTENUATION";
static const char *RAY_DEPTH = "RAY_DEPTH";
static const char *METALLIC = "METALLIC";
static const char *DIELECTRIC = "DIELECTRIC";
static const char *IOR = "IOR";


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
	(*dst)[0] = val;

	token = strtok(NULL, " ");
	RAY_VERIFY(token != NULL, "Colour line has less than 3 arguments!");
	val = atof(token);
	(*dst)[1] = val;

	token = strtok(NULL, " ");
	RAY_VERIFY(token != NULL, "Colour line has less than 3 arguments!");
	val = atof(token);
	(*dst)[2] = val;
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

typedef struct {
  vec3 property;
  vec3 intensity;
} common_light_data;

scene parse(const char *path) {
	ssize_t res;
	size_t len;
	char *line = NULL;
	FILE *f = fopen(path, "r");	

	RAY_VERIFY(f != NULL, "Problems with opening scene file (%s)", strerror(errno));

	scene s;	
  primitive *p;
  light_point *lp;
  light_directed *ld;

  common_light_data stack_tmp;
  common_light_data *current_light_common_data;

  init(&s.primitives, sizeof *p);
  init(&s.lights_point, sizeof *lp);
  init(&s.lights_directed, sizeof *ld);

  vec3 default_pos = {0, 0, 0};
  vec4 default_rotation = {0, 0, 0, 1};

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
      p = push(&s.primitives);
      p->material = MATERIAL_DIFFUSE;
      glm_vec3_copy(default_pos, p->position);
      glm_vec4_copy(default_rotation, p->rotation);
		} else if (strncmp(PLANE, line, strlen(PLANE)) == 0) {
      RAY_VERIFY(p != NULL, "primitive is  NULL!!!");
      process_vec3(line, &p->normal); 
      p->type = PRIMITIVE_PLANE;
		} else if (strncmp(ELLIPSOID, line, strlen(ELLIPSOID)) == 0) {
      RAY_VERIFY(p != NULL, "primitive is  NULL!!!");
      process_vec3(line, &p->radius); 
      p->type = PRIMITIVE_ELLIPSOID;
		} else if (strncmp(BOX, line, strlen(BOX)) == 0) {
      RAY_VERIFY(p != NULL, "primitive is  NULL!!!");
      process_vec3(line, &p->sizes); 
      p->type = PRIMITIVE_BOX;
		} else if (strncmp(POSITION, line, strlen(POSITION)) == 0) {
      RAY_VERIFY(p != NULL, "primitive is  NULL!!!");
      process_vec3(line, &p->position); 
		} else if (strncmp(ROTATION, line, strlen(ROTATION)) == 0) {
      RAY_VERIFY(p != NULL, "primitive is  NULL!!!");
      process_vec4(line, &p->rotation); 
		} else if (strncmp(COLOUR, line, strlen(COLOUR)) == 0) {
      RAY_VERIFY(p != NULL, "primitive is  NULL!!!");
      process_colour(line, &p->colour); 
    } else if (strncmp(METALLIC, line, strlen(METALLIC)) == 0) {
      RAY_VERIFY(p != NULL, "primitive is  NULL!!!");
      p->material = MATERIAL_METAL;
    } else if (strncmp(DIELECTRIC, line, strlen(DIELECTRIC)) == 0) {
      RAY_VERIFY(p != NULL, "primitive is  NULL!!!");
      p->material = MATERIAL_DIELECTRIC;
    } else if (strncmp(IOR, line, strlen(IOR)) == 0) {
      RAY_VERIFY(p != NULL, "primitive is  NULL!!!");
      process_float(line, &p->ior); 
		} else if (strncmp(NEW_LIGHT, line, strlen(NEW_LIGHT)) == 0) {
      lp = NULL;
      ld = NULL;
      stack_tmp.property[0] = 0;
      stack_tmp.property[1] = 0;
      stack_tmp.property[2] = 0;
      stack_tmp.intensity[0] = 0;
      stack_tmp.intensity[1] = 0;
      stack_tmp.intensity[2] = 0;
      
      current_light_common_data = &stack_tmp;
		} else if (strncmp(LIGHT_INTENSITY, line, strlen(LIGHT_INTENSITY)) == 0) {
      RAY_VERIFY(current_light_common_data != NULL, "primitive is  NULL!!!");
      process_vec3(line, &current_light_common_data->intensity);
		} else if (strncmp(LIGHT_POSITION, line, strlen(LIGHT_POSITION)) == 0) {
      if (lp == NULL) {
        lp = push(&s.lights_point);
      }
      current_light_common_data = (common_light_data*)lp; 
      process_vec3(line, &lp->position);
		} else if (strncmp(LIGHT_DIRECTION, line, strlen(LIGHT_DIRECTION)) == 0) {
      ld = push(&s.lights_directed);
      current_light_common_data = (common_light_data*)ld; 
      process_vec3(line, &ld->direction);
		} else if (strncmp(LIGHT_ATTENUATION, line, strlen(LIGHT_ATTENUATION)) == 0) {
      if (lp == NULL) {
        lp = push(&s.lights_point);
      }
      current_light_common_data = (common_light_data*)lp; 
      process_vec3(line, &lp->attenuation);
		} else if (strncmp(AMBIENT_LIGHT, line, strlen(AMBIENT_LIGHT)) == 0) {
      process_colour(line, &s.ambient);
		} else if (strcmp(line, "\n") != 0) {
			RAY_WARNING("Unknown command '%s'", line);  
    }
	}
	
	fclose(f);

	return s;
}

