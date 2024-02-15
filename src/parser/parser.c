#include "parser.h"

#include <errno.h>
#include <string.h>

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

scene parse(const char *path) {
	FILE *f = fopen(path, "r");	

	RAY_VERIFY(f != NULL, "Problems with opening scene file (%s)", strerror(errno));

						
}

