#include "ppm.h"

#include <stdio.h>
#include <errno.h>
#include <string.h>

#include "common/assert.h"

static const float GAMMA = 2.2;
static const float IGAMMA = 1 / GAMMA;

typedef struct {
  uint8_t r, g, b;
} RGBBytes;

void aces_tonemap(vec3 brightness, vec3 dest) {
  float a = 2.51f;
  float b = 0.03f;
  float c = 2.43f;
  float d = 0.59f;
  float e = 0.14f;

  vec3 f, s; 
  glm_vec3_scale(brightness, a, f);
  glm_vec3_adds(f, b, f);
  glm_vec3_mul(f, brightness, f);

  glm_vec3_scale(brightness, c, s);
  glm_vec3_adds(s, d, s);
  glm_vec3_mul(s, brightness, s);
  glm_vec3_adds(s, e, s);

  glm_vec3_div(f, s, dest);

  glm_vec3_clamp(dest, 0, 1);
}

void colour_processing(RGB rgb, RGBBytes *dest) {
  aces_tonemap(rgb, rgb);
  
  dest->r = (uint8_t) roundf(powf(rgb[0], IGAMMA) * 255);
  dest->g = (uint8_t) roundf(powf(rgb[1], IGAMMA) * 255);
  dest->b = (uint8_t) roundf(powf(rgb[2], IGAMMA) * 255);
}

void write_ppm(const char *path, size_t w, size_t h, const RGB *data) {
	FILE *f = fopen(path, "w");

	RAY_VERIFY(f != NULL, "Failed to open %s for output (%s)", path, strerror(errno));

	fprintf(f, "P6\n");
	fprintf(f, "%ld %ld\n", w, h);
	fprintf(f, "255\n");

  RGBBytes *bytes_data = malloc(sizeof(RGBBytes) * w * h);
  for (size_t i = 0; i < w * h; i++) {
    colour_processing(data[i], &bytes_data[i]);
  }
	fwrite(bytes_data, sizeof(*bytes_data), w * h, f);
  free(bytes_data);

	fclose(f);
}


