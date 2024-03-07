#include <stdio.h>

#include "parser/parser.h"
#include "output/ppm.h"
#include "render/ray.h"

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "Wrong number of arguments!\nUsage: ray [INPUT SCENE PATH] [OUTPUT PATH]\n");
    return 1;
  }

  char *scene_path = argv[1];
  scene s = parse(scene_path);
  print_scene(&s);

  RGB *rgb = render(&s);

  write_ppm(argv[2], s.cam.w, s.cam.h, rgb);
  free(rgb);
}
