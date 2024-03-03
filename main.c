#include <stdio.h>
#include <math.h>

#include "parser/parser.h"
#include "output/ppm.h"

int main(int argc, char **argv) {
	if (argc < 3) {
		fprintf(stderr, "Wrong number of arguments!\nUsage: ray [INPUT SCENE PATH] [OUTPUT PATH]\n");
		return 1;
	}
	
	 char *scene_path = argv[1];
	 scene s = parse(scene_path);
	 print_scene(&s);
	
	//const size_t w = 1024;
	//const size_t h = 780;
	//RGB pic[w * h];

	//for (size_t x = 0; x < w; x++) {
	//	for (size_t y = 0; y < h; y++) {
	//		RGB col;
	//		col.r = (uint8_t) (255.f * ((float)w - x) / w);
	//		col.g = (uint8_t) (255.f * ((float)h - y) / h);
	//		col.r = (uint8_t) (255.f * fabs(sin(x)));

	//		pic[x + w * y] = col;
	//	}
	//}

	//write_ppm(argv[2], w, h, pic);
}
