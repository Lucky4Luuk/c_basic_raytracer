#include <stdio.h>

#define MATHC_USE_INT32
#define MATHC_USE_DOUBLE_FLOATING_POINT

#include "mathc.h"
#include "display.h"

#include "camera.h"
#include "raytracing.h"

int main() {
	createScreen();

	struct Camera camera;
	camera.pos[0] = 0.0;
	camera.pos[1] = 0.0;
	camera.pos[2] = -5.0;

	camera.dir[0] = 0.0;
	camera.dir[1] = 0.0;
	camera.dir[2] = 1.0;

	camera.fov = 60.0;

	for (int ix = 0; ix < SCREEN_WIDTH; ++ix) {
		for (int iy = 0; iy < SCREEN_HEIGHT; ++iy) {
			mfloat_t uv[VEC2_SIZE];
			uv[0] = ((mfloat_t)ix) / ((mfloat_t)SCREEN_WIDTH);
			uv[1] = ((mfloat_t)iy) / ((mfloat_t)SCREEN_HEIGHT);

			mfloat_t aspect = ((mfloat_t)SCREEN_WIDTH) / ((mfloat_t)SCREEN_HEIGHT);

			struct Ray ray = camera_to_ray(uv, aspect, camera);
			struct Rayhit hit = trace(ray, 255);

			// uint8_t r = hit.steps;
			// uint8_t g = hit.steps;
			// uint8_t b = hit.steps;
			uint8_t r = (uint8_t)(ray.rd[0] * 255);
			uint8_t g = (uint8_t)(ray.rd[1] * 255);
			uint8_t b = (uint8_t)(ray.rd[2] * 255);
			setPixel(ix, iy, r, g, b);
		}
	}
	flushScreen();

	deleteScreen();
	return 0;
}
