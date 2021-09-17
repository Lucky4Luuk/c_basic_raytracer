#include <stdio.h>

#include "lib/mathc.h"
#include "display.h"

#include "raytracing.h"

int main() {
	createScreen();

	for (int ix = 0; ix < SCREEN_WIDTH; ++ix) {
		for (int iy = 0; iy < SCREEN_HEIGHT; ++iy) {
			uint8_t r = ix;
			uint8_t g = iy;
			uint8_t b = 0;
			setPixel(ix, iy, r, g, b);
		}
	}
	flushScreen();

	deleteScreen();
	return 0;
}
