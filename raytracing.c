#include "raytracing.h"

struct Rayhit trace(mfloat_t* ro, mfloat_t* rd, int max_steps) {
	struct Rayhit hit;
	hit.hit = false;
	hit.steps = 0;

	mfloat_t pos[VEC3_SIZE];
	pos[0] = ro[0];
	pos[1] = ro[1];
	pos[2] = ro[2];

	for (int step = 0; step < max_steps; ++step) {
		if (map(pos).distance < EPSILON) {
			hit.hit = true;
			//I hate not being able to copy arrays by just assigning
			hit.pos[0] = pos[0];
			hit.pos[1] = pos[1];
			hit.pos[2] = pos[2];
			break;
		}
		hit.steps = step; //I think copying the value of step is cheaper than incrementing the value of hit.steps but it really doesn't matter
	}

	return hit;
}
