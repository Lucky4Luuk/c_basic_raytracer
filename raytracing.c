#include "raytracing.h"

const mfloat_t EPSILON = 0.01;

struct Rayhit trace(struct Ray ray, int max_steps) {
	struct Rayhit hit;
	hit.hit = false;
	hit.steps = 0;

	mfloat_t pos[VEC3_SIZE];
	pos[0] = ray.ro[0];
	pos[1] = ray.ro[1];
	pos[2] = ray.ro[2];

	for (int step = 0; step < max_steps; ++step) {
		struct SceneInfo map_info = map(pos);
		if (map_info.distance < EPSILON) {
			hit.hit = true;
			//I hate not being able to copy arrays by just assigning
			hit.pos[0] = pos[0];
			hit.pos[1] = pos[1];
			hit.pos[2] = pos[2];
			break;
		}
		vec3_add(pos, pos, vec3_multiply_f(ray.rd, ray.rd, map_info.distance));
		hit.steps = step; //I think copying the value of step is cheaper than incrementing the value of hit.steps but it really doesn't matter
	}

	return hit;
}
