#include "raytracing.h"

const double EPSILON = 0.1;

struct Rayhit trace(struct Ray ray, int max_steps) {
	struct Rayhit hit;
	hit.hit = false;
	hit.steps = 0;

	double pos[VEC3_SIZE];
	pos[0] = ray.ro[0];
	pos[1] = ray.ro[1];
	pos[2] = ray.ro[2];

	double dir[VEC3_SIZE];
	dir[0] = ray.rd[0];
	dir[1] = ray.rd[1];
	dir[2] = ray.rd[2];

	double tmp[VEC3_SIZE];

	for (int step = 0; step < max_steps; ++step) {
		struct SceneInfo map_info = map(pos);
		if (map_info.distance < EPSILON) {
			vec3_multiply_f(tmp, dir, map_info.distance);
			vec3_add(pos, pos, tmp);

			hit.hit = true;
			//I hate not being able to copy arrays by just assigning
			hit.pos[0] = pos[0];
			hit.pos[1] = pos[1];
			hit.pos[2] = pos[2];
			hit.mat = map_info.mat;
			break;
		}
		vec3_multiply_f(tmp, dir, map_info.distance);
		vec3_add(pos, pos, tmp);
		hit.steps = step; //I think copying the value of step is cheaper than incrementing the value of hit.steps but it really doesn't matter
	}

	return hit;
}
