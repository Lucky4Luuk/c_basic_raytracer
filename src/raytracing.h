#ifndef _INCLUDE_RAYTRACING_H_
#define _INCLUDE_RAYTRACING_H_

#include "lib/mathc.h"
#include "scene.h"
#include "material.h"

struct Ray {
	double ro[VEC3_SIZE];
	double rd[VEC3_SIZE];
};

struct Rayhit {
	bool hit;
	int steps;
	double pos[VEC3_SIZE];
	double normal[VEC3_SIZE];
	struct Material mat;
};

struct Rayhit trace(struct Ray ray, int max_steps);

#endif
