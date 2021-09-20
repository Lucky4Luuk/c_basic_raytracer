#ifndef _INCLUDE_RAYTRACING_H_
#define _INCLUDE_RAYTRACING_H_

#include "mathc.h"
#include "scene.h"

struct Ray {
	mfloat_t ro[VEC3_SIZE];
	mfloat_t rd[VEC3_SIZE];
};

struct Rayhit {
	bool hit;
	int steps;
	mfloat_t pos[VEC3_SIZE];
	mfloat_t normal[VEC3_SIZE];
};

struct Rayhit trace(struct Ray ray, int max_steps);

#endif
