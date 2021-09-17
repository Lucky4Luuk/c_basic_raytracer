#ifndef _INCLUDE_RAYTRACING_H_
#define _INCLUDE_RAYTRACING_H_

#include "lib/mathc.h"
#include "scene.h"

#define EPSILON 0.01

struct Rayhit {
	bool hit;
	int steps;
	mfloat_t pos[VEC3_SIZE];
	mfloat_t normal[VEC3_SIZE];
};

struct Rayhit trace(mfloat_t* ro, mfloat_t* rd, int max_steps);

#endif
