#ifndef _INCLUDE_SCENE_H_
#define _INCLUDE_SCENE_H_

#include "lib/mathc.h"
#include "sdf.h"
#include "material.h"

struct SceneInfo {
	double distance;
	struct Material mat;
};

struct SceneInfo map(double* pos);

void get_normal(double* result, double* pos);

#endif
