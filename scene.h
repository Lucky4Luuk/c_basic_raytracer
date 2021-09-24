#ifndef _INCLUDE_SCENE_H_
#define _INCLUDE_SCENE_H_

#include "mathc.h"
#include "sdf.h"

struct SceneInfo {
	double distance;
};

struct SceneInfo map(double* pos);

#endif
