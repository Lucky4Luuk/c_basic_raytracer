#ifndef _INCLUDE_SCENE_H_
#define _INCLUDE_SCENE_H_

#include "mathc.h"
#include "sdf.h"

struct SceneInfo {
	mfloat_t distance;
};

struct SceneInfo map(mfloat_t* pos);

#endif
