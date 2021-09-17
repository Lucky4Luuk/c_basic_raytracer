#ifndef _INCLUDE_SCENE_H_
#define _INCLUDE_SCENE_H_

#include "lib/mathc.h"
#include "sdf.h"

struct SceneInfo {
	mfloat_t distance;
};

struct SceneInfo map(mfloat_t* pos) {
	//TMP: Hardcoded test scene
	//TODO: Create scene from main.c with functions
	mfloat_t d = sdfSphere(pos, 0.5);

	struct SceneInfo si;
	si.distance = d;

	return si;
}

#endif
