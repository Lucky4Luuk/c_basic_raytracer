#include "scene.h"

struct SceneInfo map(mfloat_t* pos) {
	//TMP: Hardcoded test scene
	//TODO: Create scene from main.c with functions
	mfloat_t d = sdfSphere(pos, 0.5);

	struct SceneInfo si;
	si.distance = d;

	return si;
}
