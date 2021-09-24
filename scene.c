#include "scene.h"

struct SceneInfo map(double* pos) {
	//TMP: Hardcoded test scene
	//TODO: Create scene from main.c with functions
	double d = sdfSphere(pos, 0.5);

	struct SceneInfo si;
	si.distance = d;

	return si;
}
