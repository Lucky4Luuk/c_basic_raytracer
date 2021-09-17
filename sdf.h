#ifndef _INCLUDE_SDF_H_
#define _INCLUDE_SDF_H_

#include "lib/mathc.h"

mfloat_t sdfSphere(mfloat_t* pos, mfloat_t radius) {
	return vec3_length(pos) - radius;
}

#endif
