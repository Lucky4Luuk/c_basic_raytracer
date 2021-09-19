#include "sdf.h"

mfloat_t sdfSphere(mfloat_t* pos, mfloat_t radius) {
	return vec3_length(pos) - radius;
}
