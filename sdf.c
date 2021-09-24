#include "sdf.h"

double sdfSphere(double* pos, double radius) {
	return vec3_length(pos) - radius;
}
