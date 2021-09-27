#ifndef _INCLUDE_LIGHTING_H_
#define _INCLUDE_LIGHTING_H_

#include "lib/mathc.h"

void shade_point(double* result, double* light, double* normal, double* albedo);

#endif
