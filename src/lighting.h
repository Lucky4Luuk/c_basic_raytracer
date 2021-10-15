#ifndef _INCLUDE_LIGHTING_H_
#define _INCLUDE_LIGHTING_H_

#include <stdlib.h>
#include "lib/mathc.h"
#include "raytracing.h"
#include "material.h"

void shade_point(int depth, double* hit_pos, double* result, double* light, double* normal, struct Material mat);

#endif
