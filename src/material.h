#ifndef _INCLUDE_MATERIAL_H_
#define _INCLUDE_MATERIAL_H_

#include "lib/mathc.h"

struct Material {
    double albedo[VEC3_SIZE];
    double emissive;
};

#endif
