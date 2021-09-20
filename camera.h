#ifndef _INCLUDE_CAMERA_H_
#define _INCLUDE_CAMERA_H_

#include "raytracing.h"
#include "mathc.h"

struct Camera {
    mfloat_t pos[VEC3_SIZE];
    mfloat_t dir[VEC3_SIZE];

    mfloat_t fov;
};

struct Ray camera_to_ray(mfloat_t* uv, mfloat_t aspect, struct Camera camera);

#endif
