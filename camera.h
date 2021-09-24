#ifndef _INCLUDE_CAMERA_H_
#define _INCLUDE_CAMERA_H_

#include "raytracing.h"
#include "mathc.h"

struct Camera {
    double pos[VEC3_SIZE];
    double dir[VEC3_SIZE];

    double fov;
};

struct Ray camera_to_ray(double* uv, double aspect, struct Camera camera);

#endif
