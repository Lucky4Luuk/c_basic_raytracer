#include "lighting.h"

void shade_point(double* result, double* light, double* normal, double* albedo) {
    double zero[VEC3_SIZE] = {0.0, 0.0, 0.0};
    double one[VEC3_SIZE] = {1.0, 1.0, 1.0};
    double light_copy[VEC3_SIZE]; //To avoid changing the original input light
    // vec3_subtract(light_copy, one, light);
    light_copy[0] = light[0]; light_copy[1] = light[1]; light_copy[2] = light[2];
    double n_dot_l = vec3_dot(light_copy, normal);
    vec3_multiply_f(result, albedo, n_dot_l);
    vec3_max(result, result, zero);
    vec3_min(result, result, one);
}
