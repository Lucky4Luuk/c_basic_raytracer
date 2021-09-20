#include "camera.h"

struct Ray camera_to_ray(mfloat_t* uv, mfloat_t aspect, struct Camera camera) {
    mfloat_t proj_mat[MAT4_SIZE];
    mfloat_t view_mat[MAT4_SIZE];
    mfloat_t invprojview[MAT4_SIZE];

    mfloat_t near = 0.02;
    mfloat_t far = 1024.0;

    mat4_perspective(proj_mat, to_radians(camera.fov), aspect, near, far);

    mfloat_t up[VEC3_SIZE] = {0.0, 1.0, 0.0};
    mat4_look_at(view_mat, camera.pos, camera.dir, up);

    mat4_multiply(invprojview, proj_mat, view_mat);
    mat4_inverse(invprojview, invprojview);

    vec2_multiply_f(uv, uv, 2.0);
    vec2_subtract_f(uv, uv, 1.0);
    mfloat_t adjusted_pos[VEC4_SIZE]     = {uv[0], uv[1], -1.0, 1.0};
    vec2_multiply_f(uv, uv, (far - near));
    mfloat_t adjusted_pos_alt[VEC4_SIZE] = {uv[0], uv[1], far + near, far - near};

    struct Ray ray;
    //Technically I have to determine this based on my near value
    //Formula: (invprojview * adjusted_pos * near).xyz
    ray.ro[0] = camera.pos[0];
    ray.ro[1] = camera.pos[1];
    ray.ro[2] = camera.pos[2];

    mfloat_t tmp[VEC4_SIZE];
    vec4_multiply_mat4(tmp, adjusted_pos_alt, invprojview);

    ray.rd[0] = tmp[0];
    ray.rd[1] = tmp[1];
    ray.rd[2] = tmp[2];

    return ray;
}
