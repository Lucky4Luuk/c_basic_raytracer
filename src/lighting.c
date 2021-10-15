#include "lighting.h"

#ifndef MAX_BOUNCES
#define MAX_BOUNCES 4
#endif

void shade_point(int depth, double* hit_pos, double* result, double* light, double* normal, struct Material mat) {
    if (depth > MAX_BOUNCES) {
        result[0] = 0.0;
        result[1] = 0.0;
        result[2] = 0.0;
        return;
    }

    //Direct (basically emission)
    double emission[VEC3_SIZE] = {mat.albedo[0], mat.albedo[1], mat.albedo[2]};
    vec3_multiply_f(emission, emission, mat.emissive);
    vec3_add(result, result, emission);

    //Indirect
    double indirect[VEC3_SIZE] = {0.0, 0.0, 0.0};

    //Generate random ray
    float x = ((float)rand()/(float)(RAND_MAX));
    float y = ((float)rand()/(float)(RAND_MAX));
    float z = ((float)rand()/(float)(RAND_MAX));
    double indirect_dir[VEC3_SIZE] = {x, y, z};
    vec3_normalize(indirect_dir, indirect_dir);
    double id_dot_n = vec3_dot(normal, indirect_dir);
    vec3_multiply_f(indirect_dir, indirect_dir, (id_dot_n) ? 1.0 : -1.0);

    struct Ray indirect_ray;
    indirect_ray.ro[0] = hit_pos[0];
    indirect_ray.ro[1] = hit_pos[1];
    indirect_ray.ro[2] = hit_pos[2];
    indirect_ray.rd[0] = indirect_dir[0];
    indirect_ray.rd[1] = indirect_dir[1];
    indirect_ray.rd[2] = indirect_dir[2];

    struct Rayhit indirect_hit = trace(indirect_ray, 128);
    if (indirect_hit.hit) {
        double hit_normal[VEC3_SIZE];
        get_normal(hit_normal, indirect_hit.pos);
        double hit_albedo[VEC3_SIZE] = {1.0, 1.0, 1.0};
        double ro[VEC3_SIZE] = {indirect_hit.pos[0], indirect_hit.pos[1], indirect_hit.pos[2]};
        vec3_add(ro, ro, hit_normal);
        shade_point(depth + 1, ro, indirect, indirect_dir, hit_normal, indirect_hit.mat);
    }

    double n_dot_l = vec3_dot(indirect_dir, normal);
    vec3_multiply_f(indirect, indirect, -n_dot_l);
    vec3_add(result, result, indirect);
}
