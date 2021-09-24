#include "camera.h"

struct Ray camera_to_ray(double* uv, double aspect, struct Camera camera) {
    /*
    float imageAspectRatio = imageWidth / (float)imageHeight; // assuming width > height
    float Px = (2 * ((x + 0.5) / imageWidth) - 1) * tan(fov / 2 * M_PI / 180) * imageAspectRatio;
    float Py = (1 - 2 * ((y + 0.5) / imageHeight) * tan(fov / 2 * M_PI / 180);
    Vec3f rayOrigin(0);
    Vec3f rayDirection = Vec3f(Px, Py, -1) - rayOrigin; // note that this just equal to Vec3f(Px, Py, -1);
    rayDirection = normalize(rayDirection); // it's a direction so don't forget to normalize
    */

    double adjusted_uv[VEC2_SIZE];
    adjusted_uv[0] = uv[0] * 2.0 - 1.0;
    adjusted_uv[1] = 1.0 - uv[1] * 2.0;

    double px = adjusted_uv[0] * tan(to_radians(camera.fov)) * aspect;
    double py = adjusted_uv[1] * tan(to_radians(camera.fov));

    struct Ray ray;
    ray.ro[0] = camera.pos[0];
    ray.ro[1] = camera.pos[1];
    ray.ro[2] = camera.pos[2];

    ray.rd[0] = px;
    ray.rd[1] = py;
    ray.rd[2] = 1.0;

    vec3_subtract(ray.rd, ray.rd, camera.pos);
    vec3_normalize(ray.rd, ray.rd);

    return ray;
}
