#include "scene.h"

struct SceneInfo map(double* pos) {
	//TMP: Hardcoded test scene
	//TODO: Create scene from main.c with functions
	double d = sdfSphere(pos, 1.0);

	struct SceneInfo si;
	si.distance = d;

	return si;
}

/*
vec3 calcNormal( in vec3 & p ) // for function f(p)
{
    const float h = 0.0001; // replace by an appropriate value
    const vec2 k = vec2(1,-1);
    return normalize( k.xyy*f( p + k.xyy*h ) +
                      k.yyx*f( p + k.yyx*h ) +
                      k.yxy*f( p + k.yxy*h ) +
                      k.xxx*f( p + k.xxx*h ) );
}
*/

void get_normal(double* result, double* pos) {
	const double h = 0.001;
	const double k[VEC2_SIZE] = {1.0, -1.0};

	double a[VEC3_SIZE] = {k[0], k[1], k[1]};
	double b[VEC3_SIZE] = {k[1], k[1], k[0]};
	double c[VEC3_SIZE] = {k[1], k[0], k[1]};
	double d[VEC3_SIZE] = {k[0], k[0], k[0]};

	result[0] = 0.0;
	result[1] = 0.0;
	result[2] = 0.0;
	double tmp[VEC3_SIZE] = {0.0, 0.0, 0.0};
	double tmp_pos[VEC3_SIZE] = {0.0, 0.0, 0.0};

	vec3_multiply_f(tmp_pos, a, h);
	vec3_add(tmp_pos, tmp_pos, pos);
	double dist = map(tmp_pos).distance;
	vec3_multiply_f(tmp, a, dist);
	vec3_add(result, result, tmp);

	vec3_multiply_f(tmp_pos, b, h);
	vec3_add(tmp_pos, tmp_pos, pos);
	dist = map(tmp_pos).distance;
	vec3_multiply_f(tmp, b, dist);
	vec3_add(result, result, tmp);

	vec3_multiply_f(tmp_pos, c, h);
	vec3_add(tmp_pos, tmp_pos, pos);
	dist = map(tmp_pos).distance;
	vec3_multiply_f(tmp, c, dist);
	vec3_add(result, result, tmp);

	vec3_multiply_f(tmp_pos, d, h);
	vec3_add(tmp_pos, tmp_pos, pos);
	dist = map(tmp_pos).distance;
	vec3_multiply_f(tmp, d, dist);
	vec3_add(result, result, tmp);

	vec3_normalize(result, result);
}

/*
vec3 calcNormal( in vec3 & p ) // for function f(p)
{
    const float h = 0.0001;      // replace by an appropriate value
    #define ZERO (min(iFrame,0)) // non-constant zero
    vec3 n = vec3(0.0);
    for( int i=ZERO; i<4; i++ )
    {
        vec3 e = 0.5773*(2.0*vec3((((i+3)>>1)&1),((i>>1)&1),(i&1))-1.0);
        n += e*map(pos+e*h).x;
    }
    return normalize(n);
}
*/

void old_get_normal(double* result, double* pos) {
	const double h = 0.001;
	result[0] = 0.0;
	result[1] = 0.0;
	result[2] = 0.0;

	for (int i = 0; i < 4; i++) {
		double e[VEC3_SIZE];
		int x = (((i+3)>>1)&1);
		int y = ((i>>1)&1);
		int z = (i&1);
		double tmp[VEC3_SIZE] = { (double)x, (double)y, (double)z };
		vec3_multiply_f(e, tmp, 2.0);
		vec3_subtract_f(e, e, 1.0);
		vec3_multiply_f(e, e, 0.5773);

		double tmp_pos[VEC3_SIZE] = { 0.0, 0.0, 0.0 };
		vec3_multiply_f(tmp_pos, e, h);
		vec3_add(tmp_pos, tmp_pos, pos);

		vec3_multiply_f(tmp, e, map(tmp_pos).distance);
		vec3_add(result, result, tmp);
	}

	vec3_normalize(result, result);
}
