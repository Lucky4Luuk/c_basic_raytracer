/*
Copyright © 2018 Felipe Ferreira da Silva

This software is provided 'as-is', without any express or implied warranty. In
no event will the authors be held liable for any damages arising from the use of
this software.

Permission is granted to anyone to use this software for any purpose, including
commercial applications, and to alter it and redistribute it freely, subject to
the following restrictions:

  1. The origin of this software must not be misrepresented; you must not claim
     that you wrote the original software. If you use this software in a
     product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
*/

#ifndef MATHC_H
#define MATHC_H

#include <stdbool.h>
#include <math.h>

#define MATHC_VERSION_YYYY 2019
#define MATHC_VERSION_MM 02
#define MATHC_VERSION_DD 16
#define MATHC_VERSION_MICRO 0

#if !defined(MATHC_NO_INT)
#define MATHC_USE_INT
#endif
#if !defined(MATHC_NO_FLOATING_POINT)
#define MATHC_USE_FLOATING_POINT
#endif
#if !defined(MATHC_NO_POINTER_STRUCT_FUNCTIONS)
#define MATHC_USE_POINTER_STRUCT_FUNCTIONS
#endif
#if !defined(MATHC_NO_STRUCT_FUNCTIONS)
#define MATHC_USE_STRUCT_FUNCTIONS
#endif
#if !defined(MATHC_NO_EASING_FUNCTIONS)
#define MATHC_USE_EASING_FUNCTIONS
#endif

#if defined(MATHC_USE_INT)
#include <stdint.h>
#endif
#if defined(MATHC_USE_FLOATING_POINT)
#include <float.h>
#endif

#define VEC2_SIZE 2
#define VEC3_SIZE 3
#define VEC4_SIZE 4
#define QUAT_SIZE 4
#define MAT2_SIZE 4
#define MAT3_SIZE 9
#define MAT4_SIZE 16

#if defined(MATHC_USE_INT)
#if defined(MATHC_INT_TYPE)
typedef MATHC_INT_TYPE mint_t;
#endif
#if !defined(MATHC_USE_INT8) && !defined(MATHC_USE_INT16) && !defined(MATHC_USE_INT32) && !defined(MATHC_USE_INT64)
#define MATHC_USE_INT32
#endif
#if defined(MATHC_USE_INT8)
#if !defined(MATHC_INT_TYPE)
typedef int8_t mint_t;
#endif
#define MINT_MAX INT8_MAX
#define MINT_MIN INT8_MIN
#endif
#if defined(MATHC_USE_INT16)
#if !defined(MATHC_INT_TYPE)
typedef int16_t mint_t;
#endif
#define MINT_MAX INT16_MAX
#define MINT_MIN INT16_MIN
#endif
#if defined(MATHC_USE_INT32)
#if !defined(MATHC_INT_TYPE)
typedef int32_t mint_t;
#endif
#define MINT_MAX INT32_MAX
#define MINT_MIN INT32_MIN
#endif
#if defined(MATHC_USE_INT64)
#if !defined(MATHC_INT_TYPE)
typedef int64_t mint_t;
#endif
#define MINT_MAX INT64_MAX
#define MINT_MIN INT64_MIN
#endif
#endif

#if defined(MATHC_USE_FLOATING_POINT)
#if !defined(MATHC_USE_SINGLE_FLOATING_POINT) && !defined(MATHC_USE_DOUBLE_FLOATING_POINT)
#define MATHC_USE_SINGLE_FLOATING_POINT
#endif
#if defined(MATHC_USE_SINGLE_FLOATING_POINT)
#define MPI 3.1415926536f
#define MPI_2 1.5707963268f
#define MPI_4 0.7853981634f
#define MFLT_EPSILON FLT_EPSILON
#define MFABS fabsf
#define MFMIN fminf
#define MFMAX fmaxf
#define MSQRT sqrtf
#define MSIN sinf
#define MCOS cosf
#define MACOS acosf
#define MASIN asinf
#define MTAN tanf
#define MATAN2 atan2f
#define MPOW powf
#define MFLOOR floorf
#define MCEIL ceilf
#define MROUND roundf
#define MFLOAT_C(c) c ## f
#endif
#if defined(MATHC_USE_DOUBLE_FLOATING_POINT)
#define MPI 3.14159265358979323846
#define MPI_2 1.57079632679489661923
#define MPI_4 0.78539816339744830962
#define MFLT_EPSILON DBL_EPSILON
#define MFABS fabs
#define MFMIN fmin
#define MFMAX fmax
#define MSQRT sqrt
#define MSIN sin
#define MCOS cos
#define MACOS acos
#define MASIN asin
#define MTAN tan
#define MATAN2 atan2
#define MPOW pow
#define MFLOOR floor
#define MCEIL ceil
#define MROUND round
#define MFLOAT_C(c) c
#endif
#endif

#if defined(MATHC_USE_STRUCT_FUNCTIONS) || defined(MATHC_USE_POINTER_STRUCT_FUNCTIONS)
#if defined(MATHC_USE_INT)
struct vec2i {
#if defined(MATHC_USE_UNIONS)
	union {
		struct {
			mint_t x;
			mint_t y;
		};
		mint_t v[VEC2_SIZE];
	};
#else
	mint_t x;
	mint_t y;
#endif
};

struct vec3i {
#if defined(MATHC_USE_UNIONS)
	union {
		struct {
			mint_t x;
			mint_t y;
			mint_t z;
		};
		mint_t v[VEC3_SIZE];
	};
#else
	mint_t x;
	mint_t y;
	mint_t z;
#endif
};

struct vec4i {
#if defined(MATHC_USE_UNIONS)
	union {
		struct {
			mint_t x;
			mint_t y;
			mint_t z;
			mint_t w;
		};
		mint_t v[VEC4_SIZE];
	};
#else
	mint_t x;
	mint_t y;
	mint_t z;
	mint_t w;
#endif
};
#endif

#if defined(MATHC_USE_FLOATING_POINT)
struct vec2 {
#if defined(MATHC_USE_UNIONS)
	union {
		struct {
			double x;
			double y;
		};
		double v[VEC2_SIZE];
	};
#else
	double x;
	double y;
#endif
};

struct vec3 {
#if defined(MATHC_USE_UNIONS)
	union {
		struct {
			double x;
			double y;
			double z;
		};
		double v[VEC3_SIZE];
	};
#else
	double x;
	double y;
	double z;
#endif
};

struct vec4 {
#if defined(MATHC_USE_UNIONS)
	union {
		struct {
			double x;
			double y;
			double z;
			double w;
		};
		double v[VEC4_SIZE];
	};
#else
	double x;
	double y;
	double z;
	double w;
#endif
};

struct quat {
#if defined(MATHC_USE_UNIONS)
	union {
		struct {
			double x;
			double y;
			double z;
			double w;
		};
		double v[QUAT_SIZE];
	};
#else
	double x;
	double y;
	double z;
	double w;
#endif
};

/*
Matrix 2×2 representation:
0/m11 2/m12
1/m21 3/m22
*/
struct mat2 {
#if defined(MATHC_USE_UNIONS)
	union {
		struct {
			double m11;
			double m21;
			double m12;
			double m22;
		};
		double v[MAT2_SIZE];
	};
#else
	double m11;
	double m21;
	double m12;
	double m22;
#endif
};

/*
Matrix 3×3 representation:
0/m11 3/m12 6/m13
1/m21 4/m22 7/m23
2/m31 5/m32 8/m33
*/
struct mat3 {
#if defined(MATHC_USE_UNIONS)
	union {
		struct {
			double m11;
			double m21;
			double m31;
			double m12;
			double m22;
			double m32;
			double m13;
			double m23;
			double m33;
		};
		double v[MAT3_SIZE];
	};
#else
	double m11;
	double m21;
	double m31;
	double m12;
	double m22;
	double m32;
	double m13;
	double m23;
	double m33;
#endif
};

/*
Matrix 4×4 representation:
0/m11 4/m12  8/m13 12/m14
1/m21 5/m22  9/m23 13/m24
2/m31 6/m32 10/m33 14/m34
3/m41 7/m42 11/m43 15/m44
*/
struct mat4 {
#if defined(MATHC_USE_UNIONS)
	union {
		struct {
			double m11;
			double m21;
			double m31;
			double m41;
			double m12;
			double m22;
			double m32;
			double m42;
			double m13;
			double m23;
			double m33;
			double m43;
			double m14;
			double m24;
			double m34;
			double m44;
		};
		double v[MAT4_SIZE];
	};
#else
	double m11;
	double m21;
	double m31;
	double m41;
	double m12;
	double m22;
	double m32;
	double m42;
	double m13;
	double m23;
	double m33;
	double m43;
	double m14;
	double m24;
	double m34;
	double m44;
#endif
};
#endif
#endif

#if defined(MATHC_USE_INT)
mint_t clampi(mint_t value, mint_t min, mint_t max);
#endif

#if defined(MATHC_USE_FLOATING_POINT)
#define MRADIANS(degrees) (degrees * MPI / MFLOAT_C(180.0))
#define MDEGREES(radians) (radians * MFLOAT_C(180.0) / MPI)
bool nearly_equal(double a, double b, double epsilon);
double to_radians(double degrees);
double to_degrees(double radians);
double clampf(double value, double min, double max);
#endif

#if defined(MATHC_USE_INT)
bool vec2i_is_zero(mint_t *v0);
bool vec2i_is_equal(mint_t *v0, mint_t *v1);
mint_t *vec2i(mint_t *result, mint_t x, mint_t y);
mint_t *vec2i_assign(mint_t *result, mint_t *v0);
#if defined(MATHC_USE_FLOATING_POINT)
mint_t *vec2i_assign_vec2(mint_t *result, double *v0);
#endif
mint_t *vec2i_zero(mint_t *result);
mint_t *vec2i_one(mint_t *result);
mint_t *vec2i_sign(mint_t *result, mint_t *v0);
mint_t *vec2i_add(mint_t *result, mint_t *v0, mint_t *v1);
mint_t *vec2i_add_i(mint_t *result, mint_t *v0, mint_t i);
mint_t *vec2i_subtract(mint_t *result, mint_t *v0, mint_t *v1);
mint_t *vec2i_subtract_i(mint_t *result, mint_t *v0, mint_t i);
mint_t *vec2i_multiply(mint_t *result, mint_t *v0, mint_t *v1);
mint_t *vec2i_multiply_i(mint_t *result, mint_t *v0, mint_t i);
mint_t *vec2i_divide(mint_t *result, mint_t *v0, mint_t *v1);
mint_t *vec2i_divide_i(mint_t *result, mint_t *v0, mint_t i);
mint_t *vec2i_snap(mint_t *result, mint_t *v0, mint_t *v1);
mint_t *vec2i_snap_i(mint_t *result, mint_t *v0, mint_t i);
mint_t *vec2i_negative(mint_t *result, mint_t *v0);
mint_t *vec2i_abs(mint_t *result, mint_t *v0);
mint_t *vec2i_max(mint_t *result, mint_t *v0, mint_t *v1);
mint_t *vec2i_min(mint_t *result, mint_t *v0, mint_t *v1);
mint_t *vec2i_clamp(mint_t *result, mint_t *v0, mint_t *v1, mint_t *v2);
mint_t *vec2i_tangent(mint_t *result, mint_t *v0);
bool vec3i_is_zero(mint_t *v0);
bool vec3i_is_equal(mint_t *v0, mint_t *v1);
mint_t *vec3i(mint_t *result, mint_t x, mint_t y, mint_t z);
mint_t *vec3i_assign(mint_t *result, mint_t *v0);
#if defined(MATHC_USE_FLOATING_POINT)
mint_t *vec3i_assign_vec3(mint_t *result, double *v0);
#endif
mint_t *vec3i_zero(mint_t *result);
mint_t *vec3i_one(mint_t *result);
mint_t *vec3i_sign(mint_t *result, mint_t *v0);
mint_t *vec3i_add(mint_t *result, mint_t *v0, mint_t *v1);
mint_t *vec3i_add_i(mint_t *result, mint_t *v0, mint_t i);
mint_t *vec3i_subtract(mint_t *result, mint_t *v0, mint_t *v1);
mint_t *vec3i_subtract_i(mint_t *result, mint_t *v0, mint_t i);
mint_t *vec3i_multiply(mint_t *result, mint_t *v0, mint_t *v1);
mint_t *vec3i_multiply_i(mint_t *result, mint_t *v0, mint_t i);
mint_t *vec3i_divide(mint_t *result, mint_t *v0, mint_t *v1);
mint_t *vec3i_divide_i(mint_t *result, mint_t *v0, mint_t i);
mint_t *vec3i_snap(mint_t *result, mint_t *v0, mint_t *v1);
mint_t *vec3i_snap_i(mint_t *result, mint_t *v0, mint_t i);
mint_t *vec3i_cross(mint_t *result, mint_t *v0, mint_t *v1);
mint_t *vec3i_negative(mint_t *result, mint_t *v0);
mint_t *vec3i_abs(mint_t *result, mint_t *v0);
mint_t *vec3i_max(mint_t *result, mint_t *v0, mint_t *v1);
mint_t *vec3i_min(mint_t *result, mint_t *v0, mint_t *v1);
mint_t *vec3i_clamp(mint_t *result, mint_t *v0, mint_t *v1, mint_t *v2);
bool vec4i_is_zero(mint_t *v0);
bool vec4i_is_equal(mint_t *v0, mint_t *v1);
mint_t *vec4i(mint_t *result, mint_t x, mint_t y, mint_t z, mint_t w);
mint_t *vec4i_assign(mint_t *result, mint_t *v0);
#if defined(MATHC_USE_FLOATING_POINT)
mint_t *vec4i_assign_vec4(mint_t *result, double *v0);
#endif
mint_t *vec4i_zero(mint_t *result);
mint_t *vec4i_one(mint_t *result);
mint_t *vec4i_sign(mint_t *result, mint_t *v0);
mint_t *vec4i_add(mint_t *result, mint_t *v0, mint_t *v1);
mint_t *vec4i_add_i(mint_t *result, mint_t *v0, mint_t i);
mint_t *vec4i_subtract(mint_t *result, mint_t *v0, mint_t *v1);
mint_t *vec4i_subtract_i(mint_t *result, mint_t *v0, mint_t i);
mint_t *vec4i_multiply(mint_t *result, mint_t *v0, mint_t *v1);
mint_t *vec4i_multiply_i(mint_t *result, mint_t *v0, mint_t i);
mint_t *vec4i_divide(mint_t *result, mint_t *v0, mint_t *v1);
mint_t *vec4i_divide_i(mint_t *result, mint_t *v0, mint_t i);
mint_t *vec4i_snap(mint_t *result, mint_t *v0, mint_t *v1);
mint_t *vec4i_snap_i(mint_t *result, mint_t *v0, mint_t i);
mint_t *vec4i_negative(mint_t *result, mint_t *v0);
mint_t *vec4i_abs(mint_t *result, mint_t *v0);
mint_t *vec4i_max(mint_t *result, mint_t *v0, mint_t *v1);
mint_t *vec4i_min(mint_t *result, mint_t *v0, mint_t *v1);
mint_t *vec4i_clamp(mint_t *result, mint_t *v0, mint_t *v1, mint_t *v2);
#endif
#if defined(MATHC_USE_FLOATING_POINT)
bool vec2_is_zero(double *v0);
bool vec2_is_equal(double *v0, double *v1);
double *vec2(double *result, double x, double y);
double *vec2_assign(double *result, double *v0);
#if defined(MATHC_USE_INT)
double *vec2_assign_vec2i(double *result, mint_t *v0);
#endif
double *vec2_zero(double *result);
double *vec2_one(double *result);
double *vec2_sign(double *result, double *v0);
double *vec2_add(double *result, double *v0, double *v1);
double *vec2_add_f(double *result, double *v0, double f);
double *vec2_subtract(double *result, double *v0, double *v1);
double *vec2_subtract_f(double *result, double *v0, double f);
double *vec2_multiply(double *result, double *v0, double *v1);
double *vec2_multiply_f(double *result, double *v0, double f);
double *vec2_multiply_mat2(double *result, double *v0, double *m0);
double *vec2_divide(double *result, double *v0, double *v1);
double *vec2_divide_f(double *result, double *v0, double f);
double *vec2_snap(double *result, double *v0, double *v1);
double *vec2_snap_f(double *result, double *v0, double f);
double *vec2_negative(double *result, double *v0);
double *vec2_abs(double *result, double *v0);
double *vec2_floor(double *result, double *v0);
double *vec2_ceil(double *result, double *v0);
double *vec2_round(double *result, double *v0);
double *vec2_max(double *result, double *v0, double *v1);
double *vec2_min(double *result, double *v0, double *v1);
double *vec2_clamp(double *result, double *v0, double *v1, double *v2);
double *vec2_normalize(double *result, double *v0);
double vec2_dot(double *v0, double *v1);
double *vec2_project(double *result, double *v0, double *v1);
double *vec2_slide(double *result, double *v0, double *normal);
double *vec2_reflect(double *result, double *v0, double *normal);
double *vec2_tangent(double *result, double *v0);
double *vec2_rotate(double *result, double *v0, double f);
double *vec2_lerp(double *result, double *v0, double *v1, double f);
double *vec2_bezier3(double *result, double *v0, double *v1, double *v2, double f);
double *vec2_bezier4(double *result, double *v0, double *v1, double *v2, double *v3, double f);
double vec2_angle(double *v0);
double vec2_length(double *v0);
double vec2_length_squared(double *v0);
double vec2_distance(double *v0, double *v1);
double vec2_distance_squared(double *v0, double *v1);
bool vec2_linear_independent(double *v0, double *v1);
double** vec2_orthonormalization(double result[2][2], double basis[2][2]);
bool vec3_is_zero(double *v0);
bool vec3_is_equal(double *v0, double *v1);
double *vec3(double *result, double x, double y, double z);
double *vec3_assign(double *result, double *v0);
#if defined(MATHC_USE_INT)
double *vec3_assign_vec3i(double *result, mint_t *v0);
#endif
double *vec3_zero(double *result);
double *vec3_one(double *result);
double *vec3_sign(double *result, double *v0);
double *vec3_add(double *result, double *v0, double *v1);
double *vec3_add_f(double *result, double *v0, double f);
double *vec3_subtract(double *result, double *v0, double *v1);
double *vec3_subtract_f(double *result, double *v0, double f);
double *vec3_multiply(double *result, double *v0, double *v1);
double *vec3_multiply_f(double *result, double *v0, double f);
double *vec3_multiply_mat3(double *result, double *v0, double *m0);
double *vec3_divide(double *result, double *v0, double *v1);
double *vec3_divide_f(double *result, double *v0, double f);
double *vec3_snap(double *result, double *v0, double *v1);
double *vec3_snap_f(double *result, double *v0, double f);
double *vec3_negative(double *result, double *v0);
double *vec3_abs(double *result, double *v0);
double *vec3_floor(double *result, double *v0);
double *vec3_ceil(double *result, double *v0);
double *vec3_round(double *result, double *v0);
double *vec3_max(double *result, double *v0, double *v1);
double *vec3_min(double *result, double *v0, double *v1);
double *vec3_clamp(double *result, double *v0, double *v1, double *v2);
double *vec3_cross(double *result, double *v0, double *v1);
double *vec3_normalize(double *result, double *v0);
double vec3_dot(double *v0, double *v1);
double *vec3_project(double *result, double *v0, double *v1);
double *vec3_slide(double *result, double *v0, double *normal);
double *vec3_reflect(double *result, double *v0, double *normal);
double *vec3_rotate(double *result, double *v0, double *ra, double f);
double *vec3_lerp(double *result, double *v0, double *v1, double f);
double *vec3_bezier3(double *result, double *v0, double *v1, double *v2, double f);
double *vec3_bezier4(double *result, double *v0, double *v1, double *v2, double *v3, double f);
double vec3_length(double *v0);
double vec3_length_squared(double *v0);
double vec3_distance(double *v0, double *v1);
double vec3_distance_squared(double *v0, double *v1);
bool vec3_linear_independent(double *v0, double *v1, double *v2);
double** vec3_orthonormalization(double result[3][3], double basis[3][3]);
bool vec4_is_zero(double *v0);
bool vec4_is_equal(double *v0, double *v1);
double *vec4(double *result, double x, double y, double z, double w);
double *vec4_assign(double *result, double *v0);
#if defined(MATHC_USE_INT)
double *vec4_assign_vec4i(double *result, mint_t *v0);
#endif
double *vec4_zero(double *result);
double *vec4_one(double *result);
double *vec4_sign(double *result, double *v0);
double *vec4_add(double *result, double *v0, double *v1);
double *vec4_add_f(double *result, double *v0, double f);
double *vec4_subtract(double *result, double *v0, double *v1);
double *vec4_subtract_f(double *result, double *v0, double f);
double *vec4_multiply(double *result, double *v0, double *v1);
double *vec4_multiply_f(double *result, double *v0, double f);
double *vec4_multiply_mat4(double *result, double *v0, double *m0);
double *vec4_divide(double *result, double *v0, double *v1);
double *vec4_divide_f(double *result, double *v0, double f);
double *vec4_snap(double *result, double *v0, double *v1);
double *vec4_snap_f(double *result, double *v0, double f);
double *vec4_negative(double *result, double *v0);
double *vec4_abs(double *result, double *v0);
double *vec4_floor(double *result, double *v0);
double *vec4_ceil(double *result, double *v0);
double *vec4_round(double *result, double *v0);
double *vec4_max(double *result, double *v0, double *v1);
double *vec4_min(double *result, double *v0, double *v1);
double *vec4_clamp(double *result, double *v0, double *v1, double *v2);
double *vec4_normalize(double *result, double *v0);
double *vec4_lerp(double *result, double *v0, double *v1, double f);
bool quat_is_zero(double *q0);
bool quat_is_equal(double *q0, double *q1);
double *quat(double *result, double x, double y, double z, double w);
double *quat_assign(double *result, double *q0);
double *quat_zero(double *result);
double *quat_null(double *result);
double *quat_multiply(double *result, double *q0, double *q1);
double *quat_multiply_f(double *result, double *q0, double f);
double *quat_divide(double *result, double *q0, double *q1);
double *quat_divide_f(double *result, double *q0, double f);
double *quat_negative(double *result, double *q0);
double *quat_conjugate(double *result, double *q0);
double *quat_inverse(double *result, double *q0);
double *quat_normalize(double *result, double *q0);
double quat_dot(double *q0, double *q1);
double *quat_power(double *result, double *q0, double exponent);
double *quat_from_axis_angle(double *result, double *v0, double angle);
double *quat_from_vec3(double *result, double *v0, double *v1);
double *quat_from_mat4(double *result, double *m0);
double *quat_lerp(double *result, double *q0, double *q1, double f);
double *quat_slerp(double *result, double *q0, double *q1, double f);
double quat_length(double *q0);
double quat_length_squared(double *q0);
double quat_angle(double *q0, double *q1);
double *mat2(double *result, double m11, double m12, double m21, double m22);
double *mat2_zero(double *result);
double *mat2_identity(double *result);
double mat2_determinant(double *m0);
double *mat2_assign(double *result, double *m0);
double *mat2_negative(double *result, double *m0);
double *mat2_transpose(double *result, double *m0);
double *mat2_cofactor(double *result, double *m0);
double *mat2_adjugate(double *result, double *m0);
double *mat2_multiply(double *result, double *m0, double *m1);
double *mat2_multiply_f(double *result, double *m0, double f);
double *mat2_inverse(double *result, double *m0);
double *mat2_scaling(double *result, double *v0);
double *mat2_scale(double *result, double *m0, double *v0);
double *mat2_rotation_z(double *result, double f);
double *mat2_lerp(double *result, double *m0, double *m1, double f);
double *mat3(double *result, double m11, double m12, double m13, double m21, double m22, double m23, double m31, double m32, double m33);
double *mat3_zero(double *result);
double *mat3_identity(double *result);
double mat3_determinant(double *m0);
double *mat3_assign(double *result, double *m0);
double *mat3_negative(double *result, double *m0);
double *mat3_transpose(double *result, double *m0);
double *mat3_cofactor(double *result, double *m0);
double *mat3_multiply(double *result, double *m0, double *m1);
double *mat3_multiply_f(double *result, double *m0, double f);
double *mat3_inverse(double *result, double *m0);
double *mat3_scaling(double *result, double *v0);
double *mat3_scale(double *result, double *m0, double *v0);
double *mat3_rotation_x(double *result, double f);
double *mat3_rotation_y(double *result, double f);
double *mat3_rotation_z(double *result, double f);
double *mat3_rotation_axis(double *result, double *v0, double f);
double *mat3_rotation_quat(double *result, double *q0);
double *mat3_lerp(double *result, double *m0, double *m1, double f);
double *mat4(double *result, double m11, double m12, double m13, double m14, double m21, double m22, double m23, double m24, double m31, double m32, double m33, double m34, double m41, double m42, double m43, double m44);
double *mat4_zero(double *result);
double *mat4_identity(double *result);
double mat4_determinant(double *m0);
double *mat4_assign(double *result, double *m0);
double *mat4_negative(double *result, double *m0);
double *mat4_transpose(double *result, double *m0);
double *mat4_cofactor(double *result, double *m0);
double *mat4_rotation_x(double *result, double f);
double *mat4_rotation_y(double *result, double f);
double *mat4_rotation_z(double *result, double f);
double *mat4_rotation_axis(double *result, double *v0, double f);
double *mat4_rotation_quat(double *result, double *q0);
double *mat4_translation(double *result, double *m0, double *v0);
double *mat4_translate(double *result, double *m0, double *v0);
double *mat4_scaling(double *result, double *m0, double *v0);
double *mat4_scale(double *result, double *m0, double *v0);
double *mat4_multiply(double *result, double *m0, double *m1);
double *mat4_multiply_f(double *result, double *m0, double f);
double *mat4_inverse(double *result, double *m0);
double *mat4_lerp(double *result, double *m0, double *m1, double f);
double *mat4_look_at(double *result, double *position, double *target, double *up);
double *mat4_ortho(double *result, double l, double r, double b, double t, double n, double f);
double *mat4_perspective(double *result, double fov_y, double aspect, double n, double f);
double *mat4_perspective_fov(double *result, double fov, double w, double h, double n, double f);
double *mat4_perspective_infinite(double *result, double fov_y, double aspect, double n);
#endif

#if defined(MATHC_USE_STRUCT_FUNCTIONS)
#if defined(MATHC_USE_INT)
bool svec2i_is_zero(struct vec2i v0);
bool svec2i_is_equal(struct vec2i v0, struct vec2i v1);
struct vec2i svec2i(mint_t x, mint_t y);
struct vec2i svec2i_assign(struct vec2i v0);
#if defined(MATHC_USE_FLOATING_POINT)
struct vec2i svec2i_assign_vec2(struct vec2 v0);
#endif
struct vec2i svec2i_zero(void);
struct vec2i svec2i_one(void);
struct vec2i svec2i_sign(struct vec2i v0);
struct vec2i svec2i_add(struct vec2i v0, struct vec2i v1);
struct vec2i svec2i_add_i(struct vec2i v0, mint_t i);
struct vec2i svec2i_subtract(struct vec2i v0, struct vec2i v1);
struct vec2i svec2i_subtract_i(struct vec2i v0, mint_t i);
struct vec2i svec2i_multiply(struct vec2i v0, struct vec2i v1);
struct vec2i svec2i_multiply_i(struct vec2i v0, mint_t i);
struct vec2i svec2i_divide(struct vec2i v0, struct vec2i v1);
struct vec2i svec2i_divide_i(struct vec2i v0, mint_t i);
struct vec2i svec2i_snap(struct vec2i v0, struct vec2i v1);
struct vec2i svec2i_snap_i(struct vec2i v0, mint_t i);
struct vec2i svec2i_negative(struct vec2i v0);
struct vec2i svec2i_abs(struct vec2i v0);
struct vec2i svec2i_max(struct vec2i v0, struct vec2i v1);
struct vec2i svec2i_min(struct vec2i v0, struct vec2i v1);
struct vec2i svec2i_clamp(struct vec2i v0, struct vec2i v1, struct vec2i v2);
struct vec2i svec2i_tangent(struct vec2i v0);
bool svec3i_is_zero(struct vec3i v0);
bool svec3i_is_equal(struct vec3i v0, struct vec3i v1);
struct vec3i svec3i(mint_t x, mint_t y, mint_t z);
struct vec3i svec3i_assign(struct vec3i v0);
#if defined(MATHC_USE_FLOATING_POINT)
struct vec3i svec3i_assign_vec3(struct vec3 v0);
#endif
struct vec3i svec3i_zero(void);
struct vec3i svec3i_one(void);
struct vec3i svec3i_sign(struct vec3i v0);
struct vec3i svec3i_add(struct vec3i v0, struct vec3i v1);
struct vec3i svec3i_add_i(struct vec3i v0, mint_t i);
struct vec3i svec3i_subtract(struct vec3i v0, struct vec3i v1);
struct vec3i svec3i_subtract_i(struct vec3i v0, mint_t i);
struct vec3i svec3i_multiply(struct vec3i v0, struct vec3i v1);
struct vec3i svec3i_multiply_i(struct vec3i v0, mint_t i);
struct vec3i svec3i_divide(struct vec3i v0, struct vec3i v1);
struct vec3i svec3i_divide_i(struct vec3i v0, mint_t i);
struct vec3i svec3i_snap(struct vec3i v0, struct vec3i v1);
struct vec3i svec3i_snap_i(struct vec3i v0, mint_t i);
struct vec3i svec3i_cross(struct vec3i v0, struct vec3i v1);
struct vec3i svec3i_negative(struct vec3i v0);
struct vec3i svec3i_abs(struct vec3i v0);
struct vec3i svec3i_max(struct vec3i v0, struct vec3i v1);
struct vec3i svec3i_min(struct vec3i v0, struct vec3i v1);
struct vec3i svec3i_clamp(struct vec3i v0, struct vec3i v1, struct vec3i v2);
bool svec4i_is_zero(struct vec4i v0);
bool svec4i_is_equal(struct vec4i v0, struct vec4i v1);
struct vec4i svec4i(mint_t x, mint_t y, mint_t z, mint_t w);
struct vec4i svec4i_assign(struct vec4i v0);
#if defined(MATHC_USE_FLOATING_POINT)
struct vec4i svec4i_assign_vec4(struct vec4 v0);
#endif
struct vec4i svec4i_zero(void);
struct vec4i svec4i_one(void);
struct vec4i svec4i_sign(struct vec4i v0);
struct vec4i svec4i_add(struct vec4i v0, struct vec4i v1);
struct vec4i svec4i_add_i(struct vec4i v0, mint_t i);
struct vec4i svec4i_subtract(struct vec4i v0, struct vec4i v1);
struct vec4i svec4i_subtract_i(struct vec4i v0, mint_t i);
struct vec4i svec4i_multiply(struct vec4i v0, struct vec4i v1);
struct vec4i svec4i_multiply_i(struct vec4i v0, mint_t i);
struct vec4i svec4i_divide(struct vec4i v0, struct vec4i v1);
struct vec4i svec4i_divide_i(struct vec4i v0, mint_t i);
struct vec4i svec4i_snap(struct vec4i v0, struct vec4i v1);
struct vec4i svec4i_snap_i(struct vec4i v0, mint_t i);
struct vec4i svec4i_negative(struct vec4i v0);
struct vec4i svec4i_abs(struct vec4i v0);
struct vec4i svec4i_max(struct vec4i v0, struct vec4i v1);
struct vec4i svec4i_min(struct vec4i v0, struct vec4i v1);
struct vec4i svec4i_clamp(struct vec4i v0, struct vec4i v1, struct vec4i v2);
#endif
#if defined(MATHC_USE_FLOATING_POINT)
bool svec2_is_zero(struct vec2 v0);
bool svec2_is_equal(struct vec2 v0, struct vec2 v1);
struct vec2 svec2(double x, double y);
struct vec2 svec2_assign(struct vec2 v0);
#if defined(MATHC_USE_INT)
struct vec2 svec2_assign_vec2i(struct vec2i v0);
#endif
struct vec2 svec2_zero(void);
struct vec2 svec2_one(void);
struct vec2 svec2_sign(struct vec2 v0);
struct vec2 svec2_add(struct vec2 v0, struct vec2 v1);
struct vec2 svec2_add_f(struct vec2 v0, double f);
struct vec2 svec2_subtract(struct vec2 v0, struct vec2 v1);
struct vec2 svec2_subtract_f(struct vec2 v0, double f);
struct vec2 svec2_multiply(struct vec2 v0, struct vec2 v1);
struct vec2 svec2_multiply_f(struct vec2 v0, double f);
struct vec2 svec2_multiply_mat2(struct vec2 v0, struct mat2 m0);
struct vec2 svec2_divide(struct vec2 v0, struct vec2 v1);
struct vec2 svec2_divide_f(struct vec2 v0, double f);
struct vec2 svec2_snap(struct vec2 v0, struct vec2 v1);
struct vec2 svec2_snap_f(struct vec2 v0, double f);
struct vec2 svec2_negative(struct vec2 v0);
struct vec2 svec2_abs(struct vec2 v0);
struct vec2 svec2_floor(struct vec2 v0);
struct vec2 svec2_ceil(struct vec2 v0);
struct vec2 svec2_round(struct vec2 v0);
struct vec2 svec2_max(struct vec2 v0, struct vec2 v1);
struct vec2 svec2_min(struct vec2 v0, struct vec2 v1);
struct vec2 svec2_clamp(struct vec2 v0, struct vec2 v1, struct vec2 v2);
struct vec2 svec2_normalize(struct vec2 v0);
double svec2_dot(struct vec2 v0, struct vec2 v1);
struct vec2 svec2_project(struct vec2 v0, struct vec2 v1);
struct vec2 svec2_slide(struct vec2 v0, struct vec2 normal);
struct vec2 svec2_reflect(struct vec2 v0, struct vec2 normal);
struct vec2 svec2_tangent(struct vec2 v0);
struct vec2 svec2_rotate(struct vec2 v0, double f);
struct vec2 svec2_lerp(struct vec2 v0, struct vec2 v1, double f);
struct vec2 svec2_bezier3(struct vec2 v0, struct vec2 v1, struct vec2 v2, double f);
struct vec2 svec2_bezier4(struct vec2 v0, struct vec2 v1, struct vec2 v2, struct vec2 v3, double f);
double svec2_angle(struct vec2 v0);
double svec2_length(struct vec2 v0);
double svec2_length_squared(struct vec2 v0);
double svec2_distance(struct vec2 v0, struct vec2 v1);
double svec2_distance_squared(struct vec2 v0, struct vec2 v1);
bool svec3_is_zero(struct vec3 v0);
bool svec3_is_equal(struct vec3 v0, struct vec3 v1);
struct vec3 svec3(double x, double y, double z);
struct vec3 svec3_assign(struct vec3 v0);
#if defined(MATHC_USE_INT)
struct vec3 svec3_assign_vec3i(struct vec3i v0);
#endif
struct vec3 svec3_zero(void);
struct vec3 svec3_one(void);
struct vec3 svec3_sign(struct vec3 v0);
struct vec3 svec3_add(struct vec3 v0, struct vec3 v1);
struct vec3 svec3_add_f(struct vec3 v0, double f);
struct vec3 svec3_subtract(struct vec3 v0, struct vec3 v1);
struct vec3 svec3_subtract_f(struct vec3 v0, double f);
struct vec3 svec3_multiply(struct vec3 v0, struct vec3 v1);
struct vec3 svec3_multiply_f(struct vec3 v0, double f);
struct vec3 svec3_multiply_mat3(struct vec3 v0, struct mat3 m0);
struct vec3 svec3_divide(struct vec3 v0, struct vec3 v1);
struct vec3 svec3_divide_f(struct vec3 v0, double f);
struct vec3 svec3_snap(struct vec3 v0, struct vec3 v1);
struct vec3 svec3_snap_f(struct vec3 v0, double f);
struct vec3 svec3_negative(struct vec3 v0);
struct vec3 svec3_abs(struct vec3 v0);
struct vec3 svec3_floor(struct vec3 v0);
struct vec3 svec3_ceil(struct vec3 v0);
struct vec3 svec3_round(struct vec3 v0);
struct vec3 svec3_max(struct vec3 v0, struct vec3 v1);
struct vec3 svec3_min(struct vec3 v0, struct vec3 v1);
struct vec3 svec3_clamp(struct vec3 v0, struct vec3 v1, struct vec3 v2);
struct vec3 svec3_cross(struct vec3 v0, struct vec3 v1);
struct vec3 svec3_normalize(struct vec3 v0);
double svec3_dot(struct vec3 v0, struct vec3 v1);
struct vec3 svec3_project(struct vec3 v0, struct vec3 v1);
struct vec3 svec3_slide(struct vec3 v0, struct vec3 normal);
struct vec3 svec3_reflect(struct vec3 v0, struct vec3 normal);
struct vec3 svec3_rotate(struct vec3 v0, struct vec3 ra, double f);
struct vec3 svec3_lerp(struct vec3 v0, struct vec3 v1, double f);
struct vec3 svec3_bezier3(struct vec3 v0, struct vec3 v1, struct vec3 v2, double f);
struct vec3 svec3_bezier4(struct vec3 v0, struct vec3 v1, struct vec3 v2, struct vec3 v3, double f);
double svec3_length(struct vec3 v0);
double svec3_length_squared(struct vec3 v0);
double svec3_distance(struct vec3 v0, struct vec3 v1);
double svec3_distance_squared(struct vec3 v0, struct vec3 v1);
bool svec4_is_zero(struct vec4 v0);
bool svec4_is_equal(struct vec4 v0, struct vec4 v1);
struct vec4 svec4(double x, double y, double z, double w);
struct vec4 svec4_assign(struct vec4 v0);
#if defined(MATHC_USE_INT)
struct vec4 svec4_assign_vec4i(struct vec4i v0);
#endif
struct vec4 svec4_zero(void);
struct vec4 svec4_one(void);
struct vec4 svec4_sign(struct vec4 v0);
struct vec4 svec4_add(struct vec4 v0, struct vec4 v1);
struct vec4 svec4_add_f(struct vec4 v0, double f);
struct vec4 svec4_subtract(struct vec4 v0, struct vec4 v1);
struct vec4 svec4_subtract_f(struct vec4 v0, double f);
struct vec4 svec4_multiply(struct vec4 v0, struct vec4 v1);
struct vec4 svec4_multiply_f(struct vec4 v0, double f);
struct vec4 svec4_multiply_mat4(struct vec4 v0, struct mat4 m0);
struct vec4 svec4_divide(struct vec4 v0, struct vec4 v1);
struct vec4 svec4_divide_f(struct vec4 v0, double f);
struct vec4 svec4_snap(struct vec4 v0, struct vec4 v1);
struct vec4 svec4_snap_f(struct vec4 v0, double f);
struct vec4 svec4_negative(struct vec4 v0);
struct vec4 svec4_abs(struct vec4 v0);
struct vec4 svec4_floor(struct vec4 v0);
struct vec4 svec4_ceil(struct vec4 v0);
struct vec4 svec4_round(struct vec4 v0);
struct vec4 svec4_max(struct vec4 v0, struct vec4 v1);
struct vec4 svec4_min(struct vec4 v0, struct vec4 v1);
struct vec4 svec4_clamp(struct vec4 v0, struct vec4 v1, struct vec4 v2);
struct vec4 svec4_normalize(struct vec4 v0);
struct vec4 svec4_lerp(struct vec4 v0, struct vec4 v1, double f);
bool squat_is_zero(struct quat q0);
bool squat_is_equal(struct quat q0, struct quat q1);
struct quat squat(double x, double y, double z, double w);
struct quat squat_assign(struct quat q0);
struct quat squat_zero(void);
struct quat squat_null(void);
struct quat squat_multiply(struct quat q0, struct quat q1);
struct quat squat_multiply_f(struct quat q0, double f);
struct quat squat_divide(struct quat q0, struct quat q1);
struct quat squat_divide_f(struct quat q0, double f);
struct quat squat_negative(struct quat q0);
struct quat squat_conjugate(struct quat q0);
struct quat squat_inverse(struct quat q0);
struct quat squat_normalize(struct quat q0);
double squat_dot(struct quat q0, struct quat q1);
struct quat squat_power(struct quat q0, double exponent);
struct quat squat_from_axis_angle(struct vec3 v0, double angle);
struct quat squat_from_vec3(struct vec3 v0, struct vec3 v1);
struct quat squat_from_mat4(struct mat4 m0);
struct quat squat_lerp(struct quat q0, struct quat q1, double f);
struct quat squat_slerp(struct quat q0, struct quat q1, double f);
double squat_length(struct quat q0);
double squat_length_squared(struct quat q0);
double squat_angle(struct quat q0, struct quat q1);
struct mat2 smat2(double m11, double m12, double m21, double m22);
struct mat2 smat2_zero(void);
struct mat2 smat2_identity(void);
double smat2_determinant(struct mat2 m0);
struct mat2 smat2_assign(struct mat2 m0);
struct mat2 smat2_negative(struct mat2 m0);
struct mat2 smat2_transpose(struct mat2 m0);
struct mat2 smat2_cofactor(struct mat2 m0);
struct mat2 smat2_adjugate(struct mat2 m0);
struct mat2 smat2_multiply(struct mat2 m0, struct mat2 m1);
struct mat2 smat2_multiply_f(struct mat2 m0, double f);
struct mat2 smat2_inverse(struct mat2 m0);
struct mat2 smat2_scaling(struct vec2 v0);
struct mat2 smat2_scale(struct mat2 m0, struct vec2 v0);
struct mat2 smat2_rotation_z(double f);
struct mat2 smat2_lerp(struct mat2 m0, struct mat2 m1, double f);
struct mat3 smat3(double m11, double m12, double m13, double m21, double m22, double m23, double m31, double m32, double m33);
struct mat3 smat3_zero(void);
struct mat3 smat3_identity(void);
double smat3_determinant(struct mat3 m0);
struct mat3 smat3_assign(struct mat3 m0);
struct mat3 smat3_negative(struct mat3 m0);
struct mat3 smat3_transpose(struct mat3 m0);
struct mat3 smat3_cofactor(struct mat3 m0);
struct mat3 smat3_multiply(struct mat3 m0, struct mat3 m1);
struct mat3 smat3_multiply_f(struct mat3 m0, double f);
struct mat3 smat3_inverse(struct mat3 m0);
struct mat3 smat3_scaling(struct vec3 v0);
struct mat3 smat3_scale(struct mat3 m0, struct vec3 v0);
struct mat3 smat3_rotation_x(double f);
struct mat3 smat3_rotation_y(double f);
struct mat3 smat3_rotation_z(double f);
struct mat3 smat3_rotation_axis(struct vec3 v0, double f);
struct mat3 smat3_rotation_quat(struct quat q0);
struct mat3 smat3_lerp(struct mat3 m0, struct mat3 m1, double f);
struct mat4 smat4(double m11, double m12, double m13, double m14, double m21, double m22, double m23, double m24, double m31, double m32, double m33, double m34, double m41, double m42, double m43, double m44);
struct mat4 smat4_zero(void);
struct mat4 smat4_identity(void);
double smat4_determinant(struct mat4 m0);
struct mat4 smat4_assign(struct mat4 m0);
struct mat4 smat4_negative(struct mat4 m0);
struct mat4 smat4_transpose(struct mat4 m0);
struct mat4 smat4_cofactor(struct mat4 m0);
struct mat4 smat4_rotation_x(double f);
struct mat4 smat4_rotation_y(double f);
struct mat4 smat4_rotation_z(double f);
struct mat4 smat4_rotation_axis(struct vec3 v0, double f);
struct mat4 smat4_rotation_quat(struct quat q0);
struct mat4 smat4_translation(struct mat4 m0, struct vec3 v0);
struct mat4 smat4_translate(struct mat4 m0, struct vec3 v0);
struct mat4 smat4_scaling(struct mat4 m0, struct vec3 v0);
struct mat4 smat4_scale(struct mat4 m0, struct vec3 v0);
struct mat4 smat4_multiply(struct mat4 m0, struct mat4 m1);
struct mat4 smat4_multiply_f(struct mat4 m0, double f);
struct mat4 smat4_inverse(struct mat4 m0);
struct mat4 smat4_lerp(struct mat4 m0, struct mat4 m1, double f);
struct mat4 smat4_look_at(struct vec3 position, struct vec3 target, struct vec3 up);
struct mat4 smat4_ortho(double l, double r, double b, double t, double n, double f);
struct mat4 smat4_perspective(double fov_y, double aspect, double n, double f);
struct mat4 smat4_perspective_fov(double fov, double w, double h, double n, double f);
struct mat4 smat4_perspective_infinite(double fov_y, double aspect, double n);
#endif
#endif

#if defined(MATHC_USE_POINTER_STRUCT_FUNCTIONS)
#if defined(MATHC_USE_INT)
bool psvec2i_is_zero(struct vec2i *v0);
bool psvec2i_is_equal(struct vec2i *v0, struct vec2i *v1);
struct vec2i *psvec2i(struct vec2i *result, mint_t x, mint_t y);
struct vec2i *psvec2i_assign(struct vec2i *result, struct vec2i *v0);
#if defined(MATHC_USE_FLOATING_POINT)
struct vec2i *psvec2i_assign_vec2(struct vec2i *result, struct vec2 *v0);
#endif
struct vec2i *psvec2i_zero(struct vec2i *result);
struct vec2i *psvec2i_one(struct vec2i *result);
struct vec2i *psvec2i_sign(struct vec2i *result, struct vec2i *v0);
struct vec2i *psvec2i_add(struct vec2i *result, struct vec2i *v0, struct vec2i *v1);
struct vec2i *psvec2i_add_i(struct vec2i *result, struct vec2i *v0, mint_t i);
struct vec2i *psvec2i_subtract(struct vec2i *result, struct vec2i *v0, struct vec2i *v1);
struct vec2i *psvec2i_subtract_i(struct vec2i *result, struct vec2i *v0, mint_t i);
struct vec2i *psvec2i_multiply(struct vec2i *result, struct vec2i *v0, struct vec2i *v1);
struct vec2i *psvec2i_multiply_i(struct vec2i *result, struct vec2i *v0, mint_t i);
struct vec2i *psvec2i_divide(struct vec2i *result, struct vec2i *v0, struct vec2i *v1);
struct vec2i *psvec2i_divide_i(struct vec2i *result, struct vec2i *v0, mint_t i);
struct vec2i *psvec2i_snap(struct vec2i *result, struct vec2i *v0, struct vec2i *v1);
struct vec2i *psvec2i_snap_i(struct vec2i *result, struct vec2i *v0, mint_t i);
struct vec2i *psvec2i_negative(struct vec2i *result, struct vec2i *v0);
struct vec2i *psvec2i_abs(struct vec2i *result, struct vec2i *v0);
struct vec2i *psvec2i_max(struct vec2i *result, struct vec2i *v0, struct vec2i *v1);
struct vec2i *psvec2i_min(struct vec2i *result, struct vec2i *v0, struct vec2i *v1);
struct vec2i *psvec2i_clamp(struct vec2i *result, struct vec2i *v0, struct vec2i *v1, struct vec2i *v2);
struct vec2i *psvec2i_tangent(struct vec2i *result, struct vec2i *v0);
bool psvec3i_is_zero(struct vec3i *v0);
bool psvec3i_is_equal(struct vec3i *v0, struct vec3i *v1);
struct vec3i *psvec3i(struct vec3i *result, mint_t x, mint_t y, mint_t z);
struct vec3i *psvec3i_assign(struct vec3i *result, struct vec3i *v0);
#if defined(MATHC_USE_FLOATING_POINT)
struct vec3i *psvec3i_assign_vec3(struct vec3i *result, struct vec3 *v0);
#endif
struct vec3i *psvec3i_zero(struct vec3i *result);
struct vec3i *psvec3i_one(struct vec3i *result);
struct vec3i *psvec3i_sign(struct vec3i *result, struct vec3i *v0);
struct vec3i *psvec3i_add(struct vec3i *result, struct vec3i *v0, struct vec3i *v1);
struct vec3i *psvec3i_add_i(struct vec3i *result, struct vec3i *v0, mint_t i);
struct vec3i *psvec3i_subtract(struct vec3i *result, struct vec3i *v0, struct vec3i *v1);
struct vec3i *psvec3i_subtract_i(struct vec3i *result, struct vec3i *v0, mint_t i);
struct vec3i *psvec3i_multiply(struct vec3i *result, struct vec3i *v0, struct vec3i *v1);
struct vec3i *psvec3i_multiply_i(struct vec3i *result, struct vec3i *v0, mint_t i);
struct vec3i *psvec3i_divide(struct vec3i *result, struct vec3i *v0, struct vec3i *v1);
struct vec3i *psvec3i_divide_i(struct vec3i *result, struct vec3i *v0, mint_t i);
struct vec3i *psvec3i_snap(struct vec3i *result, struct vec3i *v0, struct vec3i *v1);
struct vec3i *psvec3i_snap_i(struct vec3i *result, struct vec3i *v0, mint_t i);
struct vec3i *psvec3i_cross(struct vec3i *result, struct vec3i *v0, struct vec3i *v1);
struct vec3i *psvec3i_negative(struct vec3i *result, struct vec3i *v0);
struct vec3i *psvec3i_abs(struct vec3i *result, struct vec3i *v0);
struct vec3i *psvec3i_max(struct vec3i *result, struct vec3i *v0, struct vec3i *v1);
struct vec3i *psvec3i_min(struct vec3i *result, struct vec3i *v0, struct vec3i *v1);
struct vec3i *psvec3i_clamp(struct vec3i *result, struct vec3i *v0, struct vec3i *v1, struct vec3i *v2);
bool psvec4i_is_zero(struct vec4i *v0);
bool psvec4i_is_equal(struct vec4i *v0, struct vec4i *v1);
struct vec4i *psvec4i(struct vec4i *result, mint_t x, mint_t y, mint_t z, mint_t w);
struct vec4i *psvec4i_assign(struct vec4i *result, struct vec4i *v0);
#if defined(MATHC_USE_FLOATING_POINT)
struct vec4i *psvec4i_assign_vec4(struct vec4i *result, struct vec4 *v0);
#endif
struct vec4i *psvec4i_zero(struct vec4i *result);
struct vec4i *psvec4i_one(struct vec4i *result);
struct vec4i *psvec4i_sign(struct vec4i *result, struct vec4i *v0);
struct vec4i *psvec4i_add(struct vec4i *result, struct vec4i *v0, struct vec4i *v1);
struct vec4i *psvec4i_add_i(struct vec4i *result, struct vec4i *v0, mint_t i);
struct vec4i *psvec4i_subtract(struct vec4i *result, struct vec4i *v0, struct vec4i *v1);
struct vec4i *psvec4i_subtract_i(struct vec4i *result, struct vec4i *v0, mint_t i);
struct vec4i *psvec4i_multiply(struct vec4i *result, struct vec4i *v0, struct vec4i *v1);
struct vec4i *psvec4i_multiply_i(struct vec4i *result, struct vec4i *v0, mint_t i);
struct vec4i *psvec4i_divide(struct vec4i *result, struct vec4i *v0, struct vec4i *v1);
struct vec4i *psvec4i_divide_i(struct vec4i *result, struct vec4i *v0, mint_t i);
struct vec4i *psvec4i_snap(struct vec4i *result, struct vec4i *v0, struct vec4i *v1);
struct vec4i *psvec4i_snap_i(struct vec4i *result, struct vec4i *v0, mint_t i);
struct vec4i *psvec4i_negative(struct vec4i *result, struct vec4i *v0);
struct vec4i *psvec4i_abs(struct vec4i *result, struct vec4i *v0);
struct vec4i *psvec4i_max(struct vec4i *result, struct vec4i *v0, struct vec4i *v1);
struct vec4i *psvec4i_min(struct vec4i *result, struct vec4i *v0, struct vec4i *v1);
struct vec4i *psvec4i_clamp(struct vec4i *result, struct vec4i *v0, struct vec4i *v1, struct vec4i *v2);
#endif
#if defined(MATHC_USE_FLOATING_POINT)
bool psvec2_is_zero(struct vec2 *v0);
bool psvec2_is_equal(struct vec2 *v0, struct vec2 *v1);
struct vec2 *psvec2(struct vec2 *result, double x, double y);
struct vec2 *psvec2_assign(struct vec2 *result, struct vec2 *v0);
#if defined(MATHC_USE_INT)
struct vec2 *psvec2_assign_vec2i(struct vec2 *result, struct vec2i *v0);
#endif
struct vec2 *psvec2_zero(struct vec2 *result);
struct vec2 *psvec2_one(struct vec2 *result);
struct vec2 *psvec2_sign(struct vec2 *result, struct vec2 *v0);
struct vec2 *psvec2_add(struct vec2 *result, struct vec2 *v0, struct vec2 *v1);
struct vec2 *psvec2_add_f(struct vec2 *result, struct vec2 *v0, double f);
struct vec2 *psvec2_subtract(struct vec2 *result, struct vec2 *v0, struct vec2 *v1);
struct vec2 *psvec2_subtract_f(struct vec2 *result, struct vec2 *v0, double f);
struct vec2 *psvec2_multiply(struct vec2 *result, struct vec2 *v0, struct vec2 *v1);
struct vec2 *psvec2_multiply_f(struct vec2 *result, struct vec2 *v0, double f);
struct vec2 *psvec2_multiply_mat2(struct vec2 *result, struct vec2 *v0, struct mat2 *m0);
struct vec2 *psvec2_divide(struct vec2 *result, struct vec2 *v0, struct vec2 *v1);
struct vec2 *psvec2_divide_f(struct vec2 *result, struct vec2 *v0, double f);
struct vec2 *psvec2_snap(struct vec2 *result, struct vec2 *v0, struct vec2 *v1);
struct vec2 *psvec2_snap_f(struct vec2 *result, struct vec2 *v0, double f);
struct vec2 *psvec2_negative(struct vec2 *result, struct vec2 *v0);
struct vec2 *psvec2_abs(struct vec2 *result, struct vec2 *v0);
struct vec2 *psvec2_floor(struct vec2 *result, struct vec2 *v0);
struct vec2 *psvec2_ceil(struct vec2 *result, struct vec2 *v0);
struct vec2 *psvec2_round(struct vec2 *result, struct vec2 *v0);
struct vec2 *psvec2_max(struct vec2 *result, struct vec2 *v0, struct vec2 *v1);
struct vec2 *psvec2_min(struct vec2 *result, struct vec2 *v0, struct vec2 *v1);
struct vec2 *psvec2_clamp(struct vec2 *result, struct vec2 *v0, struct vec2 *v1, struct vec2 *v2);
struct vec2 *psvec2_normalize(struct vec2 *result, struct vec2 *v0);
double psvec2_dot(struct vec2 *v0, struct vec2 *v1);
struct vec2 *psvec2_project(struct vec2 *result, struct vec2 *v0, struct vec2 *v1);
struct vec2 *psvec2_slide(struct vec2 *result, struct vec2 *v0, struct vec2 *normal);
struct vec2 *psvec2_reflect(struct vec2 *result, struct vec2 *v0, struct vec2 *normal);
struct vec2 *psvec2_tangent(struct vec2 *result, struct vec2 *v0);
struct vec2 *psvec2_rotate(struct vec2 *result, struct vec2 *v0, double f);
struct vec2 *psvec2_lerp(struct vec2 *result, struct vec2 *v0, struct vec2 *v1, double f);
struct vec2 *psvec2_bezier3(struct vec2 *result, struct vec2 *v0, struct vec2 *v1, struct vec2 *v2, double f);
struct vec2 *psvec2_bezier4(struct vec2 *result, struct vec2 *v0, struct vec2 *v1, struct vec2 *v2, struct vec2 *v3, double f);
double psvec2_angle(struct vec2 *v0);
double psvec2_length(struct vec2 *v0);
double psvec2_length_squared(struct vec2 *v0);
double psvec2_distance(struct vec2 *v0, struct vec2 *v1);
double psvec2_distance_squared(struct vec2 *v0, struct vec2 *v1);
bool psvec3_is_zero(struct vec3 *v0);
bool psvec3_is_equal(struct vec3 *v0, struct vec3 *v1);
struct vec3 *psvec3(struct vec3 *result, double x, double y, double z);
struct vec3 *psvec3_assign(struct vec3 *result, struct vec3 *v0);
#if defined(MATHC_USE_INT)
struct vec3 *psvec3_assign_vec3i(struct vec3 *result, struct vec3i *v0);
#endif
struct vec3 *psvec3_zero(struct vec3 *result);
struct vec3 *psvec3_one(struct vec3 *result);
struct vec3 *psvec3_sign(struct vec3 *result, struct vec3 *v0);
struct vec3 *psvec3_add(struct vec3 *result, struct vec3 *v0, struct vec3 *v1);
struct vec3 *psvec3_add_f(struct vec3 *result, struct vec3 *v0, double f);
struct vec3 *psvec3_subtract(struct vec3 *result, struct vec3 *v0, struct vec3 *v1);
struct vec3 *psvec3_subtract_f(struct vec3 *result, struct vec3 *v0, double f);
struct vec3 *psvec3_multiply(struct vec3 *result, struct vec3 *v0, struct vec3 *v1);
struct vec3 *psvec3_multiply_f(struct vec3 *result, struct vec3 *v0, double f);
struct vec3 *psvec3_multiply_mat3(struct vec3 *result, struct vec3 *v0, struct mat3 *m0);
struct vec3 *psvec3_divide(struct vec3 *result, struct vec3 *v0, struct vec3 *v1);
struct vec3 *psvec3_divide_f(struct vec3 *result, struct vec3 *v0, double f);
struct vec3 *psvec3_snap(struct vec3 *result, struct vec3 *v0, struct vec3 *v1);
struct vec3 *psvec3_snap_f(struct vec3 *result, struct vec3 *v0, double f);
struct vec3 *psvec3_negative(struct vec3 *result, struct vec3 *v0);
struct vec3 *psvec3_abs(struct vec3 *result, struct vec3 *v0);
struct vec3 *psvec3_floor(struct vec3 *result, struct vec3 *v0);
struct vec3 *psvec3_ceil(struct vec3 *result, struct vec3 *v0);
struct vec3 *psvec3_round(struct vec3 *result, struct vec3 *v0);
struct vec3 *psvec3_max(struct vec3 *result, struct vec3 *v0, struct vec3 *v1);
struct vec3 *psvec3_min(struct vec3 *result, struct vec3 *v0, struct vec3 *v1);
struct vec3 *psvec3_clamp(struct vec3 *result, struct vec3 *v0, struct vec3 *v1, struct vec3 *v2);
struct vec3 *psvec3_cross(struct vec3 *result, struct vec3 *v0, struct vec3 *v1);
struct vec3 *psvec3_normalize(struct vec3 *result, struct vec3 *v0);
double psvec3_dot(struct vec3 *v0, struct vec3 *v1);
struct vec3 *psvec3_project(struct vec3 *result, struct vec3 *v0, struct vec3 *v1);
struct vec3 *psvec3_slide(struct vec3 *result, struct vec3 *v0, struct vec3 *normal);
struct vec3 *psvec3_reflect(struct vec3 *result, struct vec3 *v0, struct vec3 *normal);
struct vec3 *psvec3_rotate(struct vec3 *result, struct vec3 *v0, struct vec3 *ra, double f);
struct vec3 *psvec3_lerp(struct vec3 *result, struct vec3 *v0, struct vec3 *v1, double f);
struct vec3 *psvec3_bezier3(struct vec3 *result, struct vec3 *v0, struct vec3 *v1, struct vec3 *v2, double f);
struct vec3 *psvec3_bezier4(struct vec3 *result, struct vec3 *v0, struct vec3 *v1, struct vec3 *v2, struct vec3 *v3, double f);
double psvec3_length(struct vec3 *v0);
double psvec3_length_squared(struct vec3 *v0);
double psvec3_distance(struct vec3 *v0, struct vec3 *v1);
double psvec3_distance_squared(struct vec3 *v0, struct vec3 *v1);
bool psvec4_is_zero(struct vec4 *v0);
bool psvec4_is_equal(struct vec4 *v0, struct vec4 *v1);
struct vec4 *psvec4(struct vec4 *result, double x, double y, double z, double w);
struct vec4 *psvec4_assign(struct vec4 *result, struct vec4 *v0);
#if defined(MATHC_USE_INT)
struct vec4 *psvec4_assign_vec4i(struct vec4 *result, struct vec4i *v0);
#endif
struct vec4 *psvec4_zero(struct vec4 *result);
struct vec4 *psvec4_one(struct vec4 *result);
struct vec4 *psvec4_sign(struct vec4 *result, struct vec4 *v0);
struct vec4 *psvec4_add(struct vec4 *result, struct vec4 *v0, struct vec4 *v1);
struct vec4 *psvec4_add_f(struct vec4 *result, struct vec4 *v0, double f);
struct vec4 *psvec4_subtract(struct vec4 *result, struct vec4 *v0, struct vec4 *v1);
struct vec4 *psvec4_subtract_f(struct vec4 *result, struct vec4 *v0, double f);
struct vec4 *psvec4_multiply(struct vec4 *result, struct vec4 *v0, struct vec4 *v1);
struct vec4 *psvec4_multiply_f(struct vec4 *result, struct vec4 *v0, double f);
struct vec4 *psvec4_multiply_mat4(struct vec4 *result, struct vec4 *v0, struct mat4 *m0);
struct vec4 *psvec4_divide(struct vec4 *result, struct vec4 *v0, struct vec4 *v1);
struct vec4 *psvec4_divide_f(struct vec4 *result, struct vec4 *v0, double f);
struct vec4 *psvec4_snap(struct vec4 *result, struct vec4 *v0, struct vec4 *v1);
struct vec4 *psvec4_snap_f(struct vec4 *result, struct vec4 *v0, double f);
struct vec4 *psvec4_negative(struct vec4 *result, struct vec4 *v0);
struct vec4 *psvec4_abs(struct vec4 *result, struct vec4 *v0);
struct vec4 *psvec4_floor(struct vec4 *result, struct vec4 *v0);
struct vec4 *psvec4_ceil(struct vec4 *result, struct vec4 *v0);
struct vec4 *psvec4_round(struct vec4 *result, struct vec4 *v0);
struct vec4 *psvec4_max(struct vec4 *result, struct vec4 *v0, struct vec4 *v1);
struct vec4 *psvec4_min(struct vec4 *result, struct vec4 *v0, struct vec4 *v1);
struct vec4 *psvec4_clamp(struct vec4 *result, struct vec4 *v0, struct vec4 *v1, struct vec4 *v2);
struct vec4 *psvec4_normalize(struct vec4 *result, struct vec4 *v0);
struct vec4 *psvec4_lerp(struct vec4 *result, struct vec4 *v0, struct vec4 *v1, double f);
bool psquat_is_zero(struct quat *q0);
bool psquat_is_equal(struct quat *q0, struct quat *q1);
struct quat *psquat(struct quat *result, double x, double y, double z, double w);
struct quat *psquat_assign(struct quat *result, struct quat *q0);
struct quat *psquat_zero(struct quat *result);
struct quat *psquat_null(struct quat *result);
struct quat *psquat_multiply(struct quat *result, struct quat *q0, struct quat *q1);
struct quat *psquat_multiply_f(struct quat *result, struct quat *q0, double f);
struct quat *psquat_divide(struct quat *result, struct quat *q0, struct quat *q1);
struct quat *psquat_divide_f(struct quat *result, struct quat *q0, double f);
struct quat *psquat_negative(struct quat *result, struct quat *q0);
struct quat *psquat_conjugate(struct quat *result, struct quat *q0);
struct quat *psquat_inverse(struct quat *result, struct quat *q0);
struct quat *psquat_normalize(struct quat *result, struct quat *q0);
double psquat_dot(struct quat *q0, struct quat *q1);
struct quat *psquat_power(struct quat *result, struct quat *q0, double exponent);
struct quat *psquat_from_axis_angle(struct quat *result, struct vec3 *v0, double angle);
struct quat *psquat_from_vec3(struct quat *result, struct vec3 *v0, struct vec3 *v1);
struct quat *psquat_from_mat4(struct quat *result, struct mat4 *m0);
struct quat *psquat_lerp(struct quat *result, struct quat *q0, struct quat *q1, double f);
struct quat *psquat_slerp(struct quat *result, struct quat *q0, struct quat *q1, double f);
double psquat_length(struct quat *q0);
double psquat_length_squared(struct quat *q0);
double psquat_angle(struct quat *q0, struct quat *q1);
struct mat2 *psmat2(struct mat2 *result, double m11, double m12, double m21, double m22);
struct mat2 *psmat2_zero(struct mat2 *result);
struct mat2 *psmat2_identity(struct mat2 *result);
double psmat2_determinant(struct mat2 *m0);
struct mat2 *psmat2_assign(struct mat2 *result, struct mat2 *m0);
struct mat2 *psmat2_negative(struct mat2 *result, struct mat2 *m0);
struct mat2 *psmat2_transpose(struct mat2 *result, struct mat2 *m0);
struct mat2 *psmat2_cofactor(struct mat2 *result, struct mat2 *m0);
struct mat2 *psmat2_adjugate(struct mat2 *result, struct mat2 *m0);
struct mat2 *psmat2_multiply(struct mat2 *result, struct mat2 *m0, struct mat2 *m1);
struct mat2 *psmat2_multiply_f(struct mat2 *result, struct mat2 *m0, double f);
struct mat2 *psmat2_inverse(struct mat2 *result, struct mat2 *m0);
struct mat2 *psmat2_scaling(struct mat2 *result, struct vec2 *v0);
struct mat2 *psmat2_scale(struct mat2 *result, struct mat2 *m0, struct vec2 *v0);
struct mat2 *psmat2_rotation_z(struct mat2 *result, double f);
struct mat2 *psmat2_lerp(struct mat2 *result, struct mat2 *m0, struct mat2 *m1, double f);
struct mat3 *psmat3(struct mat3 *result, double m11, double m12, double m13, double m21, double m22, double m23, double m31, double m32, double m33);
struct mat3 *psmat3_zero(struct mat3 *result);
struct mat3 *psmat3_identity(struct mat3 *result);
double psmat3_determinant(struct mat3 *m0);
struct mat3 *psmat3_assign(struct mat3 *result, struct mat3 *m0);
struct mat3 *psmat3_negative(struct mat3 *result, struct mat3 *m0);
struct mat3 *psmat3_transpose(struct mat3 *result, struct mat3 *m0);
struct mat3 *psmat3_cofactor(struct mat3 *result, struct mat3 *m0);
struct mat3 *psmat3_multiply(struct mat3 *result, struct mat3 *m0, struct mat3 *m1);
struct mat3 *psmat3_multiply_f(struct mat3 *result, struct mat3 *m0, double f);
struct mat3 *psmat3_inverse(struct mat3 *result, struct mat3 *m0);
struct mat3 *psmat3_scaling(struct mat3 *result, struct vec3 *v0);
struct mat3 *psmat3_scale(struct mat3 *result, struct mat3 *m0, struct vec3 *v0);
struct mat3 *psmat3_rotation_x(struct mat3 *result, double f);
struct mat3 *psmat3_rotation_y(struct mat3 *result, double f);
struct mat3 *psmat3_rotation_z(struct mat3 *result, double f);
struct mat3 *psmat3_rotation_axis(struct mat3 *result, struct vec3 *v0, double f);
struct mat3 *psmat3_rotation_quat(struct mat3 *result, struct quat *q0);
struct mat3 *psmat3_lerp(struct mat3 *result, struct mat3 *m0, struct mat3 *m1, double f);
struct mat4 *psmat4(struct mat4 *result, double m11, double m12, double m13, double m14, double m21, double m22, double m23, double m24, double m31, double m32, double m33, double m34, double m41, double m42, double m43, double m44);
struct mat4 *psmat4_zero(struct mat4 *result);
struct mat4 *psmat4_identity(struct mat4 *result);
double psmat4_determinant(struct mat4 *m0);
struct mat4 *psmat4_assign(struct mat4 *result, struct mat4 *m0);
struct mat4 *psmat4_negative(struct mat4 *result, struct mat4 *m0);
struct mat4 *psmat4_transpose(struct mat4 *result, struct mat4 *m0);
struct mat4 *psmat4_cofactor(struct mat4 *result, struct mat4 *m0);
struct mat4 *psmat4_rotation_x(struct mat4 *result, double f);
struct mat4 *psmat4_rotation_y(struct mat4 *result, double f);
struct mat4 *psmat4_rotation_z(struct mat4 *result, double f);
struct mat4 *psmat4_rotation_axis(struct mat4 *result, struct vec3 *v0, double f);
struct mat4 *psmat4_rotation_quat(struct mat4 *result, struct quat *q0);
struct mat4 *psmat4_translation(struct mat4 *result, struct mat4 *m0, struct vec3 *v0);
struct mat4 *psmat4_translate(struct mat4 *result, struct mat4 *m0, struct vec3 *v0);
struct mat4 *psmat4_scaling(struct mat4 *result, struct mat4 *m0, struct vec3 *v0);
struct mat4 *psmat4_scale(struct mat4 *result, struct mat4 *m0, struct vec3 *v0);
struct mat4 *psmat4_multiply(struct mat4 *result, struct mat4 *m0, struct mat4 *m1);
struct mat4 *psmat4_multiply_f(struct mat4 *result, struct mat4 *m0, double f);
struct mat4 *psmat4_inverse(struct mat4 *result, struct mat4 *m0);
struct mat4 *psmat4_lerp(struct mat4 *result, struct mat4 *m0, struct mat4 *m1, double f);
struct mat4 *psmat4_look_at(struct mat4 *result, struct vec3 *position, struct vec3 *target, struct vec3 *up);
struct mat4 *psmat4_ortho(struct mat4 *result, double l, double r, double b, double t, double n, double f);
struct mat4 *psmat4_perspective(struct mat4 *result, double fov_y, double aspect, double n, double f);
struct mat4 *psmat4_perspective_fov(struct mat4 *result, double fov, double w, double h, double n, double f);
struct mat4 *psmat4_perspective_infinite(struct mat4 *result, double fov_y, double aspect, double n);
#endif
#endif

#if defined(MATHC_USE_FLOATING_POINT) && defined(MATHC_USE_EASING_FUNCTIONS)
double quadratic_ease_out(double f);
double quadratic_ease_in(double f);
double quadratic_ease_in_out(double f);
double cubic_ease_out(double f);
double cubic_ease_in(double f);
double cubic_ease_in_out(double f);
double quartic_ease_out(double f);
double quartic_ease_in(double f);
double quartic_ease_in_out(double f);
double quintic_ease_out(double f);
double quintic_ease_in(double f);
double quintic_ease_in_out(double f);
double sine_ease_out(double f);
double sine_ease_in(double f);
double sine_ease_in_out(double f);
double circular_ease_out(double f);
double circular_ease_in(double f);
double circular_ease_in_out(double f);
double exponential_ease_out(double f);
double exponential_ease_in(double f);
double exponential_ease_in_out(double f);
double elastic_ease_out(double f);
double elastic_ease_in(double f);
double elastic_ease_in_out(double f);
double back_ease_out(double f);
double back_ease_in(double f);
double back_ease_in_out(double f);
double bounce_ease_out(double f);
double bounce_ease_in(double f);
double bounce_ease_in_out(double f);
#endif

#endif
