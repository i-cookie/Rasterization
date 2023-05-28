#include "transforms.h"

#include "CGL/matrix3x3.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"

namespace CGL {

Vector2D operator*(const Matrix3x3 &m, const Vector2D &v) {
	Vector3D mv = m * Vector3D(v.x, v.y, 1);
	return Vector2D(mv.x / mv.z, mv.y / mv.z);
}

Matrix3x3 translate(float dx, float dy) {
	Matrix3x3 M_trans = Matrix3x3::identity();
	M_trans[2][0] = dx;
	M_trans[2][1] = dy;
	return M_trans;
}

Matrix3x3 scale(float sx, float sy) {
	// Part 3: Fill this in.
	Matrix3x3 M_scale = Matrix3x3::identity();
	M_scale[0][0] = sx;
	M_scale[1][1] = sy;

	return M_scale;
}

// The input argument is in degrees counterclockwise
Matrix3x3 rotate(float deg) {
	// Part 3: Fill this in.
	deg = deg * PI / 180;
	Matrix3x3 M_rotate = Matrix3x3::identity();
	M_rotate[0][0] = cos(deg);
	M_rotate[1][0] = sin(deg);
	M_rotate[0][1] = -sin(deg);
	M_rotate[1][1] = cos(deg);

	return M_rotate;
}

}
