#include "WGQuaternion.h"
#include <assert.h>
#include "WGUtil.h"

WGQuaternion::WGQuaternion() :
    X(0),
    Y(0),
    Z(0),
    W(1) {
}

WGQuaternion::WGQuaternion(double x, double y, double z, double w) :
    X(x),
    Y(y),
    Z(z),
    W(w) {
}

void WGQuaternion::Normalize() {
    double d = X * X + Y * Y + Z * Z + W * W;
    if (d == 0) {
        X = 0;
        Y = 0;
        Z = 0;
        W = 1;
    }
    else {
        X /= d;
        Y /= d;
        Z /= d;
        W /= d;
    }
}

bool WGQuaternion::IsNormalized() const {
    return abs(X * X + Y * Y + Z * Z + W * W - 1) <= g_double_epsilon;
}

WGQuaternion WGQuaternion::Inverse() const {
    assert(IsNormalized());
    return WGQuaternion(-X, -Y, -Z, W);
}

WGQuaternion WGQuaternion::operator*(const WGQuaternion& other) const {
    return WGQuaternion(
        Y * other.Z - other.Y * Z + X * other.W + other.X * W,
        Z * other.X - other.Z * X + Y * other.W + other.Y * W,
        X * other.Y - other.X * Y + Z * other.W + other.Z * W,
        W * other.W - X * other.X - Y * other.Y - Z * other.Z
    );
}

WGVector3d WGQuaternion::operator*(const WGVector3d& vt) const {
    WGQuaternion q = (*this) * WGQuaternion(vt.X, vt.Y, vt.Z, 0) * WGQuaternion(-X, -Y, -Z, W);
    return WGVector3d(q.X, q.Y, q.Z);
}

//Rz * Ry * Rx
WGVector3d WGQuaternion::ToEuler() const {
    double roll, pitch, yaw;
    double d = 2 * (W * Y - Z * X);
    if (d >= 1) {
        pitch = g_pi * 0.5;
        yaw = 0;
        roll = atan2(X * Y - Z * W, X * Z + W * Y);
    }
    else if (d <= -1) {
        pitch = -g_pi * 0.5;
        yaw = 0;
        roll = atan2(-X * Y + Z * W, -X * Z - W * Y);
    }
    else {
        pitch = asin(d);
        yaw = atan2(2 * (W * Z + X * Y), 1 - 2 * (Y * Y + Z * Z));
        roll = atan2(2 * (W * X + Y * Z), 1 - 2 * (X * X + Y * Y));
    }
    return WGVector3d(roll, pitch, yaw);
}

WGQuaternion WGQuaternion::BuildByAngleAxis(double angle, const WGVector3d& axis) {
    assert(axis.IsNormalized());
    double a = angle * 0.5;
    double s = sin(a);
    double c = cos(a);
    return WGQuaternion(axis.X * s, axis.Y * s, axis.Z * s, c);
}

//Rz * Ry * Rx
WGQuaternion WGQuaternion::BuildByEuler(const WGVector3d& euler_angle) {
    double r = euler_angle.X * 0.5;
    double p = euler_angle.Y * 0.5;
    double y = euler_angle.Z * 0.5;
    double cr = cos(r);
    double sr = sin(r);
    double cp = cos(p);
    double sp = sin(p);
    double cy = cos(y);
    double sy = sin(y);
    return WGQuaternion(
        sr * cp * cy - cr * sp * sy,
        cr * sp * cy + sr * cp * sy,
        cr * cp * sy - sr * sp * cy,
        cr * cp * cy + sr * sp * sy
    );
}

WGQuaternion WGQuaternion::BuildByAxis(const WGVector3d& axis_x, const WGVector3d& axis_y, const WGVector3d& axis_z) {
    assert(axis_x.IsNormalized());
    assert(axis_y.IsNormalized());
    assert(axis_z.IsNormalized());
    assert(abs(axis_x.Cross(axis_y).Dot(axis_z) - 1) <= g_unit_epsilon);
    double t = axis_x.X + axis_y.Y + axis_z.Z;
    if (t > 0) {
        double w = sqrt(1 + t) / 2;
        return WGQuaternion(
            (axis_y.Z - axis_z.Y) / (4 * w),
            (axis_z.X - axis_x.Z) / (4 * w),
            (axis_x.Y - axis_y.X) / (4 * w),
            w
        );
    }
    if (axis_x.X > axis_y.Y && axis_x.X > axis_z.Z) {
        double x = sqrt(1 + axis_x.X - axis_y.Y - axis_z.Z) / 2;
        return WGQuaternion(
            x,
            (axis_y.X + axis_x.Y) / (4 * x),
            (axis_z.X + axis_x.Z) / (4 * x),
            (axis_y.Z - axis_z.Y) / (4 * x)
        );
    }
    if (axis_y.Y > axis_z.Z) {
        double y = sqrt(1 - axis_x.X + axis_y.Y - axis_z.Z) / 2;
        return WGQuaternion(
            (axis_y.X + axis_x.Y) / (4 * y),
            y,
            (axis_z.Y + axis_y.Z) / (4 * y),
            (axis_z.X - axis_x.Z) / (4 * y)
        );
    }
    double z = sqrt(1 - axis_x.X - axis_y.Y + axis_z.Z) / 2;
    return WGQuaternion(
        (axis_z.X + axis_x.Z) / (4 * z),
        (axis_z.Y + axis_y.Z) / (4 * z),
        z,
        (axis_x.Y - axis_y.X) / (4 * z)
    );
}

