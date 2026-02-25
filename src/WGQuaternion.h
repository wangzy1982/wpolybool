#ifndef _WG_QUATERNION_
#define _WG_QUATERNION_

#include "WGVector3d.h"

class WGQuaternion {
public:
    WGQuaternion();
    WGQuaternion(double x, double y, double z, double w);
    void Normalize();
    bool IsNormalized() const;
    WGQuaternion Inverse() const;
    WGQuaternion operator*(const WGQuaternion& other) const;
    WGVector3d operator*(const WGVector3d& vt) const;
    WGVector3d ToEuler() const;
public:
    static WGQuaternion BuildByAngleAxis(double angle, const WGVector3d& axis);
    static WGQuaternion BuildByEuler(const WGVector3d& euler_angle);
    static WGQuaternion BuildByAxis(const WGVector3d& axis_x, const WGVector3d& axis_y, const WGVector3d& axis_z);
public:
    double X;
    double Y;
    double Z;
    double W;
};

#endif
