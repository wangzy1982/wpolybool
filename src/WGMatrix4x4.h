#ifndef _WG_MATRIX4X4_
#define _WG_MATRIX4X4_

#include "WGVector3d.h"
#include "WGQuaternion.h"

enum class WGMatrix4x4Type {
    Undefine,
    Identity,
    TRS,
    MultTRS,
    Unknown
};

class WGMatrix4x4 {
public:
    WGMatrix4x4();
    double GetElement(int row, int col) const;
    void SetElement(int row, int col, double d);
    WGMatrix4x4Type GetType() const;
    void SetType(WGMatrix4x4Type type);
    WGMatrix4x4 MulMatrix(const WGMatrix4x4& matrix) const;
    WGVector3d MulPoint(const WGVector3d& point) const;
    WGVector3d MulVector(const WGVector3d& vt) const;
    bool Inverse(WGMatrix4x4& inverse_matrix) const;
    WGMatrix4x4 MulTranslateMatrix(const WGVector3d& translate) const;
    WGMatrix4x4 MulScaleMatrix(const WGVector3d& scale) const;
private:
    WGMatrix4x4Type m_type;
    double m_elements[4][4];
public:
    static WGMatrix4x4 BuildIdentity();
    static WGMatrix4x4 BuildTranslation(const WGVector3d& t);
    static WGMatrix4x4 BuildRotation(const WGQuaternion& r);
    static WGMatrix4x4 BuildScale(const WGVector3d& s);
    static WGMatrix4x4 BuildTRS(const WGVector3d& t, const WGQuaternion& r, const WGVector3d& s);
    static WGMatrix4x4 BuildInverseTRS(const WGVector3d& t, const WGQuaternion& r, const WGVector3d& s);
    static WGMatrix4x4 BuildOrtho(double l, double r, double b, double t, double n, double f);
    static WGMatrix4x4 BuildFrustum(double l, double r, double b, double t, double n, double f);
public:
    static WGMatrix4x4 Identity;
};

#endif
