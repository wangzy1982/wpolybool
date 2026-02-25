#ifndef _WG_BOX_3D_
#define _WG_BOX_3D_

#include "WGVector3d.h"
#include "WGMatrix4x4.h"

class WGBox3d {
public:
    WGBox3d();
    WGBox3d(const WGVector3d& point);
    WGBox3d(const WGVector3d& min, const WGVector3d& max);
    const WGVector3d& GetMin() const;
    const WGVector3d& GetMax() const;
    void Merge(const WGVector3d& point);
    void Merge(const WGBox3d& other);
    WGBox3d MulMatrix(const WGMatrix4x4& matrix) const;
    bool IsIntersected(const WGBox3d& other, double epsilon) const;
    bool IsInner(const WGBox3d& other, double epsilon) const;
private:
    WGVector3d m_min;
    WGVector3d m_max;
};

#endif
