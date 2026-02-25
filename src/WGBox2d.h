#ifndef _WG_BOX_2D_
#define _WG_BOX_2D_

#include "WGVector2d.h"

class WGBox2d {
public:
    WGBox2d();
    WGBox2d(const WGVector2d& point);
    WGBox2d(const WGVector2d& min, const WGVector2d& max);
    const WGVector2d& GetMin() const;
    const WGVector2d& GetMax() const;
    void Merge(const WGVector2d& point);
    void Merge(const WGBox2d& other);
    bool IsIntersected(const WGVector2d& other, double epsilon) const;
    bool IsIntersected(const WGBox2d& other, double epsilon) const;
    void Move(const WGVector2d& vt);
private:
    WGVector2d m_min;
    WGVector2d m_max;
};

#endif
