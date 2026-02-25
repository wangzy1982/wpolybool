#ifndef _WG_EDGE_2D_
#define _WG_EDGE_2D_

#include "WGCurve2d.h"

class WGEdge2d : public WGObject2d {
public:
    WG_OBJECT_TYPE_DEF(WGEdge2d, WGObject2d);
public:
    WGEdge2d(WGCurve2d* curve);
    virtual ~WGEdge2d();
    const WGCurve2d* GetCurve() const;
    void Reverse();
    WGBox2d CalculateBox() const;
    WGEdge2d* Clone() const;
public:
    virtual void Move(const WGVector2d& vt); 
private:
    WGCurve2d* m_curve;
};

#endif
