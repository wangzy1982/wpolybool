#include "WGEdge2d.h"

WG_OBJECT_TYPE_IMP(WGEdge2d, WGObject2d);

WGEdge2d::WGEdge2d(WGCurve2d* curve) :
    m_curve(curve) {
}

WGEdge2d::~WGEdge2d() {
    delete m_curve;
}

const WGCurve2d* WGEdge2d::GetCurve() const {
    return m_curve;
}

void WGEdge2d::Reverse() {
    m_curve->Reverse();
}

WGBox2d WGEdge2d::CalculateBox() const {
    return m_curve->CalculateBox();
}

WGEdge2d* WGEdge2d::Clone() const {
    return new WGEdge2d(m_curve->Clone());
}

void WGEdge2d::Move(const WGVector2d& vt) {
    m_curve->Move(vt);
}