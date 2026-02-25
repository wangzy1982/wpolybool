#include "WGLoop2d.h"
#include <assert.h>

WG_OBJECT_TYPE_IMP(WGLoop2d, WGObject2d);

WGLoop2d::WGLoop2d() :
    m_wire(nullptr),
    m_sign_area(0) {
}

WGLoop2d::~WGLoop2d() {
    delete m_wire;
}

const WGWire2d* WGLoop2d::GetWire() const {
    return m_wire;
}

void WGLoop2d::SetWire(WGWire2d* wire, double sign_area) {
    assert(!wire || wire->GetClosed());
    if (m_wire) {
        delete m_wire;
    }
    m_wire = wire;
    m_sign_area = sign_area;
}

bool WGLoop2d::IsClockwise() const {
    return m_sign_area > 0;
}

double WGLoop2d::GetArea() const {
    return abs(m_sign_area);
}

double WGLoop2d::GetSignArea() const {
    return m_sign_area;
}

WGBox2d WGLoop2d::CalculateBox() const {
    return m_wire->CalculateBox();
}

double WGLoop2d::CalculateSignArea(const WGWire2d* wire) {
    double area = 0;
    for (int i = 0, edge_count = wire->GetEdgeCount(); i < edge_count; ++i) {
        const WGCurve2d* curve = wire->GetEdge(i)->GetCurve();
        for (int j = 0, piece_count = curve->GetPieceCount(); j < piece_count; ++j) {
            WSInterval piece_domain = curve->GetPieceDomain(j);
            area += curve->CalculateAreaBelow(j, piece_domain.Min, piece_domain.Max);
        }
    }
    return area;
}

WGLoop2d* WGLoop2d::Clone() const {
    WGLoop2d* loop = new WGLoop2d();
    loop->m_wire = m_wire->Clone();
    loop->m_sign_area = m_sign_area;
    return loop;
}

bool WGLoop2d::RepairSingularity(double distance_epsilon) {
    return m_wire->RepairSingularity(distance_epsilon);
}

void WGLoop2d::Move(const WGVector2d& vt) {
    m_wire->Move(vt);
}
