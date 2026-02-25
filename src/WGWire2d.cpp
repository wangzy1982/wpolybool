#include "WGWire2d.h"
#include "WGCurveIntersecter2d.h"

WG_OBJECT_TYPE_IMP(WGWire2d, WGObject2d);

WGWire2d::WGWire2d() :
    m_closed(false) {
}

WGWire2d::~WGWire2d() {
    Clear();
}

void WGWire2d::Reserve(int edge_capacity) {
    m_edges.reserve(edge_capacity);
}

int WGWire2d::AddEdge(WGEdge2d* edge) {
    int index = (int)m_edges.size();
    m_edges.push_back(edge);
    m_closed = false;
    return index;
}

void WGWire2d::RemoveEdge(int index) {
    WGEdge2d* edge = m_edges[index];
    m_edges.erase(m_edges.begin() + index);
    delete edge;
    m_closed = false;
}

void WGWire2d::Clear() {
    for (auto edge_itr = m_edges.begin(); edge_itr != m_edges.end(); ++edge_itr) {
        delete *edge_itr;
    }
    m_edges.clear();
    m_closed = false;
}

void WGWire2d::SetClosed(bool closed) {
    m_closed = closed;
}

bool WGWire2d::GetClosed() const {
    return m_closed;
}

int WGWire2d::GetEdgeCount() const {
    return (int)m_edges.size();
}

const WGEdge2d* WGWire2d::GetEdge(int index) const {
    return m_edges[index];
}

bool WGWire2d::CheckSelfIntersection() const {
    //todo
    return false;
}

void WGWire2d::Reverse(int start_index) {
    int n = (int)m_edges.size();
    if (n > start_index) {
        int m = (n + start_index) / 2;
        for (int i = start_index; i < m; ++i) {
            int j = n - i - 1 + start_index;
            WGEdge2d* edge1 = m_edges.at(i);
            WGEdge2d* edge2 = m_edges.at(j);
            m_edges.at(i) = edge2;
            m_edges.at(j) = edge1;
            edge1->Reverse();
            edge2->Reverse();
        }
        if ((n - start_index) & 1) {
            WGEdge2d* edge = m_edges.at(m);
            edge->Reverse();
        }
    }
}

void WGWire2d::Reverse() {
    Reverse(0);
}

WGBox2d WGWire2d::CalculateBox() const {
    WGBox2d box;
    for (auto edge_itr = m_edges.begin(); edge_itr != m_edges.end(); ++edge_itr) {
        box.Merge((*edge_itr)->CalculateBox());
    }
    return box;
}

WGWire2d* WGWire2d::Clone() const {
    WGWire2d* wire = new WGWire2d();
    wire->m_edges.reserve(m_edges.size());
    for (auto edge_itr = m_edges.begin(); edge_itr != m_edges.end(); ++edge_itr) {
        wire->m_edges.push_back((*edge_itr)->Clone());
    }
    wire->m_closed = m_closed;
    return wire;
}

bool WGWire2d::RepairSingularity(double distance_epsilon) {
    std::vector<WGCurve2d*> curves;
    for (int i = 0; i < (int)m_edges.size(); ++i) {
        WGEdge2d* edge = m_edges.at(i);
        const WGCurve2d* curve = edge->GetCurve();
        WGCurveIntersecter2d::RemoveSingularity(curve, distance_epsilon, curves);
        if (curves.size() != 1 || curves.at(0) != curve) {
            std::vector<WGEdge2d*> old_edges = std::move(m_edges);
            m_edges.reserve(old_edges.size());
            if (i > 0) {
                m_edges.insert(m_edges.end(), old_edges.begin(), old_edges.begin() + i);
            }
            delete edge;
            for (auto curve_itr = curves.begin(); curve_itr != curves.end(); ++curve_itr) {
                m_edges.push_back(new WGEdge2d(*curve_itr));
            }
            for (int j = i + 1; j < (int)old_edges.size(); ++j) {
                edge = old_edges.at(j);
                curve = edge->GetCurve();
                WGCurveIntersecter2d::RemoveSingularity(curve, distance_epsilon, curves);
                if (curves.size() != 1 || curves.at(0) != curve) {
                    delete edge;
                    for (auto curve_itr = curves.begin(); curve_itr != curves.end(); ++curve_itr) {
                        m_edges.push_back(new WGEdge2d(*curve_itr));
                    }
                }
                else {
                    m_edges.push_back(edge);
                }
            }
            return true;
        }
    }
    return false;
}

void WGWire2d::Move(const WGVector2d& vt) {
    for (auto edge_itr = m_edges.begin(); edge_itr != m_edges.end(); ++edge_itr) {
        WGEdge2d* edge = *edge_itr;
        edge->Move(vt);
    }
}