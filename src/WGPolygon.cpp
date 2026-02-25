#include "WGPolygon.h"

WG_OBJECT_TYPE_IMP(WGPolygon, WGObject2d);

WGPolygon::WGPolygon() {
    m_box = nullptr;
    m_joint_epsilon = nullptr;
}

WGPolygon::~WGPolygon() {
    Clear();
}

void WGPolygon::Reserve(int loop_capacity) {
    m_loops.reserve(loop_capacity);
}

void WGPolygon::AddLoop(WGLoop2d* loop) {
    m_loops.push_back(loop);
    if (m_box) {
        delete m_box;
        m_box = nullptr;
    }
    if (m_joint_epsilon) {
        delete m_joint_epsilon;
        m_joint_epsilon = nullptr;
    }
}

void WGPolygon::RemoveLoop(int index) {
    delete m_loops[index];
    m_loops.erase(m_loops.begin() + index);
    if (m_box) {
        delete m_box;
        m_box = nullptr;
    }
    if (m_joint_epsilon) {
        delete m_joint_epsilon;
        m_joint_epsilon = nullptr;
    }
}

void WGPolygon::Clear() {
    for (auto loop_itr = m_loops.begin(); loop_itr != m_loops.end(); ++loop_itr) {
        delete* loop_itr;
    }
    m_loops.clear();
    if (m_box) {
        delete m_box;
        m_box = nullptr;
    }
    if (m_joint_epsilon) {
        delete m_joint_epsilon;
        m_joint_epsilon = nullptr;
    }
}

int WGPolygon::GetLoopCount() const {
    return (int)m_loops.size();
}

const WGLoop2d* WGPolygon::GetLoop(int index) const {
    return m_loops[index];
}

const WGBox2d& WGPolygon::GetBox() const {
    if (!m_box) {
        m_box = new WGBox2d();
        for (auto loop_itr = m_loops.begin(); loop_itr != m_loops.end(); ++loop_itr) {
            m_box->Merge((*loop_itr)->CalculateBox());
        }
    }
    return *m_box;
}

WGPolygon* WGPolygon::Clone() const {
    WGPolygon* polygon = new WGPolygon();
    polygon->m_loops.reserve(m_loops.size());
    for (auto loop_itr = m_loops.begin(); loop_itr != m_loops.end(); ++loop_itr) {
        polygon->m_loops.push_back((*loop_itr)->Clone());
    }
    if (m_box) {
        polygon->m_box = new WGBox2d(m_box->GetMin(), m_box->GetMax());
    }
    if (m_joint_epsilon) {
        polygon->m_joint_epsilon = new double(*m_joint_epsilon);
    }
    return polygon;
}

double WGPolygon::GetJointEpsilon() const {
    if (!m_joint_epsilon) {
        double sqr_epsilon = 0;
        for (int i = 0; i < GetLoopCount(); ++i) {
            const WGLoop2d* loop = GetLoop(i);
            const WGWire2d* wire = loop->GetWire();
            for (int j = 0; j < wire->GetEdgeCount(); ++j) {
                int k = (j + 1) % wire->GetEdgeCount();
                double d = (wire->GetEdge(j)->GetCurve()->GetEndPoint() - wire->GetEdge(k)->GetCurve()->GetStartPoint()).SqrLength();
                if (d > sqr_epsilon) {
                    sqr_epsilon = d;
                }
            }
        }
        m_joint_epsilon = new double(sqrt(sqr_epsilon));
    }
    return *m_joint_epsilon;
}

double WGPolygon::CalculateArea() const {
    double area = 0;
    for (int i = 0; i < GetLoopCount(); ++i) {
        const WGLoop2d* loop = GetLoop(i);
        area += loop->GetSignArea();
    }
    return abs(area);
}

bool WGPolygon::RepairSingularity(double distance_epsilon) {
    bool b = false;
    for (auto loop_itr = m_loops.begin(); loop_itr != m_loops.end(); ++loop_itr) {
        WGLoop2d* loop = *loop_itr;
        if (loop->RepairSingularity(distance_epsilon)) {
            b = true;
        }
    }
    if (b) {
        if (m_joint_epsilon && distance_epsilon > *m_joint_epsilon) {
            delete m_joint_epsilon;
            m_joint_epsilon = nullptr;
        }
    }
    return b;
}

void WGPolygon::Move(const WGVector2d& vt) {
    for (auto loop_itr = m_loops.begin(); loop_itr != m_loops.end(); ++loop_itr) {
        WGLoop2d* loop = *loop_itr;
        loop->Move(vt);
    }
    if (m_box) {
        m_box->Move(vt);
    }
}

WGPolygon* WGPolygon::BuildPolygon(WGCurve2d** curves, int curve_count) {
    WGPolygon* polygon = new WGPolygon();
    WGLoop2d* loop = new WGLoop2d();
    WGWire2d* wire = new WGWire2d();
    for (int i = 0; i < curve_count; ++i) {
        wire->AddEdge(new WGEdge2d(curves[i]));
    }
    wire->SetClosed(true);
    loop->SetWire(wire, WGLoop2d::CalculateSignArea(wire));
    polygon->AddLoop(loop);
    return polygon;
}