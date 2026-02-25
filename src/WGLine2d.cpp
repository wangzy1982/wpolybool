#include "WGLine2d.h"
#include <assert.h>
#include "wsbasis.h"

WG_OBJECT_TYPE_IMP(WGLine2d, WGCurve2d);

WGLine2d::WGLine2d() :
    m_start_point(0, 0),
    m_end_point(0, 0) {
}

WGLine2d::WGLine2d(const WGVector2d& start_point, const WGVector2d& end_point) :
    m_start_point(start_point),
    m_end_point(end_point) {
}

WGVector2d WGLine2d::GetStartPoint() const {
    return m_start_point;
}

WGVector2d WGLine2d::GetEndPoint() const {
    return m_end_point;
}

void WGLine2d::SetStartPoint(const WGVector2d& point) {
    m_start_point = point;
}

void WGLine2d::SetEndPoint(const WGVector2d& point) {
    m_end_point = point;
}

int WGLine2d::GetPieceCount() const {
    return 1;
}

WSInterval WGLine2d::GetPieceDomain(int index) const {
    return WSInterval(0, 1);
}

WGVector2d WGLine2d::CalculateD0(int piece_index, double t) const {
    return m_start_point * (1 - t) + m_end_point * t;
}

WGVector2d WGLine2d::CalculateD1(int piece_index, double t) const {
    return m_end_point - m_start_point;
}

double WGLine2d::CalculateVariableEpsilon(int piece_index, double distance_epsilon) const {
    double length = (m_end_point - m_start_point).Length();
    double variable_epsilon;
    if (length <= distance_epsilon) {
        variable_epsilon = 1.1;
    }
    else {
        variable_epsilon = distance_epsilon / length;
        if (variable_epsilon < g_double_epsilon) {
            variable_epsilon = g_double_epsilon;
        }
    }
    return variable_epsilon;
}

WGBox2d WGLine2d::CalculateBox() const {
    WGBox2d box(m_start_point);
    box.Merge(m_end_point);
    return box;
}

double WGLine2d::CalculateAreaBelow(int piece_index, double start_t, double end_t) const {
    return CalculateAreaBelow(CalculateD0(piece_index, start_t), CalculateD0(piece_index, end_t));
}

WGCurve2d* WGLine2d::Clone() const {
    return new WGLine2d(m_start_point, m_end_point);
}

WGCurve2d* WGLine2d::CreateSubCurve(int start_piece_index, double start_t, int end_piece_index, double end_t) const {
    assert(end_t > start_t);
    return new WGLine2d(CalculateD0(start_piece_index, start_t), CalculateD0(end_piece_index, end_t));
}

void WGLine2d::Reverse() {
    WGVector2d temp = m_start_point;
    m_start_point = m_end_point;
    m_end_point = temp;
}

void WGLine2d::Linearize(std::vector<WGVector2d>& vertices, bool is_contiguous_with_prev,
    double limit_distance_epsilon, double limit_angle_epsilon, double distance_epsilon, double angle_epsilon) const {
    int capacity = vertices.size() + 2;
    if (vertices.capacity() < capacity) {
        vertices.reserve(capacity * 2);
    }
    if (vertices.size() == 0 || !is_contiguous_with_prev) {
        vertices.push_back(m_start_point);
    }
    vertices.push_back(m_end_point);
}

double WGLine2d::CalculateAreaBelow(const WGVector2d& start_point, const WGVector2d& end_point) {
    return 0.5 * (end_point.X - start_point.X) * (end_point.Y + start_point.Y);
}

void WGLine2d::Move(const WGVector2d& vt) {
    m_start_point = m_start_point + vt;
    m_end_point = m_end_point + vt;
}