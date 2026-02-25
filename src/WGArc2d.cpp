#include "WGArc2d.h"
#include <assert.h>
#include "wsbasis.h"

WG_OBJECT_TYPE_IMP(WGArc2d, WGCurve2d);

WGArc2d::WGArc2d() :
    m_center(0, 0),
    m_radius(1),
    m_start_angle(0),
    m_delta_angle(g_pi * 2) {
}

WGArc2d::WGArc2d(const WGVector2d& center, double radius, double start_angle, double delta_angle) :
    m_center(center),
    m_radius(radius),
    m_start_angle(start_angle),
    m_delta_angle(delta_angle) {
    assert(m_radius > 0);
    assert(m_delta_angle != 0);
}

WGVector2d WGArc2d::GetStartPoint() const {
    return CalculatePoint(m_center, m_radius, m_start_angle);
}

WGVector2d WGArc2d::GetEndPoint() const {
    return CalculatePoint(m_center, m_radius, m_start_angle + m_delta_angle);
}

const WGVector2d& WGArc2d::GetCenter() const {
    return m_center;
}

double WGArc2d::GetRadius() const {
    return m_radius;
}

double WGArc2d::GetStartAngle() const {
    return m_start_angle;
}

double WGArc2d::GetDeltaAngle() const {
    return m_delta_angle;
}

void WGArc2d::SetCenter(const WGVector2d& center) {
    m_center = center;
}

void WGArc2d::SetRadius(double radius) {
    assert(radius > 0);
    m_radius = radius;
}

void WGArc2d::SetStartAngle(double start_angle) {
    assert(start_angle >= 0 && start_angle < g_pi * 2);
    m_start_angle = start_angle;
}

void WGArc2d::SetDeltaAngle(double delta_angle) {
    assert(delta_angle != 0);
    m_delta_angle = delta_angle;
}

double WGArc2d::GetT(double angle) const {
    return CalculateT(m_start_angle, m_delta_angle, angle);
}

double WGArc2d::GetAngle(double t) const {
    double a = m_start_angle + m_delta_angle * t;
    if (a < 0) {
        a += g_pi * 2;
    }
    else if (a >= g_pi * 2) {
        a -= g_pi * 2;
    }
    return a;
}

int WGArc2d::GetPieceCount() const {
    return 1;
}

WSInterval WGArc2d::GetPieceDomain(int index) const {
    return WSInterval(0, 1);
}

WGVector2d WGArc2d::CalculateD0(int piece_index, double t) const {
    double angle = m_start_angle + t * m_delta_angle;
    return CalculatePoint(m_center, m_radius, angle);
}

WGVector2d WGArc2d::CalculateD1(int piece_index, double t) const {
    double angle = m_start_angle + t * m_delta_angle;
    return WGVector2d(-m_radius * sin(angle) * m_delta_angle, m_radius * cos(angle) * m_delta_angle);
}

double WGArc2d::CalculateVariableEpsilon(int piece_index, double distance_epsilon) const {
    double length = abs(m_delta_angle) * m_radius;
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

WGBox2d WGArc2d::CalculateBox() const {
    double start_angle = m_start_angle;
    double delta_angle = m_delta_angle;
    if (delta_angle < 0) {
        start_angle += delta_angle;
        delta_angle = -delta_angle;
    }
    if (delta_angle >= g_pi * 2) {
        return WGBox2d(
            WGVector2d(m_center.X - m_radius, m_center.Y - m_radius),
            WGVector2d(m_center.X + m_radius, m_center.Y + m_radius)
        );
    }
    start_angle -= (g_pi * 2) * (int)(start_angle / (g_pi * 2));
    WGBox2d box(CalculatePoint(m_center, m_radius, start_angle));
    box.Merge(CalculatePoint(m_center, m_radius, start_angle + delta_angle));
    if (start_angle < 0) {
        start_angle += g_pi * 2;
    }
    assert(start_angle >= 0 && start_angle < g_pi * 2);
    double end_angle = start_angle + delta_angle;
    if (start_angle < g_pi * 0.5) {
        if (end_angle > g_pi * 0.5) {
            box.Merge(WGVector2d(m_center.X, m_center.Y + m_radius));
            if (end_angle > g_pi) {
                box.Merge(WGVector2d(m_center.X - m_radius, m_center.Y));
                if (end_angle > g_pi * 1.5) {
                    box.Merge(WGVector2d(m_center.X, m_center.Y - m_radius));
                    if (end_angle > g_pi * 2) {
                        box.Merge(WGVector2d(m_center.X + m_radius, m_center.Y));
                    }
                }
            }
        }
    }
    else if (start_angle < g_pi) {
        if (end_angle > g_pi) {
            box.Merge(WGVector2d(m_center.X - m_radius, m_center.Y));
            if (end_angle > g_pi * 1.5) {
                box.Merge(WGVector2d(m_center.X, m_center.Y - m_radius));
                if (end_angle > g_pi * 2) {
                    box.Merge(WGVector2d(m_center.X + m_radius, m_center.Y));
                    if (end_angle > g_pi * 2.5) {
                        box.Merge(WGVector2d(m_center.X, m_center.Y + m_radius));
                    }
                }
            }
        }
    }
    else if (start_angle < g_pi * 1.5) {
        if (end_angle > g_pi * 1.5) {
            box.Merge(WGVector2d(m_center.X, m_center.Y - m_radius));
            if (end_angle > g_pi * 2) {
                box.Merge(WGVector2d(m_center.X + m_radius, m_center.Y));
                if (end_angle > g_pi * 2.5) {
                    box.Merge(WGVector2d(m_center.X, m_center.Y + m_radius));
                    if (end_angle > g_pi * 3) {
                        box.Merge(WGVector2d(m_center.X - m_radius, m_center.Y));
                    }
                }
            }
        }
    }
    else {
        if (end_angle > g_pi * 2) {
            box.Merge(WGVector2d(m_center.X + m_radius, m_center.Y));
            if (end_angle > g_pi * 2.5) {
                box.Merge(WGVector2d(m_center.X, m_center.Y + m_radius));
                if (end_angle > g_pi * 3) {
                    box.Merge(WGVector2d(m_center.X - m_radius, m_center.Y));
                    if (end_angle > g_pi * 3.5) {
                        box.Merge(WGVector2d(m_center.X, m_center.Y - m_radius));
                    }
                }
            }
        }
    }
    return box;
}

double WGArc2d::CalculateAreaBelow(int piece_index, double start_t, double end_t) const {
    double start_angle = m_start_angle + start_t * m_delta_angle;
    double end_angle = m_start_angle + end_t * m_delta_angle;
    return m_center.Y * m_radius * (cos(end_angle) - cos(start_angle)) - 
        m_radius * m_radius * ((end_angle - start_angle) * 0.5 - (sin(2 * end_angle) - sin(2 * start_angle)) * 0.25);
}

WGCurve2d* WGArc2d::Clone() const {
    return new WGArc2d(m_center, m_radius, m_start_angle, m_delta_angle);
}

WGCurve2d* WGArc2d::CreateSubCurve(int start_piece_index, double start_t, int end_piece_index, double end_t) const {
    assert(end_t > start_t);
    return new WGArc2d(m_center, m_radius, m_start_angle + start_t * m_delta_angle, (end_t - start_t) * m_delta_angle);
}

void WGArc2d::Reverse() {
    m_start_angle = m_start_angle + m_delta_angle;
    if (m_start_angle < 0) {
        m_start_angle += g_pi * 2;
    }
    else if (m_start_angle >= g_pi * 2) {
        m_start_angle -= g_pi * 2;
    }
    m_delta_angle = -m_delta_angle;
}

void WGArc2d::Linearize(std::vector<WGVector2d>& vertices, bool is_contiguous_with_prev,
    double limit_distance_epsilon, double limit_angle_epsilon, double distance_epsilon, double angle_epsilon) const {
    if (distance_epsilon < m_radius) {
        double a = 2 * acos(1 - distance_epsilon / m_radius);
        if (a > angle_epsilon) {
            angle_epsilon = a;
        }
    }
    else {
        if (g_pi > angle_epsilon) {
            angle_epsilon = g_pi;
        }
    }
    if (limit_angle_epsilon < angle_epsilon) {
        angle_epsilon = limit_angle_epsilon;
    }
    if (limit_distance_epsilon < m_radius) {
        double a = 2 * acos(1 - limit_distance_epsilon / m_radius);
        if (a < angle_epsilon) {
            angle_epsilon = a;
        }
    }
    if (angle_epsilon == 0) {
        angle_epsilon = 1E-12;
    }
    int count = (int)(abs(m_delta_angle) / angle_epsilon + 1 - g_double_epsilon);
    const int max_count = 1000000;
    if (count > max_count) {
        count = max_count;
    }
    if (count == 0) {
        return;
    }
    double step = m_delta_angle / count;
    double angle = m_start_angle;
    int capacity = vertices.size() + count + 1;
    if (vertices.capacity() < capacity) {
        vertices.reserve(capacity * 2);
    }
    if (vertices.size() == 0 || !is_contiguous_with_prev) {
        vertices.push_back(CalculatePoint(m_center, m_radius, angle));
    }
    for (int i = 0; i < count; ++i) {
        angle += step;
        vertices.push_back(CalculatePoint(m_center, m_radius, angle));
    }
}

void WGArc2d::Move(const WGVector2d& vt) {
    m_center = m_center + vt;
}

double WGArc2d::CalculateT(double start_angle, double delta_angle, double angle) {
    assert(angle >= 0 && angle < g_pi * 2);
    double end_angle = start_angle + delta_angle;
    if (delta_angle > 0) {
        double d = angle;
        if (d < start_angle) {
            d += g_pi * 2;
        }
        if (d <= end_angle) {
            return (d - start_angle) / delta_angle;
        }
        else {
            if (d <= start_angle + g_pi) {
                return (d - start_angle) / delta_angle;
            }
            else {
                return (d - g_pi * 2 - start_angle) / delta_angle;
            }
        }
    }
    else {
        assert(delta_angle < 0);
        double d = angle;
        if (d > start_angle) {
            d -= g_pi * 2;
        }
        if (d >= end_angle) {
            return (d - start_angle) / delta_angle;
        }
        else {
            if (d >= start_angle - g_pi) {
                return (d - start_angle) / delta_angle;
            }
            else {
                return (d + g_pi * 2 - start_angle) / delta_angle;
            }
        }
    }
}

WGVector2d WGArc2d::CalculatePoint(const WGVector2d& center, double radius, double angle) {
    return WGVector2d(center.X + radius * cos(angle), center.Y + radius * sin(angle));
}

bool WGArc2d::BuildByPoints(const WGVector2d& point0, const WGVector2d& point1, const WGVector2d& point2,
    WGVector2d& center, double& radius, double& start_angle, double& delta_angle) {
    WGVector2d vt1 = point1 - point0;
    WGVector2d vt2 = point2 - point1;
    WGVector2d pt1 = point0 + vt1 * 0.5;
    WGVector2d pt2 = point1 + vt2 * 0.5;
    vt1 = WGVector2d(vt1.Y, -vt1.X);
    vt2 = WGVector2d(vt2.Y, -vt2.X);
    double a = vt2.Cross(vt1);
    if (is_zero(a, g_double_epsilon)) {
        return false;
    }
    double t1 = vt2.Cross(pt2 - pt1) / a;
    center = pt1 + vt1 * t1;
    WGVector2d vt0 = point0 - center;
    radius = vt0.Length();
    vt0 = vt0 / radius;
    start_angle = vt0.UnitAngle();
    vt2 = (point2 - center) / radius;
    double end_angle = vt2.UnitAngle();
    if (end_angle < start_angle) {
        end_angle += g_pi * 2;
    }
    vt1 = (point1 - center) / radius;
    double inner_angle = vt1.UnitAngle();
    if (inner_angle < start_angle) {
        inner_angle += g_pi * 2;
    }
    if (inner_angle < end_angle) {
        delta_angle = end_angle - start_angle;
    }
    else {
        delta_angle = end_angle - start_angle - g_pi * 2;
    }
    return true;
}

WGArc2d* WGArc2d::BuildByPoints(const WGVector2d& point0, const WGVector2d& point1, const WGVector2d& point2) {
    WGVector2d center;
    double radius;
    double start_angle;
    double delta_angle;
    if (!BuildByPoints(point0, point1, point2, center, radius, start_angle, delta_angle)) {
        return nullptr;
    }
    return new WGArc2d(center, radius, start_angle, delta_angle);
}