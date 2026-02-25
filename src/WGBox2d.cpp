#include "WGBox2d.h"
#include <cmath>

WGBox2d::WGBox2d() : 
    m_min(INFINITY, INFINITY),
    m_max(-INFINITY, -INFINITY) {
}

WGBox2d::WGBox2d(const WGVector2d& point) :
    m_min(point),
    m_max(point) {
}

WGBox2d::WGBox2d(const WGVector2d& min, const WGVector2d& max) :
    m_min(min),
    m_max(max) {
}

const WGVector2d& WGBox2d::GetMin() const {
    return m_min;
}

const WGVector2d& WGBox2d::GetMax() const {
    return m_max;
}

void WGBox2d::Merge(const WGVector2d& point) {
    if (point.X < m_min.X) {
        m_min.X = point.X;
    }
    if (point.X > m_max.X) {
        m_max.X = point.X;
    }
    if (point.Y < m_min.Y) {
        m_min.Y = point.Y;
    }
    if (point.Y > m_max.Y) {
        m_max.Y = point.Y;
    }
}

void WGBox2d::Merge(const WGBox2d& other) {
    if (other.m_min.X < m_min.X) {
        m_min.X = other.m_min.X;
    }
    if (other.m_max.X > m_max.X) {
        m_max.X = other.m_max.X;
    }
    if (other.m_min.Y < m_min.Y) {
        m_min.Y = other.m_min.Y;
    }
    if (other.m_max.Y > m_max.Y) {
        m_max.Y = other.m_max.Y;
    }
}

bool WGBox2d::IsIntersected(const WGVector2d& other, double epsilon) const {
    if (other.X > m_max.X + epsilon) {
        return false;
    }
    if (other.X < m_min.X - epsilon) {
        return false;
    }
    if (other.Y > m_max.Y + epsilon) {
        return false;
    }
    if (other.Y < m_min.Y - epsilon) {
        return false;
    }
    return true;
}

bool WGBox2d::IsIntersected(const WGBox2d& other, double epsilon) const {
    if (other.m_min.X > m_max.X + epsilon) {
        return false;
    }
    if (other.m_max.X < m_min.X - epsilon) {
        return false;
    }
    if (other.m_min.Y > m_max.Y + epsilon) {
        return false;
    }
    if (other.m_max.Y < m_min.Y - epsilon) {
        return false;
    }
    return true;
}

void WGBox2d::Move(const WGVector2d& vt) {
    m_min = m_min + vt;
    m_max = m_max + vt;
}