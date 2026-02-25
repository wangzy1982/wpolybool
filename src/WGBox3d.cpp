#include "WGBox3d.h"
#include <cmath>

WGBox3d::WGBox3d() :
    m_min(INFINITY, INFINITY, INFINITY),
    m_max(-INFINITY, -INFINITY, -INFINITY) {
}

WGBox3d::WGBox3d(const WGVector3d& point) :
    m_min(point),
    m_max(point) {
}

WGBox3d::WGBox3d(const WGVector3d& min, const WGVector3d& max) :
    m_min(min),
    m_max(max) {
}

const WGVector3d& WGBox3d::GetMin() const {
    return m_min;
}

const WGVector3d& WGBox3d::GetMax() const {
    return m_max;
}

void WGBox3d::Merge(const WGVector3d& point) {
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
    if (point.Z < m_min.Z) {
        m_min.Z = point.Z;
    }
    if (point.Z > m_max.Z) {
        m_max.Z = point.Z;
    }
}

void WGBox3d::Merge(const WGBox3d& other) {
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
    if (other.m_min.Z < m_min.Z) {
        m_min.Z = other.m_min.Z;
    }
    if (other.m_max.Z > m_max.Z) {
        m_max.Z = other.m_max.Z;
    }
}

WGBox3d WGBox3d::MulMatrix(const WGMatrix4x4& matrix) const {
    WGBox3d box;
    box.Merge(matrix.MulPoint(m_min));
    box.Merge(matrix.MulPoint(m_max));
    box.Merge(matrix.MulPoint(WGVector3d(m_min.X, m_min.Y, m_max.Z)));
    box.Merge(matrix.MulPoint(WGVector3d(m_min.X, m_max.Y, m_max.Z)));
    box.Merge(matrix.MulPoint(WGVector3d(m_min.X, m_max.Y, m_min.Z)));
    box.Merge(matrix.MulPoint(WGVector3d(m_max.X, m_min.Y, m_max.Z)));
    box.Merge(matrix.MulPoint(WGVector3d(m_max.X, m_min.Y, m_min.Z)));
    box.Merge(matrix.MulPoint(WGVector3d(m_max.X, m_max.Y, m_min.Z)));
    return box;
}

bool WGBox3d::IsIntersected(const WGBox3d& other, double epsilon) const {
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
    if (other.m_min.Z > m_max.Z + epsilon) {
        return false;
    }
    if (other.m_max.Z < m_min.Z - epsilon) {
        return false;
    }
    return true;
}

bool WGBox3d::IsInner(const WGBox3d& other, double epsilon) const {
    return
        m_min.X >= other.m_min.X - epsilon &&
        m_min.Y >= other.m_min.Y - epsilon &&
        m_min.Z >= other.m_min.Z - epsilon &&
        m_max.X <= other.m_max.X + epsilon &&
        m_max.Y <= other.m_max.Y + epsilon &&
        m_max.Z <= other.m_max.Z + epsilon;
}