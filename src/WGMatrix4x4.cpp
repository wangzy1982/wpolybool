#include "WGMatrix4x4.h"
#include <assert.h>
#include "WGUtil.h"

WGMatrix4x4::WGMatrix4x4() :
    m_type(WGMatrix4x4Type::Undefine) {
}

double WGMatrix4x4::GetElement(int row, int col) const {
    return m_elements[row][col];
}

void WGMatrix4x4::SetElement(int row, int col, double d) {
    m_elements[row][col] = d;
    m_type = WGMatrix4x4Type::Unknown;
}

WGMatrix4x4Type WGMatrix4x4::GetType() const {
    return m_type;
}

void WGMatrix4x4::SetType(WGMatrix4x4Type type) {
    m_type = type;
}

WGMatrix4x4 WGMatrix4x4::MulMatrix(const WGMatrix4x4& matrix) const {
    assert(m_type != WGMatrix4x4Type::Undefine);
    assert(matrix.m_type != WGMatrix4x4Type::Undefine);
    if (m_type == WGMatrix4x4Type::Identity) {
        return matrix;
    }
    if (matrix.m_type == WGMatrix4x4Type::Identity) {
        return *this;
    }
    if (m_type != WGMatrix4x4Type::Unknown) {
        if (matrix.m_type != WGMatrix4x4Type::Unknown) {
            WGMatrix4x4 result;
            result.m_elements[0][0] =
                m_elements[0][0] * matrix.m_elements[0][0] +
                m_elements[0][1] * matrix.m_elements[1][0] +
                m_elements[0][2] * matrix.m_elements[2][0];
            result.m_elements[0][1] =
                m_elements[0][0] * matrix.m_elements[0][1] +
                m_elements[0][1] * matrix.m_elements[1][1] +
                m_elements[0][2] * matrix.m_elements[2][1];
            result.m_elements[0][2] =
                m_elements[0][0] * matrix.m_elements[0][2] +
                m_elements[0][1] * matrix.m_elements[1][2] +
                m_elements[0][2] * matrix.m_elements[2][2];
            result.m_elements[0][3] =
                m_elements[0][0] * matrix.m_elements[0][3] +
                m_elements[0][1] * matrix.m_elements[1][3] +
                m_elements[0][2] * matrix.m_elements[2][3] +
                m_elements[0][3];
            result.m_elements[1][0] =
                m_elements[1][0] * matrix.m_elements[0][0] +
                m_elements[1][1] * matrix.m_elements[1][0] +
                m_elements[1][2] * matrix.m_elements[2][0];
            result.m_elements[1][1] =
                m_elements[1][0] * matrix.m_elements[0][1] +
                m_elements[1][1] * matrix.m_elements[1][1] +
                m_elements[1][2] * matrix.m_elements[2][1];
            result.m_elements[1][2] =
                m_elements[1][0] * matrix.m_elements[0][2] +
                m_elements[1][1] * matrix.m_elements[1][2] +
                m_elements[1][2] * matrix.m_elements[2][2];
            result.m_elements[1][3] =
                m_elements[1][0] * matrix.m_elements[0][3] +
                m_elements[1][1] * matrix.m_elements[1][3] +
                m_elements[1][2] * matrix.m_elements[2][3] +
                m_elements[1][3];
            result.m_elements[2][0] =
                m_elements[2][0] * matrix.m_elements[0][0] +
                m_elements[2][1] * matrix.m_elements[1][0] +
                m_elements[2][2] * matrix.m_elements[2][0];
            result.m_elements[2][1] =
                m_elements[2][0] * matrix.m_elements[0][1] +
                m_elements[2][1] * matrix.m_elements[1][1] +
                m_elements[2][2] * matrix.m_elements[2][1];
            result.m_elements[2][2] =
                m_elements[2][0] * matrix.m_elements[0][2] +
                m_elements[2][1] * matrix.m_elements[1][2] +
                m_elements[2][2] * matrix.m_elements[2][2];
            result.m_elements[2][3] =
                m_elements[2][0] * matrix.m_elements[0][3] +
                m_elements[2][1] * matrix.m_elements[1][3] +
                m_elements[2][2] * matrix.m_elements[2][3] +
                m_elements[2][3];
            result.m_elements[3][0] = 0;
            result.m_elements[3][1] = 0;
            result.m_elements[3][2] = 0;
            result.m_elements[3][3] = 1;
            result.m_type = WGMatrix4x4Type::MultTRS;
            return result;
        }
        WGMatrix4x4 result;
        result.m_elements[0][0] =
            m_elements[0][0] * matrix.m_elements[0][0] +
            m_elements[0][1] * matrix.m_elements[1][0] +
            m_elements[0][2] * matrix.m_elements[2][0] +
            m_elements[0][3] * matrix.m_elements[3][0];
        result.m_elements[0][1] =
            m_elements[0][0] * matrix.m_elements[0][1] +
            m_elements[0][1] * matrix.m_elements[1][1] +
            m_elements[0][2] * matrix.m_elements[2][1] +
            m_elements[0][3] * matrix.m_elements[3][1];
        result.m_elements[0][2] =
            m_elements[0][0] * matrix.m_elements[0][2] +
            m_elements[0][1] * matrix.m_elements[1][2] +
            m_elements[0][2] * matrix.m_elements[2][2] +
            m_elements[0][3] * matrix.m_elements[3][2];
        result.m_elements[0][3] =
            m_elements[0][0] * matrix.m_elements[0][3] +
            m_elements[0][1] * matrix.m_elements[1][3] +
            m_elements[0][2] * matrix.m_elements[2][3] +
            m_elements[0][3] * matrix.m_elements[3][3];
        result.m_elements[1][0] =
            m_elements[1][0] * matrix.m_elements[0][0] +
            m_elements[1][1] * matrix.m_elements[1][0] +
            m_elements[1][2] * matrix.m_elements[2][0] +
            m_elements[1][3] * matrix.m_elements[3][0];
        result.m_elements[1][1] =
            m_elements[1][0] * matrix.m_elements[0][1] +
            m_elements[1][1] * matrix.m_elements[1][1] +
            m_elements[1][2] * matrix.m_elements[2][1] +
            m_elements[1][3] * matrix.m_elements[3][1];
        result.m_elements[1][2] =
            m_elements[1][0] * matrix.m_elements[0][2] +
            m_elements[1][1] * matrix.m_elements[1][2] +
            m_elements[1][2] * matrix.m_elements[2][2] +
            m_elements[1][3] * matrix.m_elements[3][2];
        result.m_elements[1][3] =
            m_elements[1][0] * matrix.m_elements[0][3] +
            m_elements[1][1] * matrix.m_elements[1][3] +
            m_elements[1][2] * matrix.m_elements[2][3] +
            m_elements[1][3] * matrix.m_elements[3][3];
        result.m_elements[2][0] =
            m_elements[2][0] * matrix.m_elements[0][0] +
            m_elements[2][1] * matrix.m_elements[1][0] +
            m_elements[2][2] * matrix.m_elements[2][0] +
            m_elements[2][3] * matrix.m_elements[3][0];
        result.m_elements[2][1] =
            m_elements[2][0] * matrix.m_elements[0][1] +
            m_elements[2][1] * matrix.m_elements[1][1] +
            m_elements[2][2] * matrix.m_elements[2][1] +
            m_elements[2][3] * matrix.m_elements[3][1];
        result.m_elements[2][2] =
            m_elements[2][0] * matrix.m_elements[0][2] +
            m_elements[2][1] * matrix.m_elements[1][2] +
            m_elements[2][2] * matrix.m_elements[2][2] +
            m_elements[2][3] * matrix.m_elements[3][2];
        result.m_elements[2][3] =
            m_elements[2][0] * matrix.m_elements[0][3] +
            m_elements[2][1] * matrix.m_elements[1][3] +
            m_elements[2][2] * matrix.m_elements[2][3] +
            m_elements[2][3] * matrix.m_elements[3][3];
        result.m_elements[3][0] = matrix.m_elements[3][0];
        result.m_elements[3][1] = matrix.m_elements[3][1];
        result.m_elements[3][2] = matrix.m_elements[3][2];
        result.m_elements[3][3] = matrix.m_elements[3][3];
        result.m_type = WGMatrix4x4Type::Unknown;
        return result;
    }
    if (matrix.m_type != WGMatrix4x4Type::Unknown) {
        WGMatrix4x4 result;
        result.m_elements[0][0] =
            m_elements[0][0] * matrix.m_elements[0][0] +
            m_elements[0][1] * matrix.m_elements[1][0] +
            m_elements[0][2] * matrix.m_elements[2][0];
        result.m_elements[0][1] =
            m_elements[0][0] * matrix.m_elements[0][1] +
            m_elements[0][1] * matrix.m_elements[1][1] +
            m_elements[0][2] * matrix.m_elements[2][1];
        result.m_elements[0][2] =
            m_elements[0][0] * matrix.m_elements[0][2] +
            m_elements[0][1] * matrix.m_elements[1][2] +
            m_elements[0][2] * matrix.m_elements[2][2];
        result.m_elements[0][3] =
            m_elements[0][0] * matrix.m_elements[0][3] +
            m_elements[0][1] * matrix.m_elements[1][3] +
            m_elements[0][2] * matrix.m_elements[2][3] +
            m_elements[0][3];
        result.m_elements[1][0] =
            m_elements[1][0] * matrix.m_elements[0][0] +
            m_elements[1][1] * matrix.m_elements[1][0] +
            m_elements[1][2] * matrix.m_elements[2][0];
        result.m_elements[1][1] =
            m_elements[1][0] * matrix.m_elements[0][1] +
            m_elements[1][1] * matrix.m_elements[1][1] +
            m_elements[1][2] * matrix.m_elements[2][1];
        result.m_elements[1][2] =
            m_elements[1][0] * matrix.m_elements[0][2] +
            m_elements[1][1] * matrix.m_elements[1][2] +
            m_elements[1][2] * matrix.m_elements[2][2];
        result.m_elements[1][3] =
            m_elements[1][0] * matrix.m_elements[0][3] +
            m_elements[1][1] * matrix.m_elements[1][3] +
            m_elements[1][2] * matrix.m_elements[2][3] +
            m_elements[1][3];
        result.m_elements[2][0] =
            m_elements[2][0] * matrix.m_elements[0][0] +
            m_elements[2][1] * matrix.m_elements[1][0] +
            m_elements[2][2] * matrix.m_elements[2][0];
        result.m_elements[2][1] =
            m_elements[2][0] * matrix.m_elements[0][1] +
            m_elements[2][1] * matrix.m_elements[1][1] +
            m_elements[2][2] * matrix.m_elements[2][1];
        result.m_elements[2][2] =
            m_elements[2][0] * matrix.m_elements[0][2] +
            m_elements[2][1] * matrix.m_elements[1][2] +
            m_elements[2][2] * matrix.m_elements[2][2];
        result.m_elements[2][3] =
            m_elements[2][0] * matrix.m_elements[0][3] +
            m_elements[2][1] * matrix.m_elements[1][3] +
            m_elements[2][2] * matrix.m_elements[2][3] +
            m_elements[2][3];
        result.m_elements[3][0] =
            m_elements[3][0] * matrix.m_elements[0][0] +
            m_elements[3][1] * matrix.m_elements[1][0] +
            m_elements[3][2] * matrix.m_elements[2][0];
        result.m_elements[3][1] =
            m_elements[3][0] * matrix.m_elements[0][1] +
            m_elements[3][1] * matrix.m_elements[1][1] +
            m_elements[3][2] * matrix.m_elements[2][1];
        result.m_elements[3][2] =
            m_elements[3][0] * matrix.m_elements[0][2] +
            m_elements[3][1] * matrix.m_elements[1][2] +
            m_elements[3][2] * matrix.m_elements[2][2];
        result.m_elements[3][3] =
            m_elements[3][0] * matrix.m_elements[0][3] +
            m_elements[3][1] * matrix.m_elements[1][3] +
            m_elements[3][2] * matrix.m_elements[2][3] +
            m_elements[3][3];
        result.m_type = WGMatrix4x4Type::Unknown;
        return result;
    }
    WGMatrix4x4 result;
    result.m_elements[0][0] =
        m_elements[0][0] * matrix.m_elements[0][0] +
        m_elements[0][1] * matrix.m_elements[1][0] +
        m_elements[0][2] * matrix.m_elements[2][0] +
        m_elements[0][3] * matrix.m_elements[3][0];
    result.m_elements[0][1] =
        m_elements[0][0] * matrix.m_elements[0][1] +
        m_elements[0][1] * matrix.m_elements[1][1] +
        m_elements[0][2] * matrix.m_elements[2][1] +
        m_elements[0][3] * matrix.m_elements[3][1];
    result.m_elements[0][2] =
        m_elements[0][0] * matrix.m_elements[0][2] +
        m_elements[0][1] * matrix.m_elements[1][2] +
        m_elements[0][2] * matrix.m_elements[2][2] +
        m_elements[0][3] * matrix.m_elements[3][2];
    result.m_elements[0][3] =
        m_elements[0][0] * matrix.m_elements[0][3] +
        m_elements[0][1] * matrix.m_elements[1][3] +
        m_elements[0][2] * matrix.m_elements[2][3] +
        m_elements[0][3] * matrix.m_elements[3][3];
    result.m_elements[1][0] =
        m_elements[1][0] * matrix.m_elements[0][0] +
        m_elements[1][1] * matrix.m_elements[1][0] +
        m_elements[1][2] * matrix.m_elements[2][0] +
        m_elements[1][3] * matrix.m_elements[3][0];
    result.m_elements[1][1] =
        m_elements[1][0] * matrix.m_elements[0][1] +
        m_elements[1][1] * matrix.m_elements[1][1] +
        m_elements[1][2] * matrix.m_elements[2][1] +
        m_elements[1][3] * matrix.m_elements[3][1];
    result.m_elements[1][2] =
        m_elements[1][0] * matrix.m_elements[0][2] +
        m_elements[1][1] * matrix.m_elements[1][2] +
        m_elements[1][2] * matrix.m_elements[2][2] +
        m_elements[1][3] * matrix.m_elements[3][2];
    result.m_elements[1][3] =
        m_elements[1][0] * matrix.m_elements[0][3] +
        m_elements[1][1] * matrix.m_elements[1][3] +
        m_elements[1][2] * matrix.m_elements[2][3] +
        m_elements[1][3] * matrix.m_elements[3][3];
    result.m_elements[2][0] =
        m_elements[2][0] * matrix.m_elements[0][0] +
        m_elements[2][1] * matrix.m_elements[1][0] +
        m_elements[2][2] * matrix.m_elements[2][0] +
        m_elements[2][3] * matrix.m_elements[3][0];
    result.m_elements[2][1] =
        m_elements[2][0] * matrix.m_elements[0][1] +
        m_elements[2][1] * matrix.m_elements[1][1] +
        m_elements[2][2] * matrix.m_elements[2][1] +
        m_elements[2][3] * matrix.m_elements[3][1];
    result.m_elements[2][2] =
        m_elements[2][0] * matrix.m_elements[0][2] +
        m_elements[2][1] * matrix.m_elements[1][2] +
        m_elements[2][2] * matrix.m_elements[2][2] +
        m_elements[2][3] * matrix.m_elements[3][2];
    result.m_elements[2][3] =
        m_elements[2][0] * matrix.m_elements[0][3] +
        m_elements[2][1] * matrix.m_elements[1][3] +
        m_elements[2][2] * matrix.m_elements[2][3] +
        m_elements[2][3] * matrix.m_elements[3][3];
    result.m_elements[3][0] =
        m_elements[3][0] * matrix.m_elements[0][0] +
        m_elements[3][1] * matrix.m_elements[1][0] +
        m_elements[3][2] * matrix.m_elements[2][0] +
        m_elements[3][3] * matrix.m_elements[3][0];
    result.m_elements[3][1] =
        m_elements[3][0] * matrix.m_elements[0][1] +
        m_elements[3][1] * matrix.m_elements[1][1] +
        m_elements[3][2] * matrix.m_elements[2][1] +
        m_elements[3][3] * matrix.m_elements[3][1];
    result.m_elements[3][2] =
        m_elements[3][0] * matrix.m_elements[0][2] +
        m_elements[3][1] * matrix.m_elements[1][2] +
        m_elements[3][2] * matrix.m_elements[2][2] +
        m_elements[3][3] * matrix.m_elements[3][2];
    result.m_elements[3][3] =
        m_elements[3][0] * matrix.m_elements[0][3] +
        m_elements[3][1] * matrix.m_elements[1][3] +
        m_elements[3][2] * matrix.m_elements[2][3] +
        m_elements[3][3] * matrix.m_elements[3][3];
    result.m_type = WGMatrix4x4Type::Unknown;
    return result;
}

WGVector3d WGMatrix4x4::MulPoint(const WGVector3d& point) const {
    assert(m_type != WGMatrix4x4Type::Undefine);
    if (m_type == WGMatrix4x4Type::Identity) {
        return point;
    }
    if (m_type != WGMatrix4x4Type::Unknown) {
        return WGVector3d(
            m_elements[0][0] * point.X + m_elements[0][1] * point.Y + m_elements[0][2] * point.Z + m_elements[0][3],
            m_elements[1][0] * point.X + m_elements[1][1] * point.Y + m_elements[1][2] * point.Z + m_elements[1][3],
            m_elements[2][0] * point.X + m_elements[2][1] * point.Y + m_elements[2][2] * point.Z + m_elements[2][3]
        );
    }
    double d = 1.0 / (m_elements[3][0] * point.X + m_elements[3][1] * point.Y + m_elements[3][2] * point.Z + m_elements[3][3]);
    return WGVector3d(
        (m_elements[0][0] * point.X + m_elements[0][1] * point.Y + m_elements[0][2] * point.Z + m_elements[0][3]) * d,
        (m_elements[1][0] * point.X + m_elements[1][1] * point.Y + m_elements[1][2] * point.Z + m_elements[1][3]) * d,
        (m_elements[2][0] * point.X + m_elements[2][1] * point.Y + m_elements[2][2] * point.Z + m_elements[2][3]) * d
    );
}

WGVector3d WGMatrix4x4::MulVector(const WGVector3d& vt) const {
    assert(m_type != WGMatrix4x4Type::Undefine);
    assert(m_type != WGMatrix4x4Type::Unknown);
    if (m_type == WGMatrix4x4Type::Identity) {
        return vt;
    }
    return WGVector3d(
        m_elements[0][0] * vt.X + m_elements[0][1] * vt.Y + m_elements[0][2] * vt.Z,
        m_elements[1][0] * vt.X + m_elements[1][1] * vt.Y + m_elements[1][2] * vt.Z,
        m_elements[2][0] * vt.X + m_elements[2][1] * vt.Y + m_elements[2][2] * vt.Z
    );
}

bool WGMatrix4x4::Inverse(WGMatrix4x4& inverse_matrix) const {
    assert(m_type != WGMatrix4x4Type::Undefine);
    if (m_type == WGMatrix4x4Type::Identity) {
        inverse_matrix = *this;
        return true;
    }
    if (m_type != WGMatrix4x4Type::Unknown) {
        double a00 = m_elements[1][1] * m_elements[2][2] - m_elements[1][2] * m_elements[2][1];
        double a01 = m_elements[1][0] * m_elements[2][2] - m_elements[1][2] * m_elements[2][0];
        double a02 = m_elements[1][0] * m_elements[2][1] - m_elements[1][1] * m_elements[2][0];
        double a03 = 0;
        double a10 = m_elements[0][1] * m_elements[2][2] - m_elements[0][2] * m_elements[2][1];
        double a11 = m_elements[0][0] * m_elements[2][2] - m_elements[0][2] * m_elements[2][0];
        double a12 = m_elements[0][0] * m_elements[2][1] - m_elements[0][1] * m_elements[2][0];
        double a13 = 0;
        double a20 = m_elements[0][1] * m_elements[1][2] - m_elements[0][2] * m_elements[1][1];
        double a21 = m_elements[0][0] * m_elements[1][2] - m_elements[0][2] * m_elements[1][0];
        double a22 = m_elements[0][0] * m_elements[1][1] - m_elements[0][1] * m_elements[1][0];
        double a23 = 0;
        double a30 = m_elements[0][1] * (m_elements[1][2] * m_elements[2][3] - m_elements[1][3] * m_elements[2][2]) -
            m_elements[0][2] * (m_elements[1][1] * m_elements[2][3] - m_elements[1][3] * m_elements[2][1]) +
            m_elements[0][3] * (m_elements[1][1] * m_elements[2][2] - m_elements[1][2] * m_elements[2][1]);
        double a31 = m_elements[0][0] * (m_elements[1][2] * m_elements[2][3] - m_elements[1][3] * m_elements[2][2]) -
            m_elements[0][2] * (m_elements[1][0] * m_elements[2][3] - m_elements[1][3] * m_elements[2][0]) +
            m_elements[0][3] * (m_elements[1][0] * m_elements[2][2] - m_elements[1][2] * m_elements[2][0]);
        double a32 = m_elements[0][0] * (m_elements[1][1] * m_elements[2][3] - m_elements[1][3] * m_elements[2][1]) -
            m_elements[0][1] * (m_elements[1][0] * m_elements[2][3] - m_elements[1][3] * m_elements[2][0]) +
            m_elements[0][3] * (m_elements[1][0] * m_elements[2][1] - m_elements[1][1] * m_elements[2][0]);
        double a33 = m_elements[0][0] * (m_elements[1][1] * m_elements[2][2] - m_elements[1][2] * m_elements[2][1]) -
            m_elements[0][1] * (m_elements[1][0] * m_elements[2][2] - m_elements[1][2] * m_elements[2][0]) +
            m_elements[0][2] * (m_elements[1][0] * m_elements[2][1] - m_elements[1][1] * m_elements[2][0]);
        double det = m_elements[0][0] * a00 - m_elements[0][1] * a01 + m_elements[0][2] * a02 - m_elements[0][3] * a03;
        if (det == 0) {
            return false;
        }
        double d = 1 / det;
        inverse_matrix.m_elements[0][0] = a00 * d;
        inverse_matrix.m_elements[0][1] = -a10 * d;
        inverse_matrix.m_elements[0][2] = a20 * d;
        inverse_matrix.m_elements[0][3] = -a30 * d;
        inverse_matrix.m_elements[1][0] = -a01 * d;
        inverse_matrix.m_elements[1][1] = a11 * d;
        inverse_matrix.m_elements[1][2] = -a21 * d;
        inverse_matrix.m_elements[1][3] = a31 * d;
        inverse_matrix.m_elements[2][0] = a02 * d;
        inverse_matrix.m_elements[2][1] = -a12 * d;
        inverse_matrix.m_elements[2][2] = a22 * d;
        inverse_matrix.m_elements[2][3] = -a32 * d;
        inverse_matrix.m_elements[3][0] = -a03 * d;
        inverse_matrix.m_elements[3][1] = a13 * d;
        inverse_matrix.m_elements[3][2] = -a23 * d;
        inverse_matrix.m_elements[3][3] = a33 * d;
        inverse_matrix.m_type = WGMatrix4x4Type::Unknown;
        return true;
    }
    double a00 = m_elements[1][1] * (m_elements[2][2] * m_elements[3][3] - m_elements[2][3] * m_elements[3][2]) -
        m_elements[1][2] * (m_elements[2][1] * m_elements[3][3] - m_elements[2][3] * m_elements[3][1]) +
        m_elements[1][3] * (m_elements[2][1] * m_elements[3][2] - m_elements[2][2] * m_elements[3][1]);
    double a01 = m_elements[1][0] * (m_elements[2][2] * m_elements[3][3] - m_elements[2][3] * m_elements[3][2]) -
        m_elements[1][2] * (m_elements[2][0] * m_elements[3][3] - m_elements[2][3] * m_elements[3][0]) +
        m_elements[1][3] * (m_elements[2][0] * m_elements[3][2] - m_elements[2][2] * m_elements[3][0]);
    double a02 = m_elements[1][0] * (m_elements[2][1] * m_elements[3][3] - m_elements[2][3] * m_elements[3][1]) -
        m_elements[1][1] * (m_elements[2][0] * m_elements[3][3] - m_elements[2][3] * m_elements[3][0]) +
        m_elements[1][3] * (m_elements[2][0] * m_elements[3][1] - m_elements[2][1] * m_elements[3][0]);
    double a03 = m_elements[1][0] * (m_elements[2][1] * m_elements[3][2] - m_elements[2][2] * m_elements[3][1]) -
        m_elements[1][1] * (m_elements[2][0] * m_elements[3][2] - m_elements[2][2] * m_elements[3][0]) +
        m_elements[1][2] * (m_elements[2][0] * m_elements[3][1] - m_elements[2][1] * m_elements[3][0]);
    double a10 = m_elements[0][1] * (m_elements[2][2] * m_elements[3][3] - m_elements[2][3] * m_elements[3][2]) -
        m_elements[0][2] * (m_elements[2][1] * m_elements[3][3] - m_elements[2][3] * m_elements[3][1]) +
        m_elements[0][3] * (m_elements[2][1] * m_elements[3][2] - m_elements[2][2] * m_elements[3][1]);
    double a11 = m_elements[0][0] * (m_elements[2][2] * m_elements[3][3] - m_elements[2][3] * m_elements[3][2]) -
        m_elements[0][2] * (m_elements[2][0] * m_elements[3][3] - m_elements[2][3] * m_elements[3][0]) +
        m_elements[0][3] * (m_elements[2][0] * m_elements[3][2] - m_elements[2][2] * m_elements[3][0]);
    double a12 = m_elements[0][0] * (m_elements[2][1] * m_elements[3][3] - m_elements[2][3] * m_elements[3][1]) -
        m_elements[0][1] * (m_elements[2][0] * m_elements[3][3] - m_elements[2][3] * m_elements[3][0]) +
        m_elements[0][3] * (m_elements[2][0] * m_elements[3][1] - m_elements[2][1] * m_elements[3][0]);
    double a13 = m_elements[0][0] * (m_elements[2][1] * m_elements[3][2] - m_elements[2][2] * m_elements[3][1]) -
        m_elements[0][1] * (m_elements[2][0] * m_elements[3][2] - m_elements[2][2] * m_elements[3][0]) +
        m_elements[0][2] * (m_elements[2][0] * m_elements[3][1] - m_elements[2][1] * m_elements[3][0]);
    double a20 = m_elements[0][1] * (m_elements[1][2] * m_elements[3][3] - m_elements[1][3] * m_elements[3][2]) -
        m_elements[0][2] * (m_elements[1][1] * m_elements[3][3] - m_elements[1][3] * m_elements[3][1]) +
        m_elements[0][3] * (m_elements[1][1] * m_elements[3][2] - m_elements[1][2] * m_elements[3][1]);
    double a21 = m_elements[0][0] * (m_elements[1][2] * m_elements[3][3] - m_elements[1][3] * m_elements[3][2]) -
        m_elements[0][2] * (m_elements[1][0] * m_elements[3][3] - m_elements[1][3] * m_elements[3][0]) +
        m_elements[0][3] * (m_elements[1][0] * m_elements[3][2] - m_elements[1][2] * m_elements[3][0]);
    double a22 = m_elements[0][0] * (m_elements[1][1] * m_elements[3][3] - m_elements[1][3] * m_elements[3][1]) -
        m_elements[0][1] * (m_elements[1][0] * m_elements[3][3] - m_elements[1][3] * m_elements[3][0]) +
        m_elements[0][3] * (m_elements[1][0] * m_elements[3][1] - m_elements[1][1] * m_elements[3][0]);
    double a23 = m_elements[0][0] * (m_elements[1][1] * m_elements[3][2] - m_elements[1][2] * m_elements[3][1]) -
        m_elements[0][1] * (m_elements[1][0] * m_elements[3][2] - m_elements[1][2] * m_elements[3][0]) +
        m_elements[0][2] * (m_elements[1][0] * m_elements[3][1] - m_elements[1][1] * m_elements[3][0]);
    double a30 = m_elements[0][1] * (m_elements[1][2] * m_elements[2][3] - m_elements[1][3] * m_elements[2][2]) -
        m_elements[0][2] * (m_elements[1][1] * m_elements[2][3] - m_elements[1][3] * m_elements[2][1]) +
        m_elements[0][3] * (m_elements[1][1] * m_elements[2][2] - m_elements[1][2] * m_elements[2][1]);
    double a31 = m_elements[0][0] * (m_elements[1][2] * m_elements[2][3] - m_elements[1][3] * m_elements[2][2]) -
        m_elements[0][2] * (m_elements[1][0] * m_elements[2][3] - m_elements[1][3] * m_elements[2][0]) +
        m_elements[0][3] * (m_elements[1][0] * m_elements[2][2] - m_elements[1][2] * m_elements[2][0]);
    double a32 = m_elements[0][0] * (m_elements[1][1] * m_elements[2][3] - m_elements[1][3] * m_elements[2][1]) -
        m_elements[0][1] * (m_elements[1][0] * m_elements[2][3] - m_elements[1][3] * m_elements[2][0]) +
        m_elements[0][3] * (m_elements[1][0] * m_elements[2][1] - m_elements[1][1] * m_elements[2][0]);
    double a33 = m_elements[0][0] * (m_elements[1][1] * m_elements[2][2] - m_elements[1][2] * m_elements[2][1]) -
        m_elements[0][1] * (m_elements[1][0] * m_elements[2][2] - m_elements[1][2] * m_elements[2][0]) +
        m_elements[0][2] * (m_elements[1][0] * m_elements[2][1] - m_elements[1][1] * m_elements[2][0]);
    double det = m_elements[0][0] * a00 - m_elements[0][1] * a01 + m_elements[0][2] * a02 - m_elements[0][3] * a03;
    if (det == 0) {
        return false;
    }
    double d = 1 / det;
    inverse_matrix.m_elements[0][0] = a00 * d;
    inverse_matrix.m_elements[0][1] = -a10 * d;
    inverse_matrix.m_elements[0][2] = a20 * d;
    inverse_matrix.m_elements[0][3] = -a30 * d;
    inverse_matrix.m_elements[1][0] = -a01 * d;
    inverse_matrix.m_elements[1][1] = a11 * d;
    inverse_matrix.m_elements[1][2] = -a21 * d;
    inverse_matrix.m_elements[1][3] = a31 * d;
    inverse_matrix.m_elements[2][0] = a02 * d;
    inverse_matrix.m_elements[2][1] = -a12 * d;
    inverse_matrix.m_elements[2][2] = a22 * d;
    inverse_matrix.m_elements[2][3] = -a32 * d;
    inverse_matrix.m_elements[3][0] = -a03 * d;
    inverse_matrix.m_elements[3][1] = a13 * d;
    inverse_matrix.m_elements[3][2] = -a23 * d;
    inverse_matrix.m_elements[3][3] = a33 * d;
    inverse_matrix.m_type = WGMatrix4x4Type::Unknown;
    return true;
}

WGMatrix4x4 WGMatrix4x4::MulTranslateMatrix(const WGVector3d& translate) const {
    assert(m_type != WGMatrix4x4Type::Undefine);
    if (translate.X == 0 && translate.Y == 0 && translate.Z == 0) {
        return *this;
    }
    if (m_type == WGMatrix4x4Type::Identity) {
        return BuildTranslation(translate);
    }
    WGMatrix4x4 matrix = *this;
    matrix.m_elements[0][3] += m_elements[0][0] * translate.X + m_elements[0][1] * translate.Y + m_elements[0][2] * translate.Z;
    matrix.m_elements[1][3] += m_elements[1][0] * translate.X + m_elements[1][1] * translate.Y + m_elements[1][2] * translate.Z;
    matrix.m_elements[2][3] += m_elements[2][0] * translate.X + m_elements[2][1] * translate.Y + m_elements[2][2] * translate.Z;
    if (m_type == WGMatrix4x4Type::Unknown) {
        matrix.m_elements[3][3] += m_elements[3][0] * translate.X + m_elements[3][1] * translate.Y + m_elements[3][2] * translate.Z;
    }
    return matrix;
}

WGMatrix4x4 WGMatrix4x4::MulScaleMatrix(const WGVector3d& scale) const {
    assert(m_type != WGMatrix4x4Type::Undefine);
    if (scale.X == 0 && scale.Y == 0 && scale.Z == 0) {
        return *this;
    }
    if (m_type == WGMatrix4x4Type::Identity) {
        return BuildScale(scale);
    }
    WGMatrix4x4 matrix = *this;
    matrix.m_elements[0][0] *= scale.X;
    matrix.m_elements[0][1] *= scale.Y;
    matrix.m_elements[0][2] *= scale.Z;
    matrix.m_elements[1][0] *= scale.X;
    matrix.m_elements[1][1] *= scale.Y;
    matrix.m_elements[1][2] *= scale.Z;
    matrix.m_elements[2][0] *= scale.X;
    matrix.m_elements[2][1] *= scale.Y;
    matrix.m_elements[2][2] *= scale.Z;
    if (m_type == WGMatrix4x4Type::Unknown) {
        matrix.m_elements[3][0] *= scale.X;
        matrix.m_elements[3][1] *= scale.Y;
        matrix.m_elements[3][2] *= scale.Z;
    }
    return matrix;
}

WGMatrix4x4 WGMatrix4x4::BuildIdentity() {
    WGMatrix4x4 matrix;
    matrix.m_elements[0][0] = 1;
    matrix.m_elements[0][1] = 0;
    matrix.m_elements[0][2] = 0;
    matrix.m_elements[0][3] = 0;
    matrix.m_elements[1][0] = 0;
    matrix.m_elements[1][1] = 1;
    matrix.m_elements[1][2] = 0;
    matrix.m_elements[1][3] = 0;
    matrix.m_elements[2][0] = 0;
    matrix.m_elements[2][1] = 0;
    matrix.m_elements[2][2] = 1;
    matrix.m_elements[2][3] = 0;
    matrix.m_elements[3][0] = 0;
    matrix.m_elements[3][1] = 0;
    matrix.m_elements[3][2] = 0;
    matrix.m_elements[3][3] = 1;
    matrix.m_type = WGMatrix4x4Type::Identity;
    return matrix;
}

WGMatrix4x4 WGMatrix4x4::BuildTranslation(const WGVector3d& t) {
    WGMatrix4x4 matrix;
    matrix.m_elements[0][0] = 1;
    matrix.m_elements[0][1] = 0;
    matrix.m_elements[0][2] = 0;
    matrix.m_elements[0][3] = t.X;
    matrix.m_elements[1][0] = 0;
    matrix.m_elements[1][1] = 1;
    matrix.m_elements[1][2] = 0;
    matrix.m_elements[1][3] = t.Y;
    matrix.m_elements[2][0] = 0;
    matrix.m_elements[2][1] = 0;
    matrix.m_elements[2][2] = 1;
    matrix.m_elements[2][3] = t.Z;
    matrix.m_elements[3][0] = 0;
    matrix.m_elements[3][1] = 0;
    matrix.m_elements[3][2] = 0;
    matrix.m_elements[3][3] = 1;
    matrix.m_type = WGMatrix4x4Type::TRS;
    return matrix;
}

WGMatrix4x4 WGMatrix4x4::BuildRotation(const WGQuaternion& r) {
    WGMatrix4x4 matrix;
    matrix.m_elements[0][0] = 1 - 2 * (r.Y * r.Y + r.Z * r.Z);
    matrix.m_elements[0][1] = 2 * (r.X * r.Y - r.W * r.Z);
    matrix.m_elements[0][2] = 2 * (r.X * r.Z + r.W * r.Y);
    matrix.m_elements[0][3] = 0;
    matrix.m_elements[1][0] = 2 * (r.X * r.Y + r.W * r.Z);
    matrix.m_elements[1][1] = 1 - 2 * (r.X * r.X + r.Z * r.Z);
    matrix.m_elements[1][2] = 2 * (r.Y * r.Z - r.W * r.X);
    matrix.m_elements[1][3] = 0;
    matrix.m_elements[2][0] = 2 * (r.X * r.Z - r.W * r.Y);
    matrix.m_elements[2][1] = 2 * (r.Y * r.Z + r.W * r.X);
    matrix.m_elements[2][2] = 1 - 2 * (r.X * r.X + r.Y * r.Y);
    matrix.m_elements[2][3] = 0;
    matrix.m_elements[3][0] = 0;
    matrix.m_elements[3][1] = 0;
    matrix.m_elements[3][2] = 0;
    matrix.m_elements[3][3] = 1;
    matrix.m_type = WGMatrix4x4Type::TRS;
    return matrix;
}

WGMatrix4x4 WGMatrix4x4::BuildScale(const WGVector3d& s) {
    WGMatrix4x4 matrix;
    matrix.m_elements[0][0] = s.X;
    matrix.m_elements[0][1] = 0;
    matrix.m_elements[0][2] = 0;
    matrix.m_elements[0][3] = 0;
    matrix.m_elements[1][0] = 0;
    matrix.m_elements[1][1] = s.Y;
    matrix.m_elements[1][2] = 0;
    matrix.m_elements[1][3] = 0;
    matrix.m_elements[2][0] = 0;
    matrix.m_elements[2][1] = 0;
    matrix.m_elements[2][2] = s.Z;
    matrix.m_elements[2][3] = 0;
    matrix.m_elements[3][0] = 0;
    matrix.m_elements[3][1] = 0;
    matrix.m_elements[3][2] = 0;
    matrix.m_elements[3][3] = 1;
    matrix.m_type = WGMatrix4x4Type::TRS;
    return matrix;
}

WGMatrix4x4 WGMatrix4x4::BuildTRS(const WGVector3d& t, const WGQuaternion& r, const WGVector3d& s) {
    WGMatrix4x4 matrix;
    matrix.m_elements[0][0] = (1 - 2 * (r.Y * r.Y + r.Z * r.Z)) * s.X;
    matrix.m_elements[0][1] = (2 * (r.X * r.Y - r.W * r.Z)) * s.Y;
    matrix.m_elements[0][2] = (2 * (r.X * r.Z + r.W * r.Y)) * s.Z;
    matrix.m_elements[0][3] = t.X;
    matrix.m_elements[1][0] = (2 * (r.X * r.Y + r.W * r.Z)) * s.X;
    matrix.m_elements[1][1] = (1 - 2 * (r.X * r.X + r.Z * r.Z)) * s.Y;
    matrix.m_elements[1][2] = (2 * (r.Y * r.Z - r.W * r.X)) * s.Z;
    matrix.m_elements[1][3] = t.Y;
    matrix.m_elements[2][0] = (2 * (r.X * r.Z - r.W * r.Y)) * s.X;
    matrix.m_elements[2][1] = (2 * (r.Y * r.Z + r.W * r.X)) * s.Y;
    matrix.m_elements[2][2] = (1 - 2 * (r.X * r.X + r.Y * r.Y)) * s.Z;
    matrix.m_elements[2][3] = t.Z;
    matrix.m_elements[3][0] = 0;
    matrix.m_elements[3][1] = 0;
    matrix.m_elements[3][2] = 0;
    matrix.m_elements[3][3] = 1;
    matrix.m_type = WGMatrix4x4Type::TRS;
    return matrix;
}

WGMatrix4x4 WGMatrix4x4::BuildInverseTRS(const WGVector3d& t, const WGQuaternion& r, const WGVector3d& s) {
    double r00 = 1 - 2 * (r.Y * r.Y + r.Z * r.Z);
    double r01 = 2 * (r.X * r.Y - r.W * r.Z);
    double r02 = 2 * (r.X * r.Z + r.W * r.Y);
    double r10 = 2 * (r.X * r.Y + r.W * r.Z);
    double r11 = 1 - 2 * (r.X * r.X + r.Z * r.Z);
    double r12 = 2 * (r.Y * r.Z - r.W * r.X);
    double r20 = 2 * (r.X * r.Z - r.W * r.Y);
    double r21 = 2 * (r.Y * r.Z + r.W * r.X);
    double r22 = 1 - 2 * (r.X * r.X + r.Y * r.Y);
    WGMatrix4x4 matrix;
    matrix.m_elements[0][0] = r00 / s.X;
    matrix.m_elements[0][1] = r10 / s.X;
    matrix.m_elements[0][2] = r20 / s.X;
    matrix.m_elements[0][3] = -(r00 * t.X + r10 * t.Y + r20 * t.Z) / s.X;
    matrix.m_elements[1][0] = r01 / s.Y;
    matrix.m_elements[1][1] = r11 / s.Y;
    matrix.m_elements[1][2] = r21 / s.Y;
    matrix.m_elements[1][3] = -(r01 * t.X + r11 * t.Y + r21 * t.Z) / s.Y;
    matrix.m_elements[2][0] = r02 / s.Z;
    matrix.m_elements[2][1] = r12 / s.Z;
    matrix.m_elements[2][2] = r22 / s.Z;
    matrix.m_elements[2][3] = -(r02 * t.X + r12 * t.Y + r22 * t.Z) / s.Z;
    matrix.m_elements[3][0] = 0;
    matrix.m_elements[3][1] = 0;
    matrix.m_elements[3][2] = 0;
    matrix.m_elements[3][3] = 1;
    matrix.m_type = WGMatrix4x4Type::MultTRS;
    return matrix;
}

WGMatrix4x4 WGMatrix4x4::BuildOrtho(double l, double r, double b, double t, double n, double f) {
    WGMatrix4x4 matrix;
    double w = r - l;
    double h = t - b;
    double d = f - n;
    matrix.m_elements[0][0] = 2 / w;
    matrix.m_elements[0][1] = 0;
    matrix.m_elements[0][2] = 0;
    matrix.m_elements[0][3] = -(r + l) / w;
    matrix.m_elements[1][0] = 0;
    matrix.m_elements[1][1] = 2 / h;
    matrix.m_elements[1][2] = 0;
    matrix.m_elements[1][3] = -(t + b) / h;
    matrix.m_elements[2][0] = 0;
    matrix.m_elements[2][1] = 0;
    matrix.m_elements[2][2] = -2 / d;
    matrix.m_elements[2][3] = -(f + n) / d;
    matrix.m_elements[3][0] = 0;
    matrix.m_elements[3][1] = 0;
    matrix.m_elements[3][2] = 0;
    matrix.m_elements[3][3] = 1;
    matrix.m_type = WGMatrix4x4Type::TRS;
    return matrix;
}

WGMatrix4x4 WGMatrix4x4::BuildFrustum(double l, double r, double b, double t, double n, double f) {
    WGMatrix4x4 matrix;
    double w = r - l;
    double h = t - b;
    double d = f - n;
    matrix.m_elements[0][0] = 2 * n / w;
    matrix.m_elements[0][1] = 0;
    matrix.m_elements[0][2] = (r + l) / w;
    matrix.m_elements[0][3] = 0;
    matrix.m_elements[1][0] = 0;
    matrix.m_elements[1][1] = 2 * n / h;
    matrix.m_elements[1][2] = (t + b) / h;
    matrix.m_elements[1][3] = 0;
    matrix.m_elements[2][0] = 0;
    matrix.m_elements[2][1] = 0;
    matrix.m_elements[2][2] = -(f + n) / d;
    matrix.m_elements[2][3] = -2 * f * n / d;
    matrix.m_elements[3][0] = 0;
    matrix.m_elements[3][1] = 0;
    matrix.m_elements[3][2] = -1;
    matrix.m_elements[3][3] = 0;
    matrix.m_type = WGMatrix4x4Type::Unknown;
    return matrix;
}

WGMatrix4x4 WGMatrix4x4::Identity = BuildIdentity();