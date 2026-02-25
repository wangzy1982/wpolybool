#include "WGNurbsCurve2d.h"
#include <assert.h>
#include "wsbasis.h"
#include "WGNurbs.h"
#include "WGIntegration.h"

WG_OBJECT_TYPE_IMP(WGNurbsCurve2d, WGCurve2d);

WGNurbsCurve2d::WGNurbsCurve2d() : 
    m_degree(3),
    m_control_point_count(0),
    m_knots(nullptr),
    m_control_points(nullptr),
    m_weights(nullptr) {
}

WGNurbsCurve2d::WGNurbsCurve2d(int degree, int control_point_count, bool weights) :
    m_degree(degree),
    m_control_point_count(control_point_count) {
    m_knots = new double[degree + control_point_count + 1];
    m_control_points = new WGVector2d[control_point_count];
    if (weights) {
        m_weights = new double[control_point_count];
    }
    else {
        m_weights = nullptr;
    }
}

WGNurbsCurve2d::WGNurbsCurve2d(int degree, int control_point_count, const double* knots, const WGVector2d* control_points, const double* weights) :
    m_degree(degree),
    m_control_point_count(control_point_count) {
    m_knots = new double[degree + control_point_count + 1];
    memcpy(m_knots, knots, (degree + control_point_count + 1) * sizeof(double));
    m_control_points = new WGVector2d[control_point_count];
    memcpy(m_control_points, control_points, control_point_count * sizeof(WGVector2d));
    if (weights) {
        m_weights = new double[control_point_count];
        memcpy(m_weights, weights, control_point_count * sizeof(double));
    }
    else {
        m_weights = nullptr;
    }
}

WGNurbsCurve2d::~WGNurbsCurve2d() {
    delete[] m_knots;
    delete[] m_control_points;
    delete[] m_weights;
}

WGVector2d WGNurbsCurve2d::GetStartPoint() const {
    assert(m_degree > 0);
    return CalculateD0(0, m_knots[m_degree]);
}

WGVector2d WGNurbsCurve2d::GetEndPoint() const {
    assert(m_degree > 0);
    return CalculateD0(m_control_point_count - m_degree - 1, m_knots[m_control_point_count]);
}

WGNurbsCurve2d* WGNurbsCurve2d::InsertKnot(int piece_index, double t, int insert_count) const {
    WGNurbsCurve2d* curve = new WGNurbsCurve2d(m_degree, m_control_point_count + insert_count, m_weights);
    insert_knots(curve->m_knots, curve->m_control_points, curve->m_weights, m_knots, m_control_points, m_weights, m_degree,
        m_control_point_count, piece_index, t, insert_count);
    return curve;
}

int WGNurbsCurve2d::GetDegree() const {
    return m_degree;
}

int WGNurbsCurve2d::GetControlPointCount() const {
    return m_control_point_count;
}

int WGNurbsCurve2d::GetKnotCount() const {
    return m_control_point_count + m_degree + 1;
}

const double* WGNurbsCurve2d::GetKnots() const {
    return m_knots;
}

const WGVector2d* WGNurbsCurve2d::GetControlPoints() const {
    return m_control_points;
}

const double* WGNurbsCurve2d::GetWeights() const {
    return m_weights;
}

int WGNurbsCurve2d::GetPieceCount() const {
    return m_control_point_count - m_degree;
}

WSInterval WGNurbsCurve2d::GetPieceDomain(int index) const {
    assert(m_degree > 0);
    int i = m_degree + index;
    return WSInterval(m_knots[i], m_knots[i + 1]);
}

WGVector2d WGNurbsCurve2d::CalculateD0(int piece_index, double t) const {
    assert(m_degree > 0);
    if (m_weights) {
        double basis_buffer[g_max_nurbs_basis_buffer_size];
        calculate_nurbs_basis(basis_buffer, m_knots, m_degree, m_degree, piece_index, t);
        double* basis = basis_buffer + get_nurbs_basis_buffer_size(m_degree - 1);
        WGVector2d d0(0, 0);
        double w = 0;
        for (int i = 0; i <= m_degree; ++i) {
            int j = piece_index + i;
            double d = m_weights[j] * basis[i];
            d0 = d0 + m_control_points[j] * d;
            w = w + d;
        }
        return d0 / w;
    }
    else {
        double basis_buffer[g_max_nurbs_basis_buffer_size];
        calculate_nurbs_basis(basis_buffer, m_knots, m_degree, m_degree, piece_index, t);
        double* basis = basis_buffer + get_nurbs_basis_buffer_size(m_degree - 1);
        WGVector2d d0(0, 0);
        for (int i = 0; i <= m_degree; ++i) {
            int j = piece_index + i;
            d0 = d0 + m_control_points[j] * basis[i];
        }
        return d0;
    }
}

WGVector2d WGNurbsCurve2d::CalculateD1(int piece_index, double t) const {
    assert(m_degree > 0);
    if (m_weights) {
        double basis_buffer[g_max_nurbs_basis_buffer_size];
        double basis_d1[g_max_nurbs_degree + 1];
        calculate_nurbs_basis(basis_buffer, m_knots, m_degree, m_degree, piece_index, t);
        calculate_nurbs_basis_d1(basis_d1, basis_buffer, m_knots, m_degree, piece_index);
        double* basis = basis_buffer + get_nurbs_basis_buffer_size(m_degree - 1);
        WGVector2d d0(0, 0);
        double w0 = 0;
        WGVector2d d1(0, 0);
        double w1 = 0;
        for (int i = 0; i <= m_degree; ++i) {
            int j = piece_index + i;
            double d = m_weights[j] * basis[i];
            d0 = d0 + m_control_points[j] * d;
            w0 = w0 + d;
            d = m_weights[j] * basis_d1[i];
            d1 = d1 + m_control_points[j] * d;
            w1 = w1 + d;
        }
        return (d1 * w0 - d0 * w1) / (w0 * w0);
    }
    else {
        double basis_buffer[g_max_nurbs_basis_buffer_size];
        double basis_d1[g_max_nurbs_degree + 1];
        calculate_nurbs_basis(basis_buffer, m_knots, m_degree, m_degree - 1, piece_index, t);
        calculate_nurbs_basis_d1(basis_d1, basis_buffer, m_knots, m_degree, piece_index);
        WGVector2d d1(0, 0);
        for (int i = 0; i <= m_degree; ++i) {
            int j = piece_index + i;
            d1 = d1 + m_control_points[j] * basis_d1[i];
        }
        return d1;
    }
}

double WGNurbsCurve2d::CalculateVariableEpsilon(int piece_index, double distance_epsilon) const {
    assert(m_degree > 0);
    double d = 0;
    for (int i = 0; i < m_degree; ++i) {
        int j = piece_index + i;
        d += (m_control_points[j + 1] - m_control_points[j]).Length();
    }
    int j = m_degree + piece_index;
    double dt = m_knots[j + 1] - m_knots[j];
    if (d == 0) {
        return dt;
    }
    d = distance_epsilon / d;
    if (d >= 1) {
        return dt;
    }
    return d * dt;
}

WGBox2d WGNurbsCurve2d::CalculateBox() const {
    WGBox2d box;
    for (int i = 0; i < m_control_point_count; ++i) {
        box.Merge(m_control_points[i]);
    }
    return box;
}

double WGNurbsCurve2d::CalculateAreaBelow(int piece_index, double start_t, double end_t) const {
    assert(m_degree > 0);
    class Integrand : public WGIntegrand1V {
    public:
        Integrand(const WGNurbsCurve2d* curve, int piece_index) :
            m_curve(curve),
            m_piece_index(piece_index) {
        }

        virtual double Calculate(double t) const {
            if (m_curve->m_weights) {
                double basis_buffer[g_max_nurbs_basis_buffer_size];
                double basis_d1[g_max_nurbs_degree + 1];
                calculate_nurbs_basis(basis_buffer, m_curve->m_knots, m_curve->m_degree, m_curve->m_degree, m_piece_index, t);
                calculate_nurbs_basis_d1(basis_d1, basis_buffer, m_curve->m_knots, m_curve->m_degree, m_piece_index);
                double* basis = basis_buffer + get_nurbs_basis_buffer_size(m_curve->m_degree - 1);
                double y0 = 0;
                double x0 = 0;
                double x1 = 0;
                double w0 = 0;
                double w1 = 0;
                for (int i = 0; i <= m_curve->m_degree; ++i) {
                    int j = m_piece_index + i;
                    double d = m_curve->m_weights[j] * basis[i];
                    x0 = x0 + m_curve->m_control_points[j].X * d;
                    y0 = y0 + m_curve->m_control_points[j].Y * d;
                    w0 = w0 + d;
                    d = m_curve->m_weights[j] * basis_d1[i];
                    x1 = x1 + m_curve->m_control_points[j].X * d;
                    w1 = w1 + d;
                }
                return y0 * (x1 * w0 - x0 * w1) / (w0 * w0 * w0);
            }
            else {
                double basis_buffer[g_max_nurbs_basis_buffer_size];
                double basis_d1[g_max_nurbs_degree + 1];
                calculate_nurbs_basis(basis_buffer, m_curve->m_knots, m_curve->m_degree, m_curve->m_degree, m_piece_index, t);
                calculate_nurbs_basis_d1(basis_d1, basis_buffer, m_curve->m_knots, m_curve->m_degree, m_piece_index);
                double* basis = basis_buffer + get_nurbs_basis_buffer_size(m_curve->m_degree - 1);
                double y0 = 0;
                double x1 = 0;
                for (int i = 0; i <= m_curve->m_degree; ++i) {
                    int j = m_piece_index + i;
                    y0 = y0 + m_curve->m_control_points[j].Y * basis[i];
                    x1 = x1 + m_curve->m_control_points[j].X * basis_d1[i];
                }
                return y0 * x1;
            }
        }
    private:
        const WGNurbsCurve2d* m_curve;
        int m_piece_index;
    };
    Integrand integrand(this, piece_index);
    return integrate(&integrand, start_t, end_t, 1E-6);
}

WGCurve2d* WGNurbsCurve2d::Clone() const {
    return new WGNurbsCurve2d(m_degree, m_control_point_count, m_knots, m_control_points, m_weights);
}

WGCurve2d* WGNurbsCurve2d::CreateSubCurve(int start_piece_index, double start_t, int end_piece_index, double end_t) const {
    assert(m_degree > 0);
    assert(end_piece_index > start_piece_index || (end_piece_index == start_piece_index && end_t > start_t));
    int knot_count = m_control_point_count + m_degree + 1;
    int start_insert_count = m_degree + 1;
    int knot_index = start_piece_index + m_degree;
    int i = knot_index;
    while (i >= 0) {
        if (m_knots[i] != start_t) {
            break;
        }
        --start_insert_count;
        --i;
    }
    i = knot_index + 1;
    while (i < knot_count) {
        if (m_knots[i] != start_t) {
            break;
        }
        --start_insert_count;
        ++i;
    }
    start_piece_index = i - 1 - m_degree;
    int end_insert_count = m_degree + 1;
    knot_index = end_piece_index + m_degree;
    i = knot_index;
    while (i >= 0) {
        if (m_knots[i] != end_t) {
            break;
        }
        --end_insert_count;
        --i;
    }
    i = knot_index + 1;
    while (i < knot_count) {
        if (m_knots[i] != end_t) {
            break;
        }
        --end_insert_count;
        ++i;
    }
    end_piece_index = i - 1 - m_degree;
    if (start_insert_count <= 0 && end_insert_count <= 0) {
        int control_point_count = end_piece_index + end_insert_count - start_piece_index;
        return new WGNurbsCurve2d(m_degree, control_point_count, m_knots + start_piece_index, 
            m_control_points + start_piece_index, m_weights + start_piece_index);
    }
    int n = m_control_point_count;
    if (start_insert_count > 0) {
        n += start_insert_count;
    }
    if (end_insert_count > 0) {
        n += end_insert_count;
    }
    double* knots = new double[n + m_degree + 1];
    memcpy(knots, m_knots, (m_control_point_count + m_degree + 1) * sizeof(double));
    WGVector2d* control_points = new WGVector2d[n];
    memcpy(control_points, m_control_points, m_control_point_count * sizeof(WGVector2d));
    double* weights;
    if (m_weights) {
        weights = new double[n];
        memcpy(weights, m_weights, m_control_point_count * sizeof(double));
    }
    else {
        weights = nullptr;
    }
    if (start_insert_count > 0) {
        insert_knots(knots, control_points, weights, knots, control_points, weights, m_degree, 
            m_control_point_count, start_piece_index, start_t, start_insert_count);
        end_piece_index += start_insert_count;
        start_piece_index += start_insert_count;
    }
    else {
        start_insert_count = 0;
    }
    if (end_insert_count > 0) {
        insert_knots(knots, control_points, weights, knots, control_points, weights, m_degree,
            m_control_point_count + start_insert_count, end_piece_index, end_t, end_insert_count);
        end_piece_index += end_insert_count;
    }
    int control_point_count = end_piece_index - start_piece_index;
    WGCurve2d* curve;
    if (weights) {
        curve = new WGNurbsCurve2d(m_degree, control_point_count, knots + start_piece_index,
            control_points + start_piece_index, weights + start_piece_index);
    }
    else {
        curve = new WGNurbsCurve2d(m_degree, control_point_count, knots + start_piece_index,
            control_points + start_piece_index, nullptr);
    }
    delete[] knots;
    delete[] control_points;
    delete[] weights;
    return curve;
}

void WGNurbsCurve2d::Reverse() {
    int m = m_control_point_count / 2;
    for (int i = 0; i < m; ++i) {
        int j = m_control_point_count - 1 - i;
        WGVector2d t = m_control_points[i];
        m_control_points[i] = m_control_points[j];
        m_control_points[j] = t;
    }
    if (m_weights) {
        for (int i = 0; i < m; ++i) {
            int j = m_control_point_count - 1 - i;
            double t = m_weights[i];
            m_weights[i] = m_weights[j];
            m_weights[j] = t;
        }
    }
    int n = m_control_point_count + m_degree;
    double d = m_knots[n];
    m = (n + 1) / 2;
    for (int i = 0; i < m; ++i) {
        int j = n - i;
        double t = m_knots[i];
        m_knots[i] = m_knots[j];
        m_knots[j] = t;
        m_knots[i] = d - m_knots[i];
        m_knots[j] = d - m_knots[j];
    }
    if (2 * m < n + 1) {
        m_knots[m] = d - m_knots[m];
    }
}

void WGNurbsCurve2d::Linearize(std::vector<WGVector2d>& vertices, bool is_contiguous_with_prev,
    double limit_distance_epsilon, double limit_angle_epsilon, double distance_epsilon, double angle_epsilon) const {
    assert(m_degree > 0);
    class Linearizer {
    public:
        Linearizer(const WGNurbsCurve2d* curve, std::vector<WGVector2d>& vertices,
            double limit_distance_epsilon, double limit_angle_epsilon, double distance_epsilon, double angle_epsilon) :
            m_curve(curve),
            m_vertices(vertices),
            m_limit_distance_epsilon(limit_distance_epsilon),
            m_distance_epsilon(distance_epsilon) {
            if (limit_angle_epsilon >= g_pi * 0.5) {
                m_limit_sin_epsilon = 1;
            }
            else {
                m_limit_sin_epsilon = sin(limit_angle_epsilon);
            }
            if (angle_epsilon >= g_pi * 0.5) {
                m_sin_epsilon = 1;
            }
            else {
                m_sin_epsilon = sin(angle_epsilon);
            }
        }

        void Execute(int piece_index, double start_t, double end_t, const WGVector2d& start_point, const WGVector2d& end_point) {
            WGVector2d start_dir = CalculateDir(piece_index, start_t);
            WGVector2d end_dir = CalculateDir(piece_index, end_t);
            double middle_t = (start_t + end_t) * 0.5;
            WGVector2d middle_point = m_curve->CalculateD0(piece_index, middle_t);
            WGVector2d middle_dir = CalculateDir(piece_index, middle_t);
            if (CheckLinear(start_point, start_dir, end_point, end_dir, middle_point, middle_dir)) {
                Execute(piece_index, start_t, end_t, middle_t, start_point, start_dir, end_point, end_dir, middle_point, middle_dir, 1);
            }
            else {
                Execute(piece_index, start_t, middle_t, start_point, start_dir, middle_point, middle_dir, 1);
                Execute(piece_index, middle_t, end_t, middle_point, middle_dir, end_point, end_dir, 1);
            }
        }

    private:
        WGVector2d CalculateDir(int piece_index, double t) {
            WGVector2d dir = m_curve->CalculateD1(piece_index, t);
            if (dir.Normalize(g_double_epsilon) <= g_double_epsilon) {
                dir.X = 0;
                dir.Y = 0;
            }
            return dir;
        }

        bool CheckLinear(const WGVector2d& start_point, const WGVector2d& start_dir, const WGVector2d& end_point, 
            const WGVector2d& end_dir, const WGVector2d& middle_point, const WGVector2d& middle_dir) {
            if (start_dir.Dot(end_dir) < 0) {
                return false;
            }
            if (start_dir.Dot(middle_dir) < 0) {
                return false;
            }
            if (middle_dir.Dot(end_dir) < 0) {
                return false;
            }
            double s = abs(start_dir.Cross(end_dir));
            double s2 = abs(start_dir.Cross(middle_dir));
            if (s2 > s) {
                s = s2;
            }
            s2 = abs(middle_dir.Cross(end_dir));
            if (s2 > s) {
                s = s2;
            }
            double d;
            WGVector2d vt = end_point - start_point;
            if (vt.Normalize(g_double_epsilon) <= g_double_epsilon) {
                d = (middle_point - start_point).Length();
            }
            else {
                d = abs(vt.Cross(middle_point - start_point));
            }
            if (s > m_limit_sin_epsilon) {
                return false;
            }
            if (d > m_limit_distance_epsilon) {
                return false;
            }
            return s <= m_sin_epsilon || d <= m_distance_epsilon;
        }

        void Execute(int piece_index, double start_t, double end_t, const WGVector2d& start_point, const WGVector2d& start_dir,
            const WGVector2d& end_point, const WGVector2d& end_dir, int depth) {
            double middle_t = (start_t + end_t) * 0.5;
            WGVector2d middle_point = m_curve->CalculateD0(piece_index, middle_t);
            WGVector2d middle_dir = CalculateDir(piece_index, middle_t);
            if (CheckLinear(start_point, start_dir, end_point, end_dir, middle_point, middle_dir)) {
                if (depth >= m_max_depth) {
                    m_vertices.push_back(end_point);
                }
                else {
                    Execute(piece_index, start_t, end_t, middle_t, start_point, start_dir, end_point, end_dir, middle_point, middle_dir, depth + 1);
                }
            }
            else {
                if (depth >= m_max_depth) {
                    m_vertices.push_back(middle_point);
                    m_vertices.push_back(end_point);
                }
                else {
                    Execute(piece_index, start_t, middle_t, start_point, start_dir, middle_point, middle_dir, depth + 1);
                    Execute(piece_index, middle_t, end_t, middle_point, middle_dir, end_point, end_dir, depth + 1);
                }
            }
        }

        void Execute(int piece_index, double start_t, double end_t, double middle_t, const WGVector2d& start_point, const WGVector2d& start_dir,
            const WGVector2d& end_point, const WGVector2d& end_dir, const WGVector2d& middle_point, const WGVector2d& middle_dir, int depth) {
            double middle_t1 = (start_t + middle_t) * 0.5;
            WGVector2d middle_point1 = m_curve->CalculateD0(piece_index, middle_t1);
            WGVector2d middle_dir1 = CalculateDir(piece_index, middle_t1);
            double middle_t2 = (middle_t + end_t) * 0.5;
            WGVector2d middle_point2 = m_curve->CalculateD0(piece_index, middle_t2);
            WGVector2d middle_dir2 = CalculateDir(piece_index, middle_t2);
            bool b1 = CheckLinear(start_point, start_dir, middle_point, middle_dir, middle_point1, middle_dir1);
            bool b2 = CheckLinear(middle_point, middle_dir, end_point, end_dir, middle_point2, middle_dir2);
            if (b1) {
                if (b2) {
                    m_vertices.push_back(end_point);
                }
                else {
                    if (depth >= m_max_depth) {
                        m_vertices.push_back(middle_point);
                        m_vertices.push_back(middle_point2);
                        m_vertices.push_back(end_point);
                    }
                    else {
                        Execute(piece_index, start_t, middle_t, middle_t1, start_point, start_dir, 
                            middle_point, middle_dir, middle_point1, middle_dir1, depth + 1);
                        Execute(piece_index, middle_t, middle_t2, middle_point, middle_dir, middle_point2, middle_dir2, depth + 1);
                        Execute(piece_index, middle_t2, end_t, middle_point2, middle_dir2, end_point, end_dir, depth + 1);
                    }
                }
            }
            else {
                if (b2) {
                    if (depth >= m_max_depth) {
                        m_vertices.push_back(middle_point1);
                        m_vertices.push_back(middle_point);
                        m_vertices.push_back(end_point);
                    }
                    else {
                        Execute(piece_index, start_t, middle_t1, start_point, start_dir, middle_point1, middle_dir1, depth + 1);
                        Execute(piece_index, middle_t1, middle_t, middle_point1, middle_dir1, middle_point, middle_dir, depth + 1);
                        Execute(piece_index, middle_t, middle_t2, end_t, middle_point, middle_dir,
                            middle_point2, middle_dir2, end_point, end_dir, depth + 1);
                    }
                }
                else {
                    if (depth >= m_max_depth) {
                        m_vertices.push_back(middle_point1);
                        m_vertices.push_back(middle_point);
                        m_vertices.push_back(middle_point2);
                        m_vertices.push_back(end_point);
                    }
                    else {
                        Execute(piece_index, start_t, middle_t1, start_point, start_dir, middle_point1, middle_dir1, depth + 1);
                        Execute(piece_index, middle_t1, middle_t, middle_point1, middle_dir1, middle_point, middle_dir, depth + 1);
                        Execute(piece_index, middle_t, middle_t2, middle_point, middle_dir, middle_point2, middle_dir2, depth + 1);
                        Execute(piece_index, middle_t2, end_t, middle_point2, middle_dir2, end_point, end_dir, depth + 1);
                    }
                }
            }
        }
    private:
        const WGNurbsCurve2d* m_curve;
        std::vector<WGVector2d>& m_vertices;
        double m_limit_distance_epsilon;
        double m_limit_sin_epsilon;
        double m_distance_epsilon;
        double m_sin_epsilon;
        const int m_max_depth = 10;
    };
    WGVector2d start_point = GetStartPoint();
    if (vertices.size() == 0 || !is_contiguous_with_prev) {
        vertices.push_back(start_point);
    }
    Linearizer linearizer(this, vertices, limit_distance_epsilon, limit_angle_epsilon, distance_epsilon, angle_epsilon);
    for (int i = 0; i < GetPieceCount(); ++i) {
        WSInterval t_domain = GetPieceDomain(i);
        if (t_domain.Length() <= g_double_epsilon) {
            continue;
        }
        WGVector2d end_point = CalculateD0(i, t_domain.Max);
        linearizer.Execute(i, t_domain.Min, t_domain.Max, start_point, end_point);
        start_point = end_point;
    }
}

void WGNurbsCurve2d::Move(const WGVector2d& vt) {
    for (int i = 0; i < m_control_point_count; ++i) {
        m_control_points[i] = m_control_points[i] + vt;
    }
}

WGNurbsCurve2d* WGNurbsCurve2d::BuildByControlPoints(int control_point_count, const WGVector2d* control_points) {
    if (control_point_count == 0) {
        return nullptr;
    }
    if (control_point_count == 1) {
        WGNurbsCurve2d* curve = new WGNurbsCurve2d(1, 2, false);
        curve->m_control_points[0] = control_points[0];
        curve->m_control_points[1] = control_points[0];
        curve->m_knots[0] = 0;
        curve->m_knots[1] = 0;
        curve->m_knots[2] = 1;
        curve->m_knots[3] = 1;
        return curve;
    }
    if (control_point_count == 2) {
        WGNurbsCurve2d* curve = new WGNurbsCurve2d(1, 2, false);
        curve->m_control_points[0] = control_points[0];
        curve->m_control_points[1] = control_points[1];
        curve->m_knots[0] = 0;
        curve->m_knots[1] = 0;
        curve->m_knots[2] = 1;
        curve->m_knots[3] = 1;
        return curve;
    }
    if (control_point_count == 3) {
        WGNurbsCurve2d* curve = new WGNurbsCurve2d(2, 3, false);
        curve->m_control_points[0] = control_points[0];
        curve->m_control_points[1] = control_points[1];
        curve->m_control_points[2] = control_points[2];
        curve->m_knots[0] = 0;
        curve->m_knots[1] = 0;
        curve->m_knots[2] = 0;
        curve->m_knots[3] = 1;
        curve->m_knots[4] = 1;
        curve->m_knots[5] = 1;
        return curve;
    }
    WGNurbsCurve2d* curve = new WGNurbsCurve2d(3, control_point_count, false);
    memcpy(curve->m_control_points, control_points, control_point_count * sizeof(WGVector2d));
    curve->m_knots[0] = 0;
    curve->m_knots[1] = 0;
    curve->m_knots[2] = 0;
    for (int i = 0; i <= control_point_count - 3; ++i) {
        curve->m_knots[i + 3] = i;
    }
    curve->m_knots[control_point_count + 1] = curve->m_knots[control_point_count];
    curve->m_knots[control_point_count + 2] = curve->m_knots[control_point_count];
    curve->m_knots[control_point_count + 3] = curve->m_knots[control_point_count];
    return curve;
}

WGNurbsCurve2d::BezierIterator::BezierIterator(const WGNurbsCurve2d* curve, int piece_index, int piece_count) :
    m_curve(curve) {
    m_is_bezier = false;
    int total_piece_count = curve->GetPieceCount();
    assert(piece_index < total_piece_count && piece_index + piece_count <= total_piece_count);
    m_piece_index = piece_index;
    m_piece_count = piece_count;
    if (total_piece_count == 1) {
        m_is_bezier = true;
        int degree = curve->GetDegree();
        const double* knots = curve->GetKnots();
        for (int i = 0; i < degree; ++i) {
            if (knots[i] != knots[degree]) {
                m_is_bezier = false;
                break;
            }
        }
        if (m_is_bezier) {
            for (int i = 0; i < degree; ++i) {
                if (knots[degree + 1] != knots[degree + 1 + i]) {
                    m_is_bezier = false;
                    break;
                }
            }
        }
    }
    m_offset = 0;
    int degree = curve->GetDegree();
    if (m_is_bezier) {
        m_knot_buffer = (double*)m_curve->GetKnots();
        m_control_point_buffer = (WGVector2d*)m_curve->GetControlPoints();
        m_weight_buffer = (double*)m_curve->GetWeights();
        m_control_point_count = degree + 1;
    }
    else {
        m_control_point_count = piece_count + degree;
        m_knot_buffer = new double[(piece_count + 3) * (degree + 1)];
        memcpy(m_knot_buffer, curve->GetKnots() + piece_index, (m_control_point_count + degree + 1) * sizeof(double));
        m_control_point_buffer = new WGVector2d[(piece_count + 2) * (degree + 1)];
        memcpy(m_control_point_buffer, curve->GetControlPoints() + piece_index, m_control_point_count * sizeof(WGVector2d));
        if (curve->GetWeights()) {
            m_weight_buffer = new double[(piece_count + 2) * (degree + 1)];
            memcpy(m_weight_buffer, curve->GetWeights() + piece_index, m_control_point_count * sizeof(double));
        }
        else {
            m_weight_buffer = nullptr;
        }        
        int n = 1;
        while (n <= degree) {
            if (m_knot_buffer[degree - n] != m_knot_buffer[degree]) {
                break;
            }
            ++n;
        }
        int m = 0;
        while (m + n <= degree) {
            if (m_knot_buffer[degree + m + 1] != m_knot_buffer[degree]) {
                break;
            }
            ++m;
        }
        n += m;
        if (n <= degree) {
            int insert_count = degree + 1 - n;
            insert_knots(m_knot_buffer, m_control_point_buffer, m_weight_buffer,
                m_knot_buffer, m_control_point_buffer, m_weight_buffer, degree,
                m_control_point_count, 0, m_knot_buffer[degree], insert_count);
            m_control_point_count += insert_count;
            m_offset = insert_count;
        }
    }
    First();
}

WGNurbsCurve2d::BezierIterator::~BezierIterator() {
    if (!m_is_bezier) {
        delete[] m_knot_buffer;
        delete[] m_control_point_buffer;
        delete[] m_weight_buffer;
    }
}

void WGNurbsCurve2d::BezierIterator::First() {
    m_current_knots = m_knot_buffer + m_offset;
    m_current_control_points = m_control_point_buffer + m_offset;
    if (m_weight_buffer) {
        m_current_weights = m_weight_buffer + m_offset;
    }
    else {
        m_current_weights = nullptr;
    }
    SkipZero();
    if (!m_is_bezier) {
        InsertNext();
    }
    m_current_piece_index = m_piece_index;
}

bool WGNurbsCurve2d::BezierIterator::Eof() {
    int degree = m_curve->GetDegree();
    return m_control_point_count - (int)(m_current_control_points - m_control_point_buffer) <= degree;
}

void WGNurbsCurve2d::BezierIterator::Next() {
    int n = m_curve->GetDegree() + 1;
    m_current_knots += n;
    m_current_control_points += n;
    if (m_current_weights) {
        m_current_weights += n;
    }
    SkipZero();
    InsertNext();
}

int WGNurbsCurve2d::BezierIterator::GetCurrentPieceIndex() const {
    WSInterval current_domain = GetCurrentDomain();
    while (m_current_piece_index < m_piece_index + m_piece_count) {
        WSInterval domain = m_curve->GetPieceDomain(m_current_piece_index);
        if (domain.Min == current_domain.Min && domain.Max == current_domain.Max) {
            break;
        }
        ++m_current_piece_index;
    }
    return m_current_piece_index;
}

WSInterval WGNurbsCurve2d::BezierIterator::GetCurrentDomain() const {
    int degree = m_curve->GetDegree();
    return WSInterval(m_current_knots[degree], m_current_knots[degree + 1]);
}

const WGVector2d* WGNurbsCurve2d::BezierIterator::GetCurrentColtrolPoints() const {
    return m_current_control_points;
}

const double* WGNurbsCurve2d::BezierIterator::GetCurrentWeights() const {
    return m_current_weights;
}

void WGNurbsCurve2d::BezierIterator::SkipZero() {
    while (!Eof()) {
        int degree = m_curve->GetDegree();
        if (m_current_knots[degree] != m_current_knots[degree + 1]) {
            break;
        }
        m_current_knots += 1;
        m_current_control_points += 1;
        if (m_current_weights) {
            m_current_weights += 1;
        }
    }
}

void WGNurbsCurve2d::BezierIterator::InsertNext() {
    if (!Eof()) {
        int degree = m_curve->GetDegree();
        int n = 1;
        while (n <= degree) {
            if (m_current_knots[degree + n + 1] != m_current_knots[degree + n]) {
                break;
            }
            ++n;
        }
        if (n <= degree) {
            int insert_count = degree + 1 - n;
            insert_knots(m_current_knots, m_current_control_points, m_current_weights,
                m_current_knots, m_current_control_points, m_current_weights, degree,
                m_control_point_count - (int)(m_current_control_points - m_control_point_buffer), 
                1, m_current_knots[degree + 1], insert_count);
            m_control_point_count += insert_count;
        }
    }
}