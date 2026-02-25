#ifndef _WG_NURBS_
#define _WG_NURBS_

#include <assert.h>

inline int get_nurbs_basis_buffer_size(int degree) {
    return (degree + 2) * (degree + 1) / 2;
}

const int g_max_nurbs_degree = 15;

const int g_max_nurbs_basis_buffer_size = (g_max_nurbs_degree + 2) * (g_max_nurbs_degree + 1) / 2;;

inline void calculate_nurbs_basis(double* nurbs_basis_buffer, double* knots, int degree, int dst_degree, int index, double t) {
    nurbs_basis_buffer[0] = 1;
    double* prev_basis = nurbs_basis_buffer;
    double* curr_basis = prev_basis + 1;
    for (int k = 1; k <= dst_degree; ++k) {
        int i = degree + index - k;
        curr_basis[0] = (knots[i + k + 1] - t) / (knots[i + k + 1] - knots[i + 1]) * prev_basis[0];
        for (int j = 1; j < k; ++j) {
            i = degree + index - k + j;
            curr_basis[j] = (t - knots[i]) / (knots[i + k] - knots[i]) * prev_basis[j - 1] +
                (knots[i + k + 1] - t) / (knots[i + k + 1] - knots[i + 1]) * prev_basis[j];
        }
        i = degree + index;
        curr_basis[k] = (t - knots[i]) / (knots[i + k] - knots[i]) * prev_basis[k - 1];
        prev_basis = curr_basis;
        curr_basis = prev_basis + k + 1;
    }
}

inline void calculate_nurbs_basis_d1(double* nurbs_basis_d1_buffer, double* nurbs_basis_buffer, double* knots, int degree, int index) {
    if (degree == 1) {
        nurbs_basis_d1_buffer[0] = 0;
        nurbs_basis_d1_buffer[1] = 0;
        return;
    }
    double* prev_basis = nurbs_basis_buffer + get_nurbs_basis_buffer_size(degree - 2);
    nurbs_basis_d1_buffer[0] = -degree / (knots[index + degree + 1] - knots[index + 1]) * prev_basis[0];
    for (int j = 1; j < degree; ++j) {
        int i = index + j;
        nurbs_basis_d1_buffer[j] = degree / (knots[i + degree] - knots[i]) * prev_basis[j - 1] -
            degree / (knots[i + degree + 1] - knots[i + 1]) * prev_basis[j];
    }
    int i = degree + index;
    nurbs_basis_d1_buffer[degree] = degree / (knots[i + degree] - knots[i]) * prev_basis[degree - 1];
}

inline void insert_knots(double* knot_buffer, WGVector2d* control_point_buffer, double* weight_buffer, 
    const double* knots, const WGVector2d* control_points, const double* weights,
    int degree, int control_point_count, int index, double t, int insert_count) {
    int knot_count = control_point_count + degree + 1;
    int knot_index = index + degree;
    memcpy(knot_buffer + knot_index + 1 + insert_count, knots + knot_index + 1, (knot_count - knot_index - 1) * sizeof(double));
    if (knot_buffer != knots) {
        memcpy(knot_buffer, knots, (knot_index + 1) * sizeof(double));
    }
    for (int i = 0; i < insert_count; ++i) {
        knot_buffer[knot_index + 1 + i] = t;
    }
    memcpy(control_point_buffer + index + 1 + insert_count, control_points + index + 1, (control_point_count - index - 1) * sizeof(WGVector2d));
    if (control_point_buffer != control_points) {
        memcpy(control_point_buffer, control_points, (index + 1) * sizeof(WGVector2d));
    }
    if (weights) {
        memcpy(weight_buffer + index + 1 + insert_count, weights + index + 1, (control_point_count - index - 1) * sizeof(double));
        if (weight_buffer != weights) {
            memcpy(weight_buffer, weights, (index + 1) * sizeof(double));
        }
    }
    if (weights) {
        for (int i = 0; i < insert_count; ++i) {
            int k = index + i;
            double prev_weight = weight_buffer[k];
            WGVector2d prev_point = control_point_buffer[k] * prev_weight;
            k = index + 1 + insert_count;
            int k1 = index + i + 1;
            int k2 = k1 + degree + insert_count - i;
            double b = knot_buffer[k2] - knot_buffer[k1];
            double a;
            if (b <= g_double_epsilon) {
                a = 0;
            }
            else {
                a = (t - knot_buffer[k1]) / b;
            }
            double curr_weight = weight_buffer[k];
            WGVector2d curr_point = control_point_buffer[k] * curr_weight;
            int k3 = index + i + 1;
            weight_buffer[k3] = curr_weight * a + prev_weight * (1 - a);
            control_point_buffer[k3] = (curr_point * a + prev_point * (1 - a)) / weight_buffer[k3];
            prev_weight = curr_weight;
            prev_point = curr_point;
            for (int j = index + i + 2; j <= knot_index + i; ++j) {
                int k = j + insert_count - i;
                int k1 = j;
                int k2 = k1 + degree + insert_count - i;
                double b = knot_buffer[k2] - knot_buffer[k1];
                double a;
                if (b <= g_double_epsilon) {
                    a = 0;
                }
                else {
                    a = (t - knot_buffer[k1]) / b;
                }
                double curr_weight = weight_buffer[k];
                WGVector2d curr_point = control_point_buffer[k] * curr_weight;
                int k3 = k - 1;
                weight_buffer[k3] = curr_weight * a + prev_weight * (1 - a);
                control_point_buffer[k3] = (curr_point * a + prev_point * (1 - a)) / weight_buffer[k3];
                prev_weight = curr_weight;
                prev_point = curr_point;
            }
        }
    }
    else {
        for (int i = 0; i < insert_count; ++i) {
            int k = index + i;
            WGVector2d prev_point = control_point_buffer[k];
            k = index + 1 + insert_count;
            int k1 = index + i + 1;
            int k2 = k1 + degree + insert_count - i;
            double b = knot_buffer[k2] - knot_buffer[k1];
            double a;
            if (b <= g_double_epsilon) {
                a = 0;
            }
            else {
                a = (t - knot_buffer[k1]) / b;
            }
            WGVector2d curr_point = control_point_buffer[k];
            control_point_buffer[index + i + 1] = curr_point * a + prev_point * (1 - a);
            prev_point = curr_point;
            for (int j = index + i + 2; j <= knot_index + i; ++j) {
                int k = j + insert_count - i;
                int k1 = j;
                int k2 = k1 + degree + insert_count - i;
                double b = knot_buffer[k2] - knot_buffer[k1];
                double a;
                if (b <= g_double_epsilon) {
                    a = 0;
                }
                else {
                    a = (t - knot_buffer[k1]) / b;
                }
                WGVector2d curr_point = control_point_buffer[k];
                control_point_buffer[k - 1] = curr_point * a + prev_point * (1 - a);
                prev_point = curr_point;
            }
        }
    }
}

inline void sub_bezier_min(int degree, double t, WGVector2d* control_points) {
    if (t == 1) {
        return;
    }
    if (t == 0) {
        for (int i = 1; i <= degree; ++i) {
            control_points[i] = control_points[0];
        }
    }
    else {
        for (int i = 0; i < degree; ++i) {
            for (int j = degree; j > i; --j) {
                control_points[j] = control_points[j - 1] * (1 - t) + control_points[j] * t;
            }
        }
    }
}

inline void sub_bezier_max(int degree, double t, WGVector2d* control_points) {
    if (t == 0) {
        return;
    }
    if (t == 1) {
        for (int i = 0; i < degree; ++i) {
            control_points[i] = control_points[degree];
        }
    }
    else {
        for (int i = degree; i > 0; --i) {
            for (int j = 0; j < i; ++j) {
                control_points[j] = control_points[j] * (1 - t) + control_points[j + 1] * t;
            }
        }
    }
}

inline void sub_bezier(int degree, const WSInterval& domain, WGVector2d* control_points) {
    assert(domain.Min >= 0 && domain.Max <= 1);
    if (domain.Min == domain.Max) {
        if (domain.Min == 0) {
            for (int i = 1; i <= degree; ++i) {
                control_points[i] = control_points[0];
            }
        }
        else if (domain.Min == 1) {
            for (int i = 0; i < degree; ++i) {
                control_points[i] = control_points[degree];
            }
        }
        else {
            WSReal t = domain.Min;
            sub_bezier_max(degree, t, control_points);
            for (int i = 1; i <= degree; ++i) {
                control_points[i] = control_points[0];
            }
        }
    }
    else {
        if (domain.Min > 0) {
            WSReal t = domain.Min;
            sub_bezier_max(degree, t, control_points);
        }
        if (domain.Max < 1) {
            WSReal t = (domain.Max - domain.Min) / (1 - domain.Min);
            sub_bezier_min(degree, t, control_points);
        }
    }
}

inline void sub_bezier_min(int degree, double t, double* weights) {
    if (t == 1) {
        return;
    }
    if (t == 0) {
        for (int i = 1; i <= degree; ++i) {
            weights[i] = weights[0];
        }
    }
    else {
        for (int i = 0; i < degree; ++i) {
            for (int j = degree; j > i; --j) {
                weights[j] = weights[j - 1] * (1 - t) + weights[j] * t;
            }
        }
    }
}

inline void sub_bezier_max(int degree, double t, double* weights) {
    if (t == 0) {
        return;
    }
    if (t == 1) {
        for (int i = 0; i < degree; ++i) {
            weights[i] = weights[degree];
        }
    }
    else {
        for (int i = degree; i > 0; --i) {
            for (int j = 0; j < i; ++j) {
                weights[j] = weights[j] * (1 - t) + weights[j + 1] * t;
            }
        }
    }
}

inline void sub_bezier(int degree, const WSInterval& domain, double* weights) {
    assert(domain.Min >= 0 && domain.Max <= 1);
    if (domain.Min == domain.Max) {
        if (domain.Min == 0) {
            for (int i = 1; i <= degree; ++i) {
                weights[i] = weights[0];
            }
        }
        else if (domain.Min == 1) {
            for (int i = 0; i < degree; ++i) {
                weights[i] = weights[degree];
            }
        }
        else {
            WSReal t = domain.Min;
            sub_bezier_max(degree, t, weights);
            for (int i = 1; i <= degree; ++i) {
                weights[i] = weights[0];
            }
        }
    }
    else {
        if (domain.Min > 0) {
            WSReal t = domain.Min;
            sub_bezier_max(degree, t, weights);
        }
        if (domain.Max < 1) {
            WSReal t = (domain.Max - domain.Min) / (1 - domain.Min);
            sub_bezier_min(degree, t, weights);
        }
    }
}

#endif
