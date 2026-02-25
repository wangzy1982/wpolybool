# README.md

### Introduction
This open-source library is a high-performance and highly robust tool library designed specifically for 2D polygon Boolean operations. It is fully compatible with complex polygons composed of curved edges such as <b>lines, arcs, and NURBS curves</b>, and can stably complete core Boolean operations including union, intersection, and difference. It solves the industry pain points of error-prone and unstable results in curved polygon operations in engineering scenarios.
<b>The core advantage of this algorithm library is its ultimate robustness</b>, so it can also be used as a fallback for other libraries, called when other libraries fail to perform operations.

### Core Features

#### 1. Support for Curved Edges
Supports lines, arcs, and Nurbs curves.
Curved edges participate in calculations in their native geometric form during the operation process, rather than being discretized into polyline segments, fundamentally ensuring the geometric accuracy of the operation results.

#### 2. Ultimate Robustness
Robustness is the core advantage of this library. In-depth optimizations have been made for the most prone boundary issues in Boolean operations, including:
- A stable curve intersection library with consistent topology at endpoints.
- A highly reliable positional relationship detection algorithm.
- Targeted processing of extreme cases such as near-coincidence and multiple complex coincidences.
- Various topological ambiguity analysis methods during the process to ensure correct results.
- Guided fine-tuning.

#### 3. Extremely Lightweight
Provides a concise API interface, easy to encapsulate into multiple languages such as C/C++/Python/C#, and ready to use out of the box.

### Usage Example

```cpp
int main() {
    //new first polygon
    void* poly0 = new_polygon();
    void* loop = new_loop();
    void* edge = new_line2d_edge(-1192.3076923076922, 642.01183431952666, -1136.0946745562126, 281.06508875739627);
    add_edge(loop, edge);
    edge = new_arc2d_edge(-1253.3499873368173, -151.97545645016470, 448.63450844639567, 1.3063649616864890, 3.3970043887824035);
    add_edge(loop, edge);
    double knots0[] = { 0, 0, 0, 0, 1, 2, 2, 2, 2 };
    double control_points0[] = {
        -1257.3964497041418, -600.59171597633133, -630.17751479289939, -544.37869822485197,
        -266.27218934911235, -44.378698224852073, -334.31952662721881, 639.05325443786967,
        -1192.3076923076922, 642.01183431952666 
    };
    edge = new_nurbs2d_edge(3, 5, knots0, control_points0, nullptr);
    add_edge(loop, edge);
    add_loop(poly0, loop);

    //new second polygon
    void* poly1 = new_polygon();
    loop = new_loop();
    double knots1[] = { 0, 0, 0, 0, 1, 1, 1, 1 };
    double control_points1[] = {
        -1079.1420118343194, 7.7662721893490527, -926.59023668639065, -491.49408284023662,
        2.5887573964497079, -588.57248520710050, -561.39053254437886, -98.557692307692378 
    };
    edge = new_nurbs2d_edge(3, 4, knots1, control_points1, nullptr);
    add_edge(loop, edge);
    edge = new_arc2d_edge(-338.33989851885724, 160.37112728687720, 341.75388652629005, 4.0012930208426623, -2.5901072848103430);
    add_edge(loop, edge);
    edge = new_arc2d_edge(-601.09874085543993, 449.86098636939960, 320.67574801580520, 0.14999651125525995, 2.7534215162831286);
    add_edge(loop, edge);
    edge = new_line2d_edge(-912.72189349112432, 525.51775147929004, -1079.1420118343194, 7.7662721893490527);
    add_edge(loop, edge);
    add_loop(poly1, loop);

    //intersect
    void* poly = intersect(poly0, poly1, 1E-6);

    //print result
    int loop_count = get_loop_count(poly);
    for (int i = 0; i < loop_count; ++i) {
        std::cout << "Loop - " << i << std::endl;
        loop = get_loop(poly, i);
        int edge_count = get_edge_count(loop);
        for (int j = 0; j < edge_count; ++j) {
            edge = get_edge(loop, j);
            std::cout << "\tEdge - " << j << ": ";
            if (is_line2d_edge(edge)) {
                double start_x;
                double start_y;
                double end_x;
                double end_y;
                get_line2d_edge_data(edge, start_x, start_y, end_x, end_y);
                std::cout << "Line" << std::endl;
                std::cout << "\t\tStart point: " << "(" << start_x << "," << start_y << ")" << std::endl;
                std::cout << "\t\tEnd point: " << "(" << end_x << "," << end_y << ")" << std::endl;
            }
            else if (is_arc2d_edge(edge)) {
                double center_x;
                double center_y;
                double radius;
                double start_angle;
                double delta_angle;
                get_arc2d_edge_data(edge, center_x, center_y, radius, start_angle, delta_angle);
                std::cout << "Arc" << std::endl;
                std::cout << "\t\tCenter: " << "(" << center_x << "," << center_y << ")" << std::endl;
                std::cout << "\t\tRadius: " << radius << std::endl;
                std::cout << "\t\tStart angle: " << start_angle << std::endl;
                std::cout << "\t\tDelta angle: " << delta_angle << std::endl;
            }
            else if (is_nurbs2d_edge(edge)) {
                int degree;
                int control_point_count;
                const double* knots;
                const double* control_points;
                const double* weights;
                get_nurbs2d_edge_data(edge, degree, control_point_count, knots, control_points, weights);
                std::cout << "Nurbs" << std::endl;
                std::cout << "\t\tDegree: " << degree << std::endl;
                std::cout << "\t\tControl point count: " << control_point_count << std::endl;
                std::cout << "\t\tKnots: [ ";
                for (int k = 0; k < control_point_count + degree + 1; ++k) {
                    if (k == 0) {
                        std::cout << knots[k];
                    }
                    else {
                        std::cout << " ," << knots[k];
                    }
                }
                std::cout << " ]" << std::endl;
                std::cout << "\t\tControl points: [ ";
                for (int k = 0; k < control_point_count * 2; k += 2) {
                    if (k == 0) {
                        std::cout << "(" << control_points[k] << "," << control_points[k + 1] << ")";
                    }
                    else {
                        std::cout << " ," << "(" << control_points[k] << "," << control_points[k + 1] << ")";
                    }
                }
                std::cout << " ]" << std::endl;
                if (weights) {
                    std::cout << "\t\tWeights: [ ";
                    for (int k = 0; k < control_point_count; ++k) {
                        if (k == 0) {
                            std::cout << weights[k];
                        }
                        else {
                            std::cout << " ," << weights[k];
                        }
                    }
                    std::cout << " ]" << std::endl;
                }
            }
        }
    }

    //free
    free(poly);
    free_polygon(poly0);
    free_polygon(poly1);

    return 0;
}
```

### About the Professional Version
We also provide corresponding professional version services. The professional version uses a more efficient curve intersection algorithm library (with about 5 times improvement in Nurbs curve intersection). If you have strong requirements for the efficiency of Boolean operations involving Nurbs curves, please feel free to contact us.

### Contact Information
1179422870@qq.com