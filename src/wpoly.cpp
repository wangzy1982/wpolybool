#include "wpoly.h"
#include "WGPolygon.h"
#include "WGLine2d.h"
#include "WGArc2d.h"
#include "WGNurbsCurve2d.h"

WPOLY_API_C void* new_polygon() {
    return new WGPolygon();
}

WPOLY_API_C void free_polygon(void* polygon) {
    delete (WGPolygon*)polygon;
}

WPOLY_API_C void* new_loop() {
    return new WGWire2d();
}

WPOLY_API_C void add_loop(void* polygon, void* loop) {
    WGLoop2d* loop2 = new WGLoop2d();
    ((WGWire2d*)loop)->SetClosed(true);
    loop2->SetWire((WGWire2d*)loop, WGLoop2d::CalculateSignArea((WGWire2d*)loop));
    ((WGPolygon*)polygon)->AddLoop(loop2);
}

WPOLY_API_C void* new_line2d_edge(double start_x, double start_y, double end_x, double end_y) {
    return new WGEdge2d(new WGLine2d(WGVector2d(start_x, start_y), WGVector2d(end_x, end_y)));
}

WPOLY_API_C void* new_arc2d_edge(double center_x, double center_y, double radius, double start_angle, double delta_angle) {
    return new WGEdge2d(new WGArc2d(WGVector2d(center_x, center_y), radius, start_angle, delta_angle));
}

WPOLY_API_C void* new_nurbs2d_edge(int degree, int control_point_count, const double* knots, const double* control_points, const double* weights) {
    return new WGEdge2d(new WGNurbsCurve2d(degree, control_point_count, knots, (WGVector2d*)control_points, weights));
}

WPOLY_API_C void add_edge(void* loop, void* edge) {
    ((WGWire2d*)loop)->AddEdge((WGEdge2d*)edge);
}

WPOLY_API_C int get_loop_count(void* polygon) {
    return ((WGPolygon*)polygon)->GetLoopCount();
}

WPOLY_API_C void* get_loop(void* polygon, int loop_index) {
    return (void*)((WGPolygon*)polygon)->GetLoop(loop_index)->GetWire();
}

WPOLY_API_C int get_edge_count(void* loop) {
    return ((WGWire2d*)loop)->GetEdgeCount();
}

WPOLY_API_C void* get_edge(void* loop, int edge_index) {
    return (void*)((WGWire2d*)loop)->GetEdge(edge_index);
}

WPOLY_API_C bool is_line2d_edge(void* edge) {
    return ((WGEdge2d*)edge)->GetCurve()->GetType() == WGLine2d::Type::Instance();
}

WPOLY_API_C bool is_arc2d_edge(void* edge) {
    return ((WGEdge2d*)edge)->GetCurve()->GetType() == WGArc2d::Type::Instance();
}

WPOLY_API_C bool is_nurbs2d_edge(void* edge) {
    return ((WGEdge2d*)edge)->GetCurve()->GetType() == WGNurbsCurve2d::Type::Instance();
}

WPOLY_API_C void get_line2d_edge_data(void* edge, double& start_x, double& start_y, double& end_x, double& end_y) {
    WGLine2d* line = (WGLine2d*)((WGEdge2d*)edge)->GetCurve();
    start_x = line->GetStartPoint().X;
    start_y = line->GetStartPoint().Y;
    end_x = line->GetEndPoint().X;
    end_y = line->GetEndPoint().Y;
}

WPOLY_API_C void get_arc2d_edge_data(void* edge, double& center_x, double& center_y, double& radius, double& start_angle, double& delta_angle) {
    WGArc2d* arc = (WGArc2d*)((WGEdge2d*)edge)->GetCurve();
    center_x = arc->GetCenter().X;
    center_y = arc->GetCenter().Y;
    radius = arc->GetRadius();
    start_angle = arc->GetStartAngle();
    delta_angle = arc->GetDeltaAngle();
}

WPOLY_API_C void get_nurbs2d_edge_data(void* edge, int& degree, int& control_point_count, const double*& knots, const double*& control_points, const double*& weights) {
    WGNurbsCurve2d* nurbs = (WGNurbsCurve2d*)((WGEdge2d*)edge)->GetCurve();
    degree = nurbs->GetDegree();
    control_point_count = nurbs->GetControlPointCount();
    knots = nurbs->GetKnots();
    control_points = (const double*)nurbs->GetControlPoints();
    weights = nurbs->GetWeights();
}