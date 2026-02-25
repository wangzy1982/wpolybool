#ifndef _WG_CURVE_INTERSECTER_2D_
#define _WG_CURVE_INTERSECTER_2D_

#include "WGCurve2d.h"
#include "WGIntersectHelper2d.h"

class WGPointCurveIntersection2d {
public:
    int PieceIndex;
    double T;
};

class WGCurveCurveIntersection2d {
public:
    int PieceIndex0;
    int PieceIndex1;
    double Ts0[2];
    double Ts1[2];
    bool IsSamples[2];
    bool IsOverlap;
};

class WGCurveIntersecter2d {
public:
    static void Intersect(const WGVector2d& point, const WGCurve2d* curve, int piece_index, int piece_count, 
        double distance_epsilon, WGIntersectHelper2d::IntersectCache* intersect_cache, bool& is_singularity,
        std::vector<WGPointCurveIntersection2d>& intersections);
public:
    static void Intersect(const WGCurve2d* curve0, const WGCurve2d* curve1, double distance_epsilon,
        WGIntersectHelper2d::IntersectCache* intersect_cache, bool& is_singularity0, bool& is_singularity1, 
        std::vector<WGCurveCurveIntersection2d>& intersections);
    static void Intersect(const WGCurve2d* curve0, int piece_index0, int piece_count0, 
        const WGCurve2d* curve1, int piece_index1, int piece_count1, double distance_epsilon, 
        WGIntersectHelper2d::IntersectCache* intersect_cache, bool& is_singularity0, bool& is_singularity1, 
        std::vector<WGCurveCurveIntersection2d>& intersections);
public:
    static void RemoveSingularity(const WGCurve2d* curve, double distance_epsilon, std::vector<WGCurve2d*>& curves);
private:
    static void AddIntersections(std::vector<WGPointCurveIntersection2d>& intersections,
        WGIntersectHelper2d::PointCurveIntersection* base_intersections, int base_intersection_count, int piece_index);
    static void AddIntersections(std::vector<WGPointCurveIntersection2d>& intersections,
        const std::vector<WGIntersectHelper2d::PointCurveIntersection>& base_intersections, const WSInterval& domain, int piece_index);
    static void AddIntersections(std::vector<WGCurveCurveIntersection2d>& intersections,
        WGIntersectHelper2d::CurveCurveIntersection* base_intersections, int base_intersection_count, 
        int piece_index0, int piece_index1, bool exchanged);
    static void AddIntersections(std::vector<WGCurveCurveIntersection2d>& intersections, 
        const std::vector<WGIntersectHelper2d::CurveCurveIntersection>& base_intersections,
        int piece_index0, int piece_index1, bool exchanged);
    static void AddIntersections(std::vector<WGCurveCurveIntersection2d>& intersections,
        const std::vector<WGIntersectHelper2d::CurveCurveIntersection>& base_intersections,
        const WSInterval* domain0, const WSInterval* domain1, int piece_index0, int piece_index1, bool exchanged);
    static void MergeDomains(std::vector<WSInterval>& domains);
};

#endif