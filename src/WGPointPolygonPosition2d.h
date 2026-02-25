#ifndef _WG_POINT_POLYGON_POSITION_2D_
#define _WG_POINT_POLYGON_POSITION_2D_

#include "WGPolygon.h"
#include "WGIntersectHelper2d.h"
#include "WGTopoHelper2d.h"

class WGPointPolygonPosition2d {
public:
    enum class Result {
        Inner,
        Outter,
        On,
        Unknown,
        Singularity
    };
public:
    static Result Execute(const WGVector2d& point, const WGPolygon* polygon, double distance_epsilon);
    static Result Execute(const WGVector2d& point, const WGPolygon* polygon, double distance_epsilon,
        WGIntersectHelper2d::IntersectCache* intersect_cache);
private:
    static bool QuickGetPointPolygonPosition(const WGVector2d& point, const WGPolygon* polygon, 
        double distance_epsilon, const WGWire2d* ray_wire, WGTopoHelper2d::SplitterContainer* container, 
        WGIntersectHelper2d::IntersectCache* intersect_cache, WGPointPolygonPosition2d::Result& result, int& confidence);
    static int SampleX(double base_x, double epsilon, int sample_count_per_piece,
        const WGTopoHelper2d::SplitPoint* split_point0, const WGTopoHelper2d::SplitPoint* split_point1);
};

#endif
