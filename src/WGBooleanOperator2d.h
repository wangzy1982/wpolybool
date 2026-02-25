#ifndef _WG_BOOLEAN_OPERATOR_2D_
#define _WG_BOOLEAN_OPERATOR_2D_

#include "WGPolygon.h"
#include "WGTopoHelper2d.h"

/*
struct WGPolygonPolygonBooleanResult {
    WGPolygonPolygonBooleanResult(bool success, double distance_epsilon) :
        Success(success),
        DistanceEpsilon(distance_epsilon),
        VectorMoved(0, 0) {
    }

    WGPolygonPolygonBooleanResult(bool success, double distance_epsilon, const WGVector2d& vector_moved) :
        Success(success),
        DistanceEpsilon(distance_epsilon),
        VectorMoved(vector_moved) {
    }

    bool Success;
    double DistanceEpsilon;
    WGVector2d VectorMoved;
};

struct WGPolygonPolygonBooleanSetting {
    WGPolygonPolygonBooleanSetting() :
        SlightMoveEnable(true),
        DistanceEpsilon(1E-6) {
    }

    WGPolygonPolygonBooleanSetting(bool slight_move_enable, double distance_epsilon) :
        SlightMoveEnable(slight_move_enable),
        DistanceEpsilon(distance_epsilon) {
    }

    bool SlightMoveEnable;
    double DistanceEpsilon;
};

WGPolygonPolygonBooleanResult polygon_polygon_subtract(const WGPolygon* polygon0, const WGPolygon* polygon1,
    const WGPolygonPolygonBooleanSetting* setting, WGPolygon*& result_polygon);
*/

class WGPolygonPolygonBoolean {
public:
    struct Result {
        Result(bool success, bool is_singularity0, bool is_singularity1, const WGVector2d& recommended_moving_vector);
        bool Success;
        bool IsSingularity0;
        bool IsSingularity1;
        WGVector2d RecommendedMovingVector;
    };
public:
    static bool Intersect(const WGPolygon* polygon0, const WGPolygon* polygon1, double distance_epsilon, WGPolygon*& result_polygon);
    static bool Unite(const WGPolygon* polygon0, const WGPolygon* polygon1, double distance_epsilon, WGPolygon*& result_polygon);
    static bool Subtract(const WGPolygon* polygon0, const WGPolygon* polygon1, double distance_epsilon, WGPolygon*& result_polygon);
public:
    static bool Execute(const WGPolygon* polygon0, const WGPolygon* polygon1, double distance_epsilon, 
        WGIntersectHelper2d::IntersectCache* intersect_cache, bool adjust_distance_epsilon, double max_adjust_polygon_distance, 
        int result_count, WGTopoHelper2d::SegmentSelector** selectors, WGPolygon** result_polygons);
private:
    static Result Execute(const WGPolygon* polygon0, const WGPolygon* polygon1, double distance_epsilon, 
        WGIntersectHelper2d::IntersectCache* intersect_cache, int result_count, 
        WGTopoHelper2d::SegmentSelector** selectors, WGPolygon** result_polygons);
    static void ExecuteSeparated(WGTopoHelper2d::SplitterContainer* container, const WGPolygon* polygon0, const WGPolygon* polygon1,
        bool is_box_separated, double distance_epsilon, int result_count, WGTopoHelper2d::SegmentSelector** selectors, 
        bool& is_singularity0, bool& is_singularity1, WGPolygon** result_polygons);
};

#endif
