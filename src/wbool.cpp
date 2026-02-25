#include "wbool.h"
#include "WGBooleanOperator2d.h"

WPOLY_API_C void* intersect(void* polygon0, void* polygon1, double distance_epsilon) {
    WGPolygon* result = nullptr;
    if (WGPolygonPolygonBoolean::Intersect((WGPolygon*)polygon0, (WGPolygon*)polygon1, distance_epsilon, result)) {
        return result;
    }
    return nullptr;
}

WPOLY_API_C void* unite(void* polygon0, void* polygon1, double distance_epsilon) {
    WGPolygon* result = nullptr;
    if (WGPolygonPolygonBoolean::Unite((WGPolygon*)polygon0, (WGPolygon*)polygon1, distance_epsilon, result)) {
        return result;
    }
    return nullptr;
}

WPOLY_API_C void* subtract(void* polygon0, void* polygon1, double distance_epsilon) {
    WGPolygon* result = nullptr;
    if (WGPolygonPolygonBoolean::Subtract((WGPolygon*)polygon0, (WGPolygon*)polygon1, distance_epsilon, result)) {
        return result;
    }
    return nullptr;
}