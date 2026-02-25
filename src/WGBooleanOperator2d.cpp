#include "WGBooleanOperator2d.h"
#include <assert.h>

/*
void polygon_polygon_boolean_separated(WGTopoHelper2d::SplitterContainer* container, const WGPolygon* polygon0, const WGPolygon* polygon1,
    bool is_box_separated, double distance_epsilon, int result_count, WGTopoHelper2d::SegmentSelector** selectors, WGPolygon** result_polygons) {
    for (int result_index = 0; result_index < result_count; ++result_index) {
        result_polygons[result_index] = new WGPolygon();
        WGTopoHelper2d::AddPolygonSeparatedLoops(container, polygon0, polygon1, is_box_separated, 
            selectors[result_index], distance_epsilon, result_polygons[result_index]);
    }
}

WGPolygonPolygonBooleanResult polygon_polygon_boolean(const WGPolygon* polygon0, const WGPolygon* polygon1, WGVector2d* vector_moved,
    const WGPolygonPolygonBooleanSetting* setting, int result_count, WGTopoHelper2d::SegmentSelector** selectors, WGPolygon** result_polygons) {
    WGTopoHelper2d::SplitterContainer container;
    std::vector<WGCurveCurveIntersection2d> intersections;
    intersections.reserve(16);
    std::vector<WGTopoHelper2d::Splitter> splitters;
    splitters.reserve(16);
    for (int loop_index0 = 0; loop_index0 < polygon0->GetLoopCount(); ++loop_index0) {
        const WGLoop2d* loop0 = polygon0->GetLoop(loop_index0);
        const WGWire2d* wire0 = loop0->GetWire();
        for (int edge_index0 = 0; edge_index0 < wire0->GetEdgeCount(); ++edge_index0) {
            const WGEdge2d* edge0 = wire0->GetEdge(edge_index0);
            const WGCurve2d* curve0 = edge0->GetCurve();
            for (int loop_index1 = 0; loop_index1 < polygon1->GetLoopCount(); ++loop_index1) {
                const WGLoop2d* loop1 = polygon1->GetLoop(loop_index1);
                const WGWire2d* wire1 = loop1->GetWire();
                for (int edge_index1 = 0; edge_index1 < wire1->GetEdgeCount(); ++edge_index1) {
                    const WGEdge2d* edge1 = wire1->GetEdge(edge_index1);
                    const WGCurve2d* curve1 = edge1->GetCurve();
                    curve_curve_intersect(curve0, curve1, setting->DistanceEpsilon, true, intersections);
                    for (auto intersection_itr = intersections.begin(); intersection_itr != intersections.end(); ++intersection_itr) {
                        WGCurveCurveIntersection2d& intersection = *intersection_itr;
                        splitters.push_back(WGTopoHelper2d::NewSplitter(&container, wire0, edge_index0, wire1, edge_index1, &intersection));
                    }
                    intersections.clear();
                }
            }
        }
    }
    if (splitters.size() == 0) {
        polygon_polygon_boolean_separated(&container, polygon0, polygon1, false, setting->DistanceEpsilon, result_count, selectors, result_polygons);
        return WGPolygonPolygonBooleanResult(true, setting->DistanceEpsilon, vector_moved ? *vector_moved : WGVector2d(0, 0));
    }
    WGTopoHelper2d::AddSplitters(&container, splitters);
    WGTopoHelper2d::SetSplittersIsPoint(&container, setting->DistanceEpsilon);
    WGTopoHelper2d::AddCoincideSideSplitter(&container);
    for (int iterate_count = 0; iterate_count <= 2; ++iterate_count) {
        int splitter_count = (int)container.Splitters.size();
        for (int index = 0; index <= 1; ++index) {
            WGTopoHelper2d::SplitCoincideSplitterError* error1 = WGTopoHelper2d::SplitCoincideSplitters(&container, index, 
                setting->DistanceEpsilon, iterate_count == 2);
            if (error1) {
                if (setting->SlightMoveEnable) {
                    WGVector2d vt = WGTopoHelper2d::GetSlightMoveDirectionByError(&container, index, error1) * (setting->DistanceEpsilon * 2);
                    delete error1;
                    WGPolygonPolygonBooleanSetting current_setting = *setting;
                    WGPolygon* new_polygon1;
                    if (!vector_moved) {
                        new_polygon1 = polygon1->Clone();
                    }
                    else {
                        new_polygon1 = (WGPolygon*)polygon1;
                    }
                    new_polygon1->Move(vt);
                    if (vector_moved) {
                        vt = *vector_moved + vt;
                    }
                    WGPolygonPolygonBooleanResult result = polygon_polygon_boolean(
                        polygon0, new_polygon1, vector_moved, &current_setting, result_count, selectors, result_polygons);
                    if (new_polygon1 != polygon1) {
                        delete new_polygon1;
                    }
                    return result;
                }
                else {
                    delete error1;
                }
                return WGPolygonPolygonBooleanResult(false, setting->DistanceEpsilon, vector_moved ? *vector_moved : WGVector2d(0, 0));
            }
        }
        if (splitter_count == (int)container.Splitters.size()) {
            break;
        }
    }
    WGTopoHelper2d::BuildMergedSplitters(&container);
    WGTopoHelper2d::BuildMergedSplitterPositions(&container, true, true);
    WGTopoHelper2d::InitializeIntersections(&container);
    WGTopoHelper2d::MergeIntersections(&container, 0);
    WGTopoHelper2d::MergeIntersections(&container, 1);
    WGTopoHelper2d::MergeIntersectionParts(&container, 0);
    WGTopoHelper2d::MergeIntersectionParts(&container, 1);
    WGTopoHelper2d::BuildSortedIntersectionPartIndices(&container, 0);
    WGTopoHelper2d::BuildSortedIntersectionPartIndices(&container, 1);
    WGTopoHelper2d::BuildSegments(&container, 0);
    WGTopoHelper2d::BuildSegments(&container, 1);
    WGTopoHelper2d::CalculateSegmentsPosition(&container, 0, polygon1, setting->DistanceEpsilon);
    WGTopoHelper2d::CalculateSegmentsPosition(&container, 1, polygon0, setting->DistanceEpsilon);
    for (int index = 0; index <= 1; ++index) {
        int error_segment_index = -1;
        for (int i = 0; i < (int)container.Segments[index].size(); ++i) {
            if (!WGTopoHelper2d::CheckSegment(&container, index, container.Segments[index].at(i))) {
                error_segment_index = i;
                break;
            }
        }
        if (error_segment_index != -1) {
            if (setting->SlightMoveEnable) {
                WGVector2d vt = WGTopoHelper2d::GetSlightMoveDirectionByErrorSegmentIndex(&container, index, error_segment_index) * (setting->DistanceEpsilon * 2);
                WGPolygonPolygonBooleanSetting current_setting = *setting;
                WGPolygon* new_polygon1;
                if (!vector_moved) {
                    new_polygon1 = polygon1->Clone();
                }
                else {
                    new_polygon1 = (WGPolygon*)polygon1;
                }
                new_polygon1->Move(vt);
                if (vector_moved) {
                    vt = *vector_moved + vt;
                }
                WGPolygonPolygonBooleanResult result = polygon_polygon_boolean(
                    polygon0, new_polygon1, vector_moved, &current_setting, result_count, selectors, result_polygons);
                if (new_polygon1 != polygon1) {
                    delete new_polygon1;
                }
                return result;
            }
            return WGPolygonPolygonBooleanResult(false, setting->DistanceEpsilon, vector_moved ? *vector_moved : WGVector2d(0, 0));
        }
    }
    WGTopoHelper2d::BuildCoincideSegments(&container);
    int error_coincide_segment_index = -1;
    for (int i = 0; i < (int)container.CoincideSegments.size(); ++i) {
        if (!WGTopoHelper2d::CheckCoincideSegmentDirectionFlags(container.CoincideSegments.at(i))) {
            error_coincide_segment_index = i;
            break;
        }
    }
    if (error_coincide_segment_index != -1) {
        if (setting->SlightMoveEnable) {
            WGVector2d vt = WGTopoHelper2d::GetSlightMoveDirectionByErrorCoincideSegmentIndex(
                &container, error_coincide_segment_index) * (setting->DistanceEpsilon * 2);
            WGPolygonPolygonBooleanSetting current_setting = *setting;
            WGPolygon* new_polygon1;
            if (!vector_moved) {
                new_polygon1 = polygon1->Clone();
            }
            else {
                new_polygon1 = (WGPolygon*)polygon1;
            }
            new_polygon1->Move(vt);
            if (vector_moved) {
                vt = *vector_moved + vt;
            }
            WGPolygonPolygonBooleanResult result = polygon_polygon_boolean(
                polygon0, new_polygon1, vector_moved, &current_setting, result_count, selectors, result_polygons);
            if (new_polygon1 != polygon1) {
                delete new_polygon1;
            }
            return result;
        }
        return WGPolygonPolygonBooleanResult(false, setting->DistanceEpsilon, vector_moved ? *vector_moved : WGVector2d(0, 0));
    }
    WGTopoHelper2d::BuildCoincideSegmentFlags(&container, polygon0, polygon1, setting->DistanceEpsilon);
    WGTopoHelper2d::BuildSegmentsCoincideSegmentIndex(&container, 0);
    WGTopoHelper2d::BuildSegmentsCoincideSegmentIndex(&container, 1);
    for (int index = 0; index <= 1; ++index) {
        int error_segment_index = -1;
        for (int i = 0; i < (int)container.Segments[index].size(); ++i) {
            if (!WGTopoHelper2d::CheckCoincideSegmentIndex(container.Segments[index].at(i))) {
                error_segment_index = i;
                break;
            }
        }
        if (error_segment_index != -1) {
            if (setting->SlightMoveEnable) {
                WGVector2d vt = WGTopoHelper2d::GetSlightMoveDirectionByErrorSegmentIndex(&container, index, error_segment_index) * (setting->DistanceEpsilon * 2);
                WGPolygonPolygonBooleanSetting current_setting = *setting;
                WGPolygon* new_polygon1;
                if (!vector_moved) {
                    new_polygon1 = polygon1->Clone();
                }
                else {
                    new_polygon1 = (WGPolygon*)polygon1;
                }
                new_polygon1->Move(vt);
                if (vector_moved) {
                    vt = *vector_moved + vt;
                }
                WGPolygonPolygonBooleanResult result = polygon_polygon_boolean(
                    polygon0, new_polygon1, vector_moved, &current_setting, result_count, selectors, result_polygons);
                if (new_polygon1 != polygon1) {
                    delete new_polygon1;
                }
                return result;
            }
            return WGPolygonPolygonBooleanResult(false, setting->DistanceEpsilon, vector_moved ? *vector_moved : WGVector2d(0, 0));
        }
    }
    bool success = true;
    for (int result_index = 0; result_index < result_count; ++result_index) {
        selectors[result_index]->SelectSegments(&container);
        std::vector<WGTopoHelper2d::Loop> loops;
        WGTopoHelper2d::SearchLoops(&container, loops);
        int error_loop_index = -1;
        for (int i = 0; i < (int)loops.size(); ++i) {
            if (!WGTopoHelper2d::CheckLoop(loops.at(i))) {
                error_loop_index = i;
                break;
            }
        }
        if (error_loop_index != -1) {
            for (int i = 0; i < result_index; ++i) {
                delete result_polygons[i];
            }
            success = false;
            break;
        }
        WGTopoHelper2d::RepairEdgeStartEnds(&container);
        result_polygons[result_index] = new WGPolygon();
        WGTopoHelper2d::AddPolygonLoops(&container, loops, result_polygons[result_index]);
        WGTopoHelper2d::AddPolygonSeparatedLoops(&container, polygon0, polygon1, false, selectors[result_index], 
            setting->DistanceEpsilon, result_polygons[result_index]);
    }
    if (!success) {
        error_coincide_segment_index = -1;
        for (int i = 0; i < (int)container.CoincideSegments.size(); ++i) {
            if (!WGTopoHelper2d::CheckCoincideSegmentFlag(container.CoincideSegments.at(i))) {
                error_coincide_segment_index = i;
                break;
            }
        }
        if (error_coincide_segment_index != -1) {
            if (setting->SlightMoveEnable) {
                WGVector2d vt = WGTopoHelper2d::GetSlightMoveDirectionByErrorCoincideSegmentIndex(
                    &container, error_coincide_segment_index) * (setting->DistanceEpsilon * 2);
                WGPolygonPolygonBooleanSetting current_setting = *setting;
                WGPolygon* new_polygon1;
                if (!vector_moved) {
                    new_polygon1 = polygon1->Clone();
                }
                else {
                    new_polygon1 = (WGPolygon*)polygon1;
                }
                new_polygon1->Move(vt);
                if (vector_moved) {
                    vt = *vector_moved + vt;
                }
                WGPolygonPolygonBooleanResult result = polygon_polygon_boolean(
                    polygon0, new_polygon1, vector_moved, &current_setting, result_count, selectors, result_polygons);
                if (new_polygon1 != polygon1) {
                    delete new_polygon1;
                }
                return result;
            }
            return WGPolygonPolygonBooleanResult(false, setting->DistanceEpsilon, vector_moved ? *vector_moved : WGVector2d(0, 0));
        }
        return WGPolygonPolygonBooleanResult(false, setting->DistanceEpsilon, vector_moved ? *vector_moved : WGVector2d(0, 0));
    }
    return WGPolygonPolygonBooleanResult(true, setting->DistanceEpsilon, vector_moved ? *vector_moved : WGVector2d(0, 0));
}

WGPolygonPolygonBooleanResult polygon_polygon_subtract(const WGPolygon* polygon0, const WGPolygon* polygon1,
    const WGPolygonPolygonBooleanSetting* setting, WGPolygon*& result_polygon) {
    WGPolygonPolygonBooleanSetting current_setting = *setting;
    double joint_epsilon = polygon0->GetJointEpsilon();
    if (joint_epsilon > current_setting.DistanceEpsilon) {
        current_setting.DistanceEpsilon = joint_epsilon;
    }
    joint_epsilon = polygon1->GetJointEpsilon();
    if (joint_epsilon > current_setting.DistanceEpsilon) {
        current_setting.DistanceEpsilon = joint_epsilon;
    }
    WGTopoHelper2d::PolygonPolygonSubtractionSegmentSelector selector;
    WGTopoHelper2d::SegmentSelector* selectors[1] = { &selector };
    if (!polygon0->GetBox().IsIntersected(polygon1->GetBox(), current_setting.DistanceEpsilon)) {
        polygon_polygon_boolean_separated(nullptr, polygon0, polygon1, true, current_setting.DistanceEpsilon, 1, selectors, &result_polygon);
        return WGPolygonPolygonBooleanResult(true, current_setting.DistanceEpsilon, WGVector2d(0, 0));
    }
    return polygon_polygon_boolean(polygon0, polygon1, nullptr, &current_setting, 1, selectors, &result_polygon);
}

*/

WGPolygonPolygonBoolean::Result::Result(bool success, bool is_singularity0, bool is_singularity1, const WGVector2d& recommended_moving_vector) :
    Success(success),
    IsSingularity0(is_singularity0),
    IsSingularity1(is_singularity1),
    RecommendedMovingVector(recommended_moving_vector) {
}

bool WGPolygonPolygonBoolean::Intersect(const WGPolygon* polygon0, const WGPolygon* polygon1, double distance_epsilon, WGPolygon*& result_polygon) {
    PolygonPolygonIntersectionSegmentSelector selector;
    WGTopoHelper2d::SegmentSelector* selectors[1] = { &selector };
    WGIntersectHelper2d::IntersectCache intersect_cache;
    return Execute(polygon0, polygon1, distance_epsilon, &intersect_cache, true, distance_epsilon * 10, 1, selectors, &result_polygon);
}

bool WGPolygonPolygonBoolean::Unite(const WGPolygon* polygon0, const WGPolygon* polygon1, double distance_epsilon, WGPolygon*& result_polygon) {
    PolygonPolygonUnionSegmentSelector selector;
    WGTopoHelper2d::SegmentSelector* selectors[1] = { &selector };
    WGIntersectHelper2d::IntersectCache intersect_cache;
    return Execute(polygon0, polygon1, distance_epsilon, &intersect_cache, true, distance_epsilon * 10, 1, selectors, &result_polygon);
}

bool WGPolygonPolygonBoolean::Subtract(const WGPolygon* polygon0, const WGPolygon* polygon1, double distance_epsilon, WGPolygon*& result_polygon) {
    PolygonPolygonSubtractionSegmentSelector selector;
    WGTopoHelper2d::SegmentSelector* selectors[1] = { &selector };
    WGIntersectHelper2d::IntersectCache intersect_cache;
    return Execute(polygon0, polygon1, distance_epsilon, &intersect_cache, true, distance_epsilon * 10, 1, selectors, &result_polygon);
}

bool WGPolygonPolygonBoolean::Execute(const WGPolygon* polygon0, const WGPolygon* polygon1, double distance_epsilon,
    WGIntersectHelper2d::IntersectCache* intersect_cache, bool adjust_distance_epsilon, double max_adjust_polygon_distance, 
    int result_count, WGTopoHelper2d::SegmentSelector** selectors, WGPolygon** result_polygons) {
    if (adjust_distance_epsilon) {
        double joint_epsilon = polygon0->GetJointEpsilon();
        if (joint_epsilon > distance_epsilon) {
            distance_epsilon = joint_epsilon;
        }
        joint_epsilon = polygon1->GetJointEpsilon();
        if (joint_epsilon > distance_epsilon) {
            distance_epsilon = joint_epsilon;
        }
    }
    WGPolygon* current_polygon0 = (WGPolygon*)polygon0;
    WGPolygon* current_polygon1 = (WGPolygon*)polygon1;
    WGVector2d adjust_polygon_vt(0, 0);
    while (true) {
        Result r = Execute(current_polygon0, current_polygon1, distance_epsilon, intersect_cache, result_count, selectors, result_polygons);
        if (r.Success) {
            if (current_polygon0 != polygon0) {
                delete current_polygon0;
            }
            if (current_polygon1 != polygon1) {
                delete current_polygon1;
            }
            return true;
        }
        bool b = false;
        if (r.IsSingularity0) {
            if (current_polygon0 == polygon0) {
                current_polygon0 = polygon0->Clone();
            }
            if (!current_polygon0->RepairSingularity(distance_epsilon * 0.1)) {
                if (current_polygon0 != polygon0) {
                    delete current_polygon0;
                }
                if (current_polygon1 != polygon1) {
                    delete current_polygon1;
                }
                return false;
            }
            b = true;
        }
        if (r.IsSingularity1) {
            if (current_polygon1 == polygon1) {
                current_polygon1 = polygon1->Clone();
            }
            if (!current_polygon1->RepairSingularity(distance_epsilon * 0.1)) {
                if (current_polygon0 != polygon0) {
                    delete current_polygon0;
                }
                if (current_polygon1 != polygon1) {
                    delete current_polygon1;
                }
                return false;
            }
            b = true;
        }
        if (b) {
            continue;
        }
        double d = r.RecommendedMovingVector.SqrLength();
        if (d == 0) {
            if (current_polygon0 != polygon0) {
                delete current_polygon0;
            }
            if (current_polygon1 != polygon1) {
                delete current_polygon1;
            }
            return false;
        }
        adjust_polygon_vt = adjust_polygon_vt + r.RecommendedMovingVector;
        d = adjust_polygon_vt.SqrLength();
        if (d > max_adjust_polygon_distance * max_adjust_polygon_distance) {
            if (current_polygon0 != polygon0) {
                delete current_polygon0;
            }
            if (current_polygon1 != polygon1) {
                delete current_polygon1;
            }
            return false;
        }
        if (current_polygon1 == polygon1) {
            current_polygon1 = polygon1->Clone();
        }
        current_polygon1->Move(r.RecommendedMovingVector);
    }
}

WGPolygonPolygonBoolean::Result WGPolygonPolygonBoolean::Execute(const WGPolygon* polygon0, const WGPolygon* polygon1, double distance_epsilon,
    WGIntersectHelper2d::IntersectCache* intersect_cache, int result_count,
    WGTopoHelper2d::SegmentSelector** selectors, WGPolygon** result_polygons) {
    WGTopoHelper2d::SplitterContainer container;
    std::vector<WGCurveCurveIntersection2d> intersections;
    intersections.reserve(16);
    std::vector<WGTopoHelper2d::Splitter> splitters;
    splitters.reserve(16);
    for (int loop_index0 = 0; loop_index0 < polygon0->GetLoopCount(); ++loop_index0) {
        const WGLoop2d* loop0 = polygon0->GetLoop(loop_index0);
        const WGWire2d* wire0 = loop0->GetWire();
        for (int edge_index0 = 0; edge_index0 < wire0->GetEdgeCount(); ++edge_index0) {
            const WGEdge2d* edge0 = wire0->GetEdge(edge_index0);
            const WGCurve2d* curve0 = edge0->GetCurve();
            for (int loop_index1 = 0; loop_index1 < polygon1->GetLoopCount(); ++loop_index1) {
                const WGLoop2d* loop1 = polygon1->GetLoop(loop_index1);
                const WGWire2d* wire1 = loop1->GetWire();
                for (int edge_index1 = 0; edge_index1 < wire1->GetEdgeCount(); ++edge_index1) {
                    const WGEdge2d* edge1 = wire1->GetEdge(edge_index1);
                    const WGCurve2d* curve1 = edge1->GetCurve();
                    bool is_singularity0;
                    bool is_singularity1;
                    WGCurveIntersecter2d::Intersect(curve0, curve1, distance_epsilon, intersect_cache, is_singularity0, is_singularity1, intersections);
                    if (is_singularity0 || is_singularity1) {
                        return Result(false, is_singularity0, is_singularity1, WGVector2d(0, 0));
                    }
                    for (auto intersection_itr = intersections.begin(); intersection_itr != intersections.end(); ++intersection_itr) {
                        WGCurveCurveIntersection2d& intersection = *intersection_itr;
                        splitters.push_back(WGTopoHelper2d::NewSplitter(&container, wire0, edge_index0, wire1, edge_index1, &intersection));
                    }
                    intersections.clear();
                }
            }
        }
    }
    if (splitters.size() == 0) {
        bool is_singularity0;
        bool is_singularity1;
        ExecuteSeparated(&container, polygon0, polygon1, false, distance_epsilon, result_count, selectors, 
            is_singularity0, is_singularity1, result_polygons);
        if (is_singularity0 || is_singularity1) {
            return Result(false, is_singularity0, is_singularity1, WGVector2d(0, 0));
        }
        return Result(true, false, false, WGVector2d(0, 0));
    }
    WGTopoHelper2d::AddSplitters(&container, splitters);
    WGTopoHelper2d::SetSplittersIsPoint(&container, distance_epsilon);
    WGTopoHelper2d::AddCoincideSideSplitter(&container);
    for (int iterate_count = 0; iterate_count <= 2; ++iterate_count) {
        int splitter_count = (int)container.Splitters.size();
        for (int index = 0; index <= 1; ++index) {
            bool is_singularity;
            WGTopoHelper2d::SplitCoincideSplitterError* error1 = WGTopoHelper2d::SplitCoincideSplitters(&container, index,
                distance_epsilon, intersect_cache, is_singularity, iterate_count == 2);
            if (is_singularity) {
                return Result(false, index == 1, index == 0, WGVector2d(0, 0));
            }
            if (error1) {
                WGVector2d vt = WGTopoHelper2d::GetSlightMoveDirectionByError(&container, index, error1) * (distance_epsilon * 2);
                delete error1;
                return Result(false, false, false, vt);
            }
        }
        if (splitter_count == (int)container.Splitters.size()) {
            break;
        }
    }
    WGTopoHelper2d::BuildMergedSplitters(&container);
    WGTopoHelper2d::BuildMergedSplitterPositions(&container, true, true);
    WGTopoHelper2d::InitializeIntersections(&container);
    WGTopoHelper2d::MergeIntersections(&container, 0);
    WGTopoHelper2d::MergeIntersections(&container, 1);
    WGTopoHelper2d::MergeIntersectionParts(&container, 0);
    WGTopoHelper2d::MergeIntersectionParts(&container, 1);
    WGTopoHelper2d::BuildSortedIntersectionPartIndices(&container, 0);
    WGTopoHelper2d::BuildSortedIntersectionPartIndices(&container, 1);
    WGTopoHelper2d::BuildSegments(&container, 0);
    WGTopoHelper2d::BuildSegments(&container, 1);
    bool is_singularity;
    WGTopoHelper2d::CalculateSegmentsPosition(&container, 0, polygon1, distance_epsilon, intersect_cache, is_singularity);
    if (is_singularity) {
        return Result(false, false, true, WGVector2d(0, 0));
    }
    WGTopoHelper2d::CalculateSegmentsPosition(&container, 1, polygon0, distance_epsilon, intersect_cache, is_singularity);
    if (is_singularity) {
        return Result(false, true, false, WGVector2d(0, 0));
    }
    for (int index = 0; index <= 1; ++index) {
        int error_segment_index = -1;
        for (int i = 0; i < (int)container.Segments[index].size(); ++i) {
            if (!WGTopoHelper2d::CheckSegment(&container, index, container.Segments[index].at(i))) {
                error_segment_index = i;
                break;
            }
        }
        if (error_segment_index != -1) {
            WGVector2d vt = WGTopoHelper2d::GetSlightMoveDirectionByErrorSegmentIndex(&container, index, error_segment_index) * (distance_epsilon * 2);
            return Result(false, false, false, vt);
        }
    }
    WGTopoHelper2d::BuildCoincideSegments(&container);
    int error_coincide_segment_index = -1;
    for (int i = 0; i < (int)container.CoincideSegments.size(); ++i) {
        if (!WGTopoHelper2d::CheckCoincideSegmentDirectionFlags(container.CoincideSegments.at(i))) {
            error_coincide_segment_index = i;
            break;
        }
    }
    if (error_coincide_segment_index != -1) {
        WGVector2d vt = WGTopoHelper2d::GetSlightMoveDirectionByErrorCoincideSegmentIndex(
            &container, error_coincide_segment_index) * (distance_epsilon * 2);
        return Result(false, false, false, vt);
    }
    bool is_singularity0;
    bool is_singularity1;
    WGTopoHelper2d::BuildCoincideSegmentFlags(&container, polygon0, polygon1, distance_epsilon, intersect_cache, is_singularity0, is_singularity1);
    if (is_singularity0 || is_singularity1) {
        return Result(false, is_singularity0, is_singularity1, WGVector2d(0, 0));
    }
    WGTopoHelper2d::BuildSegmentsCoincideSegmentIndex(&container, 0);
    WGTopoHelper2d::BuildSegmentsCoincideSegmentIndex(&container, 1);
    for (int index = 0; index <= 1; ++index) {
        int error_segment_index = -1;
        for (int i = 0; i < (int)container.Segments[index].size(); ++i) {
            if (!WGTopoHelper2d::CheckCoincideSegmentIndex(container.Segments[index].at(i))) {
                error_segment_index = i;
                break;
            }
        }
        if (error_segment_index != -1) {
            WGVector2d vt = WGTopoHelper2d::GetSlightMoveDirectionByErrorSegmentIndex(&container, index, error_segment_index) * (distance_epsilon * 2);
            return Result(false, false, false, vt);
        }
    }
    bool success = true;
    for (int result_index = 0; result_index < result_count; ++result_index) {
        selectors[result_index]->SelectSegments(&container);
        std::vector<WGTopoHelper2d::Loop> loops;
        WGTopoHelper2d::SearchLoops(&container, loops);
        int error_loop_index = -1;
        for (int i = 0; i < (int)loops.size(); ++i) {
            if (!WGTopoHelper2d::CheckLoop(loops.at(i))) {
                error_loop_index = i;
                break;
            }
        }
        if (error_loop_index != -1) {
            for (int i = 0; i < result_index; ++i) {
                delete result_polygons[i];
            }
            success = false;
            break;
        }
        WGTopoHelper2d::RepairEdgeStartEnds(&container);
        result_polygons[result_index] = new WGPolygon();
        WGTopoHelper2d::AddPolygonLoops(&container, loops, result_polygons[result_index]);
        bool is_singularity0;
        bool is_singularity1;
        WGTopoHelper2d::AddPolygonSeparatedLoops(&container, polygon0, polygon1, false, selectors[result_index],
            distance_epsilon, is_singularity0, is_singularity1, result_polygons[result_index]);
        if (is_singularity0 || is_singularity1) {
            return Result(false, is_singularity0, is_singularity1, WGVector2d(0, 0));
        }
    }
    if (!success) {
        error_coincide_segment_index = -1;
        for (int i = 0; i < (int)container.CoincideSegments.size(); ++i) {
            if (!WGTopoHelper2d::CheckCoincideSegmentFlag(container.CoincideSegments.at(i))) {
                error_coincide_segment_index = i;
                break;
            }
        }
        if (error_coincide_segment_index != -1) {
            WGVector2d vt = WGTopoHelper2d::GetSlightMoveDirectionByErrorCoincideSegmentIndex(
                &container, error_coincide_segment_index) * (distance_epsilon * 2);
            return Result(false, false, false, vt);
        }
        return Result(false, false, false, WGVector2d(0, 0));
    }
    return Result(true, false, false, WGVector2d(0, 0));
}

void WGPolygonPolygonBoolean::ExecuteSeparated(WGTopoHelper2d::SplitterContainer* container, const WGPolygon* polygon0, const WGPolygon* polygon1,
    bool is_box_separated, double distance_epsilon, int result_count, WGTopoHelper2d::SegmentSelector** selectors, 
    bool& is_singularity0, bool& is_singularity1, WGPolygon** result_polygons) {
    for (int result_index = 0; result_index < result_count; ++result_index) {
        result_polygons[result_index] = new WGPolygon();
        WGTopoHelper2d::AddPolygonSeparatedLoops(container, polygon0, polygon1, is_box_separated,
            selectors[result_index], distance_epsilon, is_singularity0, is_singularity1, result_polygons[result_index]);
    }
}