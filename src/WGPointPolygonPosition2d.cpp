#include "WGPointPolygonPosition2d.h"
#include "WGLine2d.h"
#include "WGCurveIntersecter2d.h"

/*
class WGPointPolygonPosition2dTool {
public:
    static bool QuickGetPointPolygonPosition(const WGVector2d& point, const WGPolygon* polygon, double distance_epsilon,
        const WGWire2d* ray_wire, WGTopoHelper2d::SplitterContainer* container, WGPointPolygonPosition2dType& result, int& confidence) {
        std::vector<WGTopoHelper2d::Splitter> splitters;
        splitters.reserve(16);
        std::vector<WGCurveCurveIntersection2d> intersections;
        const WGCurve2d* curve0 = ray_wire->GetEdge(0)->GetCurve();
        for (int i = 0; i < polygon->GetLoopCount(); ++i) {
            const WGWire2d* wire1 = polygon->GetLoop(i)->GetWire();
            for (int j = 0; j < wire1->GetEdgeCount(); ++j) {
                const WGEdge2d* edge1 = wire1->GetEdge(j);
                const WGCurve2d* curve1 = edge1->GetCurve();
                curve_curve_intersect(curve0, curve1, distance_epsilon, true, intersections);
                for (auto intersection_itr = intersections.begin(); intersection_itr != intersections.end(); ++intersection_itr) {
                    WGCurveCurveIntersection2d& intersection = *intersection_itr;
                    if (intersection.StartPoint.Ts[0] == 0) {
                        result = WGPointPolygonPosition2dType::On;
                        confidence = WGTopoHelper2d::HighestConfidence;
                        return true;
                    }
                    splitters.push_back(WGTopoHelper2d::NewSplitter(container, ray_wire, 0, wire1, j, &intersection));
                }
                intersections.clear();
            }
        }
        if (splitters.size() == 0) {
            result = WGPointPolygonPosition2dType::Outter;
            confidence = WGTopoHelper2d::HighestConfidence;
            return true;
        }
        WGTopoHelper2d::AddSplitters(container, splitters);
        WGTopoHelper2d::AddCoincideSideSplitter(container);
        WGTopoHelper2d::SetSplittersIsPoint(container, distance_epsilon);
        WGTopoHelper2d::BuildMergedSplitters(container);
        WGTopoHelper2d::BuildMergedSplitterPositions(container, true, true);
        int confidence0 = WGTopoHelper2d::HighestConfidence;
        int count0 = 0;
        int count1 = 0;
        for (int i = 0; i < container->MergedSplitters.size(); ++i) {
            const WGTopoHelper2d::MergedSplitterPosition& position = container->MergedSplitterPositions[1].at(i);
            if (position.PrevPosition == WGTopoHelper2d::SegmentPosition::Inner) {
                ++count0;
            }
            else if (position.PrevPosition == WGTopoHelper2d::SegmentPosition::Outter) {
                ++count1;
            }
            if (position.NextPosition == WGTopoHelper2d::SegmentPosition::Inner) {
                ++count0;
            }
            else if (position.NextPosition == WGTopoHelper2d::SegmentPosition::Outter) {
                ++count1;
            }
            if (position.PrevConfidence < confidence0) {
                confidence0 = position.PrevConfidence;
            }
            if (position.NextConfidence < confidence0) {
                confidence0 = position.NextConfidence;
            }
        }
        count0 = count0 & 1;
        count1 = count1 & 1;
        WGPointPolygonPosition2dType result0;
        if (count0 != count1) {
            result0 = WGPointPolygonPosition2dType::Unknown;
            confidence0 = 0;
        }
        else if (count0 == 0) {
            result0 = WGPointPolygonPosition2dType::Outter;
        }
        else {
            result0 = WGPointPolygonPosition2dType::Inner;
        }
        if (confidence0 >= WGTopoHelper2d::HighConfidence) {
            result = result0;
            confidence = confidence0;
            return true;
        }
        int confidence1 = confidence0;
        WGTopoHelper2d::InitializeIntersections(container);
        WGTopoHelper2d::Intersection intersection = WGTopoHelper2d::BuildMinIntersection(container, 0);
        WGTopoHelper2d::MergeIntersectionPart(container, 0, intersection);
        if (intersection.Parts[0].size() > 0 || intersection.Parts[1].size() > 0) {
            result = result0;
            confidence = confidence0;
            return false;
        }
        const WGTopoHelper2d::MergedSplitterPosition& position = intersection.Part[0].Position;
        WGPointPolygonPosition2dType result2;
        if (position.PrevPosition == WGTopoHelper2d::SegmentPosition::Inner) {
            result2 = WGPointPolygonPosition2dType::Inner;
        }
        else if (position.PrevPosition == WGTopoHelper2d::SegmentPosition::Outter) {
            result2 = WGPointPolygonPosition2dType::Outter;
        }
        else {
            result2 = WGPointPolygonPosition2dType::Unknown;
        }
        if (result2 == WGPointPolygonPosition2dType::Unknown) {
            result = result0;
            confidence = confidence0;
            return false;
        }
        int confidence2 = position.PrevConfidence;
        if (confidence2 >= WGTopoHelper2d::HighConfidence) {
            result = result2;
            confidence = confidence2;
            return true;
        }
        if (result0 == result2) {
            result = result0;
            confidence = confidence0 + confidence1 + confidence2;
            return confidence >= WGTopoHelper2d::HighConfidence;
        }
        result = WGPointPolygonPosition2dType::Unknown;
        confidence = 0;
        return false;
    }

    static int SampleX(double base_x, double epsilon, int sample_count_per_piece,
        const WGTopoHelper2d::SplitPoint& split_point0, const WGTopoHelper2d::SplitPoint& split_point1) {
        double sum = 0;
        WGTopoHelper2d::SampleIterator iterator(split_point0.Key, split_point1.Key, sample_count_per_piece);
        while (!iterator.Eof()) {
            WGVector2d point = iterator.GetCurrentPoint();
            double d = point.X - base_x;
            if (d > epsilon) {
                return 1;
            }
            if (d < -epsilon) {
                return -1;
            }
            sum += d;
        }
        return sum > 0 ? 1 : -1;
    }
};

WGPointPolygonPosition2dType get_point_polygon_position(const WGVector2d& point, const WGPolygon* polygon, double distance_epsilon) {
    double joint_epsilon = polygon->GetJointEpsilon();
    if (joint_epsilon > distance_epsilon) {
        distance_epsilon = joint_epsilon;
    }
    const WGBox2d& box = polygon->GetBox();
    if (!box.IsIntersected(point, distance_epsilon)) {
        return WGPointPolygonPosition2dType::Outter;
    }
    double y_length = (box.GetMax().Y - box.GetMin().Y) * 2;
    WGTopoHelper2d::SplitterContainer container0;
    WGPointPolygonPosition2dType result0;
    int confidence0;
    WGWire2d ray_wire0;
    ray_wire0.AddEdge(new WGEdge2d(new WGLine2d(point, WGVector2d(point.X, y_length))));
    if (WGPointPolygonPosition2dTool::QuickGetPointPolygonPosition(point, polygon, distance_epsilon, &ray_wire0, &container0, result0, confidence0)) {
        return result0;
    }
    WGTopoHelper2d::SplitterContainer container1;
    WGPointPolygonPosition2dType result1;
    int confidence1;
    WGWire2d ray_wire1;
    ray_wire1.AddEdge(new WGEdge2d(new WGLine2d(point, WGVector2d(point.X, -y_length))));
    if (WGPointPolygonPosition2dTool::QuickGetPointPolygonPosition(point, polygon, distance_epsilon, &ray_wire1, &container1, result1, confidence1)) {
        return result1;
    }
    if (result0 == result1 && result0 != WGPointPolygonPosition2dType::Unknown) {
        if (confidence0 + confidence1 > WGTopoHelper2d::HighConfidence) {
            return result0;
        }
    }
    std::vector<WGTopoHelper2d::SplitPoint> split_points;
    split_points.reserve(container0.SplitPoints[1].size() + container1.SplitPoints[1].size());
    split_points.insert(split_points.end(), container0.SplitPoints[1].begin(), container0.SplitPoints[1].end());
    split_points.insert(split_points.end(), container1.SplitPoints[1].begin(), container1.SplitPoints[1].end());
    std::sort(split_points.begin(), split_points.end());
    int count = 0;
    int s = 0;
    while (s < (int)split_points.size()) {
        const WGWire2d* wire = split_points.at(s).Key.Wire;
        int e = s + 1;
        while (e < (int)split_points.size()) {
            if (split_points.at(e).Key.Wire != wire) {
                break;
            }
            ++e;
        }
        if (e - s > 1) {
            for (int i = s; i < e; ++i) {
                int j = i + 1;
                if (j == e) {
                    j = s;
                }
                const WGTopoHelper2d::SplitPoint& split_point0 = split_points.at(i);
                const WGTopoHelper2d::SplitPoint& split_point1 = split_points.at(j);
                if ((split_point0.Point.Y - point.Y) * (split_point1.Point.Y - point.Y) < 0) {
                    int dx = WGPointPolygonPosition2dTool::SampleX(point.X, distance_epsilon * 2, 8, split_point0, split_point1);
                    if (dx > 0) {
                        ++count;
                    }
                }
            }
        }
        s = e;
    }
    count = count & 1;
    return count == 0 ? WGPointPolygonPosition2dType::Outter : WGPointPolygonPosition2dType::Inner;
}
*/

WGPointPolygonPosition2d::Result WGPointPolygonPosition2d::Execute(const WGVector2d& point, const WGPolygon* polygon, double distance_epsilon) {
    WGIntersectHelper2d::IntersectCache intersect_cache;
    return Execute(point, polygon, distance_epsilon, &intersect_cache);
}

WGPointPolygonPosition2d::Result WGPointPolygonPosition2d::Execute(const WGVector2d& point, const WGPolygon* polygon, double distance_epsilon,
    WGIntersectHelper2d::IntersectCache* intersect_cache) {
    double joint_epsilon = polygon->GetJointEpsilon();
    if (joint_epsilon > distance_epsilon) {
        distance_epsilon = joint_epsilon;
    }
    const WGBox2d& box = polygon->GetBox();
    if (!box.IsIntersected(point, distance_epsilon)) {
        return WGPointPolygonPosition2d::Result::Outter;
    }
    double y_length = (box.GetMax().Y - box.GetMin().Y) * 2;
    WGTopoHelper2d::SplitterContainer container0;
    WGPointPolygonPosition2d::Result result0;
    int confidence0;
    WGWire2d ray_wire0;
    ray_wire0.AddEdge(new WGEdge2d(new WGLine2d(point, WGVector2d(point.X, y_length))));
    if (QuickGetPointPolygonPosition(point, polygon, distance_epsilon, &ray_wire0, &container0, intersect_cache, result0, confidence0)) {
        return result0;
    }
    WGTopoHelper2d::SplitterContainer container1;
    WGPointPolygonPosition2d::Result result1;
    int confidence1;
    WGWire2d ray_wire1;
    ray_wire1.AddEdge(new WGEdge2d(new WGLine2d(point, WGVector2d(point.X, -y_length))));
    if (QuickGetPointPolygonPosition(point, polygon, distance_epsilon, &ray_wire1, &container1, intersect_cache, result1, confidence1)) {
        return result1;
    }
    if (result0 == result1 && result0 != WGPointPolygonPosition2d::Result::Unknown) {
        if (confidence0 + confidence1 > WGTopoHelper2d::HighConfidence) {
            return result0;
        }
    }
    std::vector<WGTopoHelper2d::SplitPoint> split_points;
    split_points.reserve(container0.SplitPoints[1].size() + container1.SplitPoints[1].size());
    split_points.insert(split_points.end(), container0.SplitPoints[1].begin(), container0.SplitPoints[1].end());
    split_points.insert(split_points.end(), container1.SplitPoints[1].begin(), container1.SplitPoints[1].end());
    std::sort(split_points.begin(), split_points.end());
    int count = 0;
    int s = 0;
    while (s < (int)split_points.size()) {
        const WGWire2d* wire = split_points.at(s).Key.Wire;
        int e = s + 1;
        while (e < (int)split_points.size()) {
            if (split_points.at(e).Key.Wire != wire) {
                break;
            }
            ++e;
        }
        if (e - s > 1) {
            for (int i = s; i < e; ++i) {
                int j = i + 1;
                if (j == e) {
                    j = s;
                }
                const WGTopoHelper2d::SplitPoint& split_point0 = split_points.at(i);
                const WGTopoHelper2d::SplitPoint& split_point1 = split_points.at(j);
                if ((split_point0.Point.Y - point.Y) * (split_point1.Point.Y - point.Y) < 0) {
                    int dx = SampleX(point.X, distance_epsilon * 2, 8, &split_point0, &split_point1);
                    if (dx > 0) {
                        ++count;
                    }
                }
            }
        }
        s = e;
    }
    count = count & 1;
    return count == 0 ? WGPointPolygonPosition2d::Result::Outter : WGPointPolygonPosition2d::Result::Inner;
}

bool WGPointPolygonPosition2d::QuickGetPointPolygonPosition(const WGVector2d& point, const WGPolygon* polygon,
    double distance_epsilon, const WGWire2d* ray_wire, WGTopoHelper2d::SplitterContainer* container,
    WGIntersectHelper2d::IntersectCache* intersect_cache, WGPointPolygonPosition2d::Result& result, int& confidence) {
    std::vector<WGTopoHelper2d::Splitter> splitters;
    splitters.reserve(16);
    std::vector<WGCurveCurveIntersection2d> intersections;
    const WGCurve2d* curve0 = ray_wire->GetEdge(0)->GetCurve();
    for (int i = 0; i < polygon->GetLoopCount(); ++i) {
        const WGWire2d* wire1 = polygon->GetLoop(i)->GetWire();
        for (int j = 0; j < wire1->GetEdgeCount(); ++j) {
            const WGEdge2d* edge1 = wire1->GetEdge(j);
            const WGCurve2d* curve1 = edge1->GetCurve();
            bool is_singularity0;
            bool is_singularity1;
            WGCurveIntersecter2d::Intersect(curve0, curve1, distance_epsilon, intersect_cache, is_singularity0, is_singularity1, intersections);
            assert(!is_singularity0);
            if (is_singularity1) {
                result = WGPointPolygonPosition2d::Result::Singularity;
                confidence = WGTopoHelper2d::HighestConfidence;
                return true;
            }
            for (auto intersection_itr = intersections.begin(); intersection_itr != intersections.end(); ++intersection_itr) {
                WGCurveCurveIntersection2d& intersection = *intersection_itr;
                if (intersection.Ts0[0] == 0) {
                    result = WGPointPolygonPosition2d::Result::On;
                    confidence = WGTopoHelper2d::HighestConfidence;
                    return true;
                }
                splitters.push_back(WGTopoHelper2d::NewSplitter(container, ray_wire, 0, wire1, j, &intersection));
            }
            intersections.clear();
        }
    }
    if (splitters.size() == 0) {
        result = WGPointPolygonPosition2d::Result::Outter;
        confidence = WGTopoHelper2d::HighestConfidence;
        return true;
    }
    WGTopoHelper2d::AddSplitters(container, splitters);
    WGTopoHelper2d::AddCoincideSideSplitter(container);
    WGTopoHelper2d::SetSplittersIsPoint(container, distance_epsilon);
    WGTopoHelper2d::BuildMergedSplitters(container);
    WGTopoHelper2d::BuildMergedSplitterPositions(container, true, true);
    int confidence0 = WGTopoHelper2d::HighestConfidence;
    int count0 = 0;
    int count1 = 0;
    for (int i = 0; i < container->MergedSplitters.size(); ++i) {
        const WGTopoHelper2d::MergedSplitterPosition& position = container->MergedSplitterPositions[1].at(i);
        if (position.PrevPosition == WGTopoHelper2d::SegmentPosition::Inner) {
            ++count0;
        }
        else if (position.PrevPosition == WGTopoHelper2d::SegmentPosition::Outter) {
            ++count1;
        }
        if (position.NextPosition == WGTopoHelper2d::SegmentPosition::Inner) {
            ++count0;
        }
        else if (position.NextPosition == WGTopoHelper2d::SegmentPosition::Outter) {
            ++count1;
        }
        if (position.PrevConfidence < confidence0) {
            confidence0 = position.PrevConfidence;
        }
        if (position.NextConfidence < confidence0) {
            confidence0 = position.NextConfidence;
        }
    }
    count0 = count0 & 1;
    count1 = count1 & 1;
    WGPointPolygonPosition2d::Result result0;
    if (count0 != count1) {
        result0 = WGPointPolygonPosition2d::Result::Unknown;
        confidence0 = 0;
    }
    else if (count0 == 0) {
        result0 = WGPointPolygonPosition2d::Result::Outter;
    }
    else {
        result0 = WGPointPolygonPosition2d::Result::Inner;
    }
    if (confidence0 >= WGTopoHelper2d::HighConfidence) {
        result = result0;
        confidence = confidence0;
        return true;
    }
    int confidence1 = confidence0;
    WGTopoHelper2d::InitializeIntersections(container);
    WGTopoHelper2d::Intersection intersection = WGTopoHelper2d::BuildMinIntersection(container, 0);
    WGTopoHelper2d::MergeIntersectionPart(container, 0, intersection);
    if (intersection.Parts[0].size() > 0 || intersection.Parts[1].size() > 0) {
        result = result0;
        confidence = confidence0;
        return false;
    }
    const WGTopoHelper2d::MergedSplitterPosition& position = intersection.Part[0].Position;
    WGPointPolygonPosition2d::Result result2;
    if (position.PrevPosition == WGTopoHelper2d::SegmentPosition::Inner) {
        result2 = WGPointPolygonPosition2d::Result::Inner;
    }
    else if (position.PrevPosition == WGTopoHelper2d::SegmentPosition::Outter) {
        result2 = WGPointPolygonPosition2d::Result::Outter;
    }
    else {
        result2 = WGPointPolygonPosition2d::Result::Unknown;
    }
    if (result2 == WGPointPolygonPosition2d::Result::Unknown) {
        result = result0;
        confidence = confidence0;
        return false;
    }
    int confidence2 = position.PrevConfidence;
    if (confidence2 >= WGTopoHelper2d::HighConfidence) {
        result = result2;
        confidence = confidence2;
        return true;
    }
    if (result0 == result2) {
        result = result0;
        confidence = confidence0 + confidence1 + confidence2;
        return confidence >= WGTopoHelper2d::HighConfidence;
    }
    result = WGPointPolygonPosition2d::Result::Unknown;
    confidence = 0;
    return false;
}

int WGPointPolygonPosition2d::SampleX(double base_x, double epsilon, int sample_count_per_piece,
    const WGTopoHelper2d::SplitPoint* split_point0, const WGTopoHelper2d::SplitPoint* split_point1) {
    double sum = 0;
    WGTopoHelper2d::SampleIterator iterator(split_point0->Key, split_point1->Key, sample_count_per_piece);
    while (!iterator.Eof()) {
        WGVector2d point = iterator.GetCurrentPoint();
        double d = point.X - base_x;
        if (d > epsilon) {
            return 1;
        }
        if (d < -epsilon) {
            return -1;
        }
        sum += d;
    }
    return sum > 0 ? 1 : -1;
}