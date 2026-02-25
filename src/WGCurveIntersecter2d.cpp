#include "WGCurveIntersecter2d.h"
#include "WGLine2d.h"
#include "WGArc2d.h"
#include "WGNurbsCurve2d.h"
#include "WGNurbs.h"
#include <assert.h>
#include <algorithm>

void WGCurveIntersecter2d::Intersect(const WGVector2d& point, const WGCurve2d* curve, int piece_index, int piece_count,
    double distance_epsilon, WGIntersectHelper2d::IntersectCache* intersect_cache, bool& is_singularity,
    std::vector<WGPointCurveIntersection2d>& intersections) {
    if (curve->GetType() == WGLine2d::Type::Instance()) {
        assert(piece_index == 0 && piece_count == 1);
        WGLine2d* line = (WGLine2d*)curve;
        WGIntersectHelper2d::PointCurveIntersection base_intersections[1];
        int n = WGIntersectHelper2d::PointLineIntersect(point, line->GetStartPoint(), line->GetEndPoint(),
            distance_epsilon, is_singularity, base_intersections);
        if (is_singularity) {
            return;
        }
        AddIntersections(intersections, base_intersections, n, piece_index);
    }
    else if (curve->GetType() == WGArc2d::Type::Instance()) {
        assert(piece_index == 0 && piece_count == 1);
        WGArc2d* arc = (WGArc2d*)curve;
        WGIntersectHelper2d::PointCurveIntersection base_intersections[2];
        int n = WGIntersectHelper2d::PointArcIntersect(point, arc->GetCenter(), arc->GetRadius(), arc->GetStartAngle(), 
            arc->GetDeltaAngle(), distance_epsilon, is_singularity, base_intersections);
        if (is_singularity) {
            return;
        }
        AddIntersections(intersections, base_intersections, n, piece_index);
    }
    else if (curve->GetType() == WGNurbsCurve2d::Type::Instance()) {
        WGNurbsCurve2d* nurbs = (WGNurbsCurve2d*)curve;
        if (nurbs->GetWeights()) {
            std::vector<WGIntersectHelper2d::PointCurveIntersection> base_intersections;
            WGNurbsCurve2d::BezierIterator bezier_iterator(nurbs, piece_index, piece_count);
            while (!bezier_iterator.Eof()) {
                WGIntersectHelper2d::PointRationalBezierCurveIntersect(point,
                    nurbs->GetDegree(), bezier_iterator.GetCurrentColtrolPoints(), bezier_iterator.GetCurrentWeights(),
                    distance_epsilon, intersect_cache, is_singularity, base_intersections);
                if (is_singularity) {
                    return;
                }
                if (base_intersections.size() > 0) {
                    int current_piece_index = bezier_iterator.GetCurrentPieceIndex();
                    WSInterval current_domain = bezier_iterator.GetCurrentDomain();
                    AddIntersections(intersections, base_intersections, current_domain, current_piece_index);
                    base_intersections.clear();
                }
                bezier_iterator.Next();
            }
        }
        else {
            std::vector<WGIntersectHelper2d::PointCurveIntersection> base_intersections;
            WGNurbsCurve2d::BezierIterator bezier_iterator(nurbs, piece_index, piece_count);
            while (!bezier_iterator.Eof()) {
                WGIntersectHelper2d::PointBezierCurveIntersect(point, nurbs->GetDegree(), bezier_iterator.GetCurrentColtrolPoints(),
                    distance_epsilon, intersect_cache, is_singularity, base_intersections);
                if (is_singularity) {
                    return;
                }
                if (base_intersections.size() > 0) {
                    int current_piece_index = bezier_iterator.GetCurrentPieceIndex();
                    WSInterval current_domain = bezier_iterator.GetCurrentDomain();
                    AddIntersections(intersections, base_intersections, current_domain, current_piece_index);
                    base_intersections.clear();
                }
                bezier_iterator.Next();
            }
        }
    }
}

void WGCurveIntersecter2d::Intersect(const WGCurve2d* curve0, const WGCurve2d* curve1, double distance_epsilon,
    WGIntersectHelper2d::IntersectCache* intersect_cache, bool& is_singularity0, bool& is_singularity1, 
    std::vector<WGCurveCurveIntersection2d>& intersections) {
    Intersect(curve0, 0, curve0->GetPieceCount(), curve1, 0, curve1->GetPieceCount(), distance_epsilon, 
        intersect_cache, is_singularity0, is_singularity1, intersections);
}

void WGCurveIntersecter2d::Intersect(const WGCurve2d* curve0, int piece_index0, int piece_count0,
    const WGCurve2d* curve1, int piece_index1, int piece_count1, double distance_epsilon,
    WGIntersectHelper2d::IntersectCache* intersect_cache, bool& is_singularity0, bool& is_singularity1, 
    std::vector<WGCurveCurveIntersection2d>& intersections) {
    if (curve0->GetType() == WGLine2d::Type::Instance()) {
        assert(piece_index0 == 0 && piece_count0 == 1);
        WGLine2d* line0 = (WGLine2d*)curve0;
        if (curve1->GetType() == WGLine2d::Type::Instance()) {
            assert(piece_index1 == 0 && piece_count1 == 1);
            WGLine2d* line1 = (WGLine2d*)curve1;
            WGIntersectHelper2d::CurveCurveIntersection base_intersections[2];
            int n = WGIntersectHelper2d::LineLineIntersect(line0->GetStartPoint(), line0->GetEndPoint(),
                line1->GetStartPoint(), line1->GetEndPoint(), distance_epsilon, is_singularity0, is_singularity1, base_intersections);
            if (is_singularity0 || is_singularity1) {
                return;
            }
            AddIntersections(intersections, base_intersections, n, piece_index0, piece_index1, false);
        }
        else if (curve1->GetType() == WGArc2d::Type::Instance()) {
            assert(piece_index1 == 0 && piece_count1 == 1);
            WGArc2d* arc1 = (WGArc2d*)curve1;
            WGIntersectHelper2d::CurveCurveIntersection base_intersections[10];
            int n = WGIntersectHelper2d::LineArcIntersect(line0->GetStartPoint(), line0->GetEndPoint(),
                arc1->GetCenter(), arc1->GetRadius(), arc1->GetStartAngle(), arc1->GetDeltaAngle(), 
                distance_epsilon, is_singularity0, is_singularity1, base_intersections);
            if (is_singularity0 || is_singularity1) {
                return;
            }
            AddIntersections(intersections, base_intersections, n, piece_index0, piece_index1, false);
        }
        else if (curve1->GetType() == WGNurbsCurve2d::Type::Instance()) {
            WGNurbsCurve2d* nurbs1 = (WGNurbsCurve2d*)curve1;
            if (nurbs1->GetWeights()) {
                std::vector<WGIntersectHelper2d::CurveCurveIntersection> base_intersections;
                WGNurbsCurve2d::BezierIterator bezier_iterator1(nurbs1, piece_index1, piece_count1);
                while (!bezier_iterator1.Eof()) {
                    WGIntersectHelper2d::LineRationalBezierCurveIntersect(line0->GetStartPoint(), line0->GetEndPoint(),
                        nurbs1->GetDegree(), bezier_iterator1.GetCurrentColtrolPoints(), bezier_iterator1.GetCurrentWeights(), 
                        distance_epsilon, intersect_cache, is_singularity0, is_singularity1, base_intersections);
                    if (is_singularity0 || is_singularity1) {
                        return;
                    }
                    if (base_intersections.size() > 0) {
                        int current_piece_index1 = bezier_iterator1.GetCurrentPieceIndex();
                        WSInterval current_domain1 = bezier_iterator1.GetCurrentDomain();
                        AddIntersections(intersections, base_intersections, nullptr, &current_domain1, piece_index0, current_piece_index1, false);
                        base_intersections.clear();
                    }
                    bezier_iterator1.Next();
                }
            }
            else {
                std::vector<WGIntersectHelper2d::CurveCurveIntersection> base_intersections;
                WGNurbsCurve2d::BezierIterator bezier_iterator1(nurbs1, piece_index1, piece_count1);
                while (!bezier_iterator1.Eof()) {
                    WGIntersectHelper2d::LineBezierCurveIntersect(line0->GetStartPoint(), line0->GetEndPoint(),
                        nurbs1->GetDegree(), bezier_iterator1.GetCurrentColtrolPoints(), distance_epsilon,
                        intersect_cache, is_singularity0, is_singularity1, base_intersections);
                    if (is_singularity0 || is_singularity1) {
                        return;
                    }
                    if (base_intersections.size() > 0) {
                        int current_piece_index1 = bezier_iterator1.GetCurrentPieceIndex();
                        WSInterval current_domain1 = bezier_iterator1.GetCurrentDomain();
                        AddIntersections(intersections, base_intersections, nullptr, &current_domain1, piece_index0, current_piece_index1, false);
                        base_intersections.clear();
                    }
                    bezier_iterator1.Next();
                }
            }
        }
    }
    else if (curve0->GetType() == WGArc2d::Type::Instance()) {
        assert(piece_index0 == 0 && piece_count0 == 1);
        WGArc2d* arc0 = (WGArc2d*)curve0;
        if (curve1->GetType() == WGLine2d::Type::Instance()) {
            assert(piece_index1 == 0 && piece_count1 == 1);
            WGLine2d* line1 = (WGLine2d*)curve1;
            WGIntersectHelper2d::CurveCurveIntersection base_intersections[10];
            int n = WGIntersectHelper2d::LineArcIntersect(line1->GetStartPoint(), line1->GetEndPoint(),
                arc0->GetCenter(), arc0->GetRadius(), arc0->GetStartAngle(), arc0->GetDeltaAngle(),
                distance_epsilon, is_singularity1, is_singularity0, base_intersections);
            if (is_singularity0 || is_singularity1) {
                return;
            }
            AddIntersections(intersections, base_intersections, n, piece_index1, piece_index0, true);
        }
        else if (curve1->GetType() == WGArc2d::Type::Instance()) {
            assert(piece_index1 == 0 && piece_count1 == 1);
            WGArc2d* arc1 = (WGArc2d*)curve1;
            WGIntersectHelper2d::CurveCurveIntersection base_intersections[10];
            int n = WGIntersectHelper2d::ArcArcIntersect(arc0->GetCenter(), arc0->GetRadius(), arc0->GetStartAngle(), arc0->GetDeltaAngle(),
                arc1->GetCenter(), arc1->GetRadius(), arc1->GetStartAngle(), arc1->GetDeltaAngle(),
                distance_epsilon, is_singularity0, is_singularity1, base_intersections);
            if (is_singularity0 || is_singularity1) {
                return;
            }
            AddIntersections(intersections, base_intersections, n, piece_index0, piece_index1, false);
        }
        else if (curve1->GetType() == WGNurbsCurve2d::Type::Instance()) {
            WGNurbsCurve2d* nurbs1 = (WGNurbsCurve2d*)curve1;
            if (nurbs1->GetWeights()) {
                std::vector<WGIntersectHelper2d::CurveCurveIntersection> base_intersections;
                WGNurbsCurve2d::BezierIterator bezier_iterator1(nurbs1, piece_index1, piece_count1);
                while (!bezier_iterator1.Eof()) {
                    WGIntersectHelper2d::ArcRationalBezierCurveIntersect(arc0->GetCenter(), arc0->GetRadius(), arc0->GetStartAngle(), arc0->GetDeltaAngle(),
                        nurbs1->GetDegree(), bezier_iterator1.GetCurrentColtrolPoints(), bezier_iterator1.GetCurrentWeights(),
                        distance_epsilon, intersect_cache, is_singularity0, is_singularity1, base_intersections);
                    if (is_singularity0 || is_singularity1) {
                        return;
                    }
                    if (base_intersections.size() > 0) {
                        int current_piece_index1 = bezier_iterator1.GetCurrentPieceIndex();
                        WSInterval current_domain1 = bezier_iterator1.GetCurrentDomain();
                        AddIntersections(intersections, base_intersections, nullptr, &current_domain1, piece_index0, current_piece_index1, false);
                        base_intersections.clear();
                    }
                    bezier_iterator1.Next();
                }
            }
            else {
                std::vector<WGIntersectHelper2d::CurveCurveIntersection> base_intersections;
                WGNurbsCurve2d::BezierIterator bezier_iterator1(nurbs1, piece_index1, piece_count1);
                while (!bezier_iterator1.Eof()) {
                    WGIntersectHelper2d::ArcBezierCurveIntersect(arc0->GetCenter(), arc0->GetRadius(), arc0->GetStartAngle(), arc0->GetDeltaAngle(),
                        nurbs1->GetDegree(), bezier_iterator1.GetCurrentColtrolPoints(), distance_epsilon,
                        intersect_cache, is_singularity0, is_singularity1, base_intersections);
                    if (is_singularity0 || is_singularity1) {
                        return;
                    }
                    if (base_intersections.size() > 0) {
                        int current_piece_index1 = bezier_iterator1.GetCurrentPieceIndex();
                        WSInterval current_domain1 = bezier_iterator1.GetCurrentDomain();
                        AddIntersections(intersections, base_intersections, nullptr, &current_domain1, piece_index0, current_piece_index1, false);
                        base_intersections.clear();
                    }
                    bezier_iterator1.Next();
                }
            }
        }
    }
    else if (curve0->GetType() == WGNurbsCurve2d::Type::Instance()) {
        WGNurbsCurve2d* nurbs0 = (WGNurbsCurve2d*)curve0;
        if (nurbs0->GetWeights()) {
            if (curve1->GetType() == WGLine2d::Type::Instance()) {
                assert(piece_index1 == 0 && piece_count1 == 1);
                WGLine2d* line1 = (WGLine2d*)curve1;
                std::vector<WGIntersectHelper2d::CurveCurveIntersection> base_intersections;
                WGNurbsCurve2d::BezierIterator bezier_iterator0(nurbs0, piece_index0, piece_count0);
                while (!bezier_iterator0.Eof()) {
                    WGIntersectHelper2d::LineRationalBezierCurveIntersect(line1->GetStartPoint(), line1->GetEndPoint(),
                        nurbs0->GetDegree(), bezier_iterator0.GetCurrentColtrolPoints(), bezier_iterator0.GetCurrentWeights(),
                        distance_epsilon, intersect_cache, is_singularity1, is_singularity0, base_intersections);
                    if (is_singularity0 || is_singularity1) {
                        return;
                    }
                    if (base_intersections.size() > 0) {
                        int current_piece_index0 = bezier_iterator0.GetCurrentPieceIndex();
                        WSInterval current_domain0 = bezier_iterator0.GetCurrentDomain();
                        AddIntersections(intersections, base_intersections, nullptr, &current_domain0, piece_index1, current_piece_index0, true);
                        base_intersections.clear();
                    }
                    bezier_iterator0.Next();
                }
            }
            else if (curve1->GetType() == WGArc2d::Type::Instance()) {
                assert(piece_index1 == 0 && piece_count1 == 1);
                WGArc2d* arc1 = (WGArc2d*)curve1;
                std::vector<WGIntersectHelper2d::CurveCurveIntersection> base_intersections;
                WGNurbsCurve2d::BezierIterator bezier_iterator0(nurbs0, piece_index0, piece_count0);
                while (!bezier_iterator0.Eof()) {
                    WGIntersectHelper2d::ArcRationalBezierCurveIntersect(arc1->GetCenter(), arc1->GetRadius(), arc1->GetStartAngle(), arc1->GetDeltaAngle(),
                        nurbs0->GetDegree(), bezier_iterator0.GetCurrentColtrolPoints(), bezier_iterator0.GetCurrentWeights(),
                        distance_epsilon, intersect_cache, is_singularity1, is_singularity0, base_intersections);
                    if (is_singularity0 || is_singularity1) {
                        return;
                    }
                    if (base_intersections.size() > 0) {
                        int current_piece_index0 = bezier_iterator0.GetCurrentPieceIndex();
                        WSInterval current_domain0 = bezier_iterator0.GetCurrentDomain();
                        AddIntersections(intersections, base_intersections, nullptr, &current_domain0, piece_index1, current_piece_index0, true);
                        base_intersections.clear();
                    }
                    bezier_iterator0.Next();
                }
            }
            else if (curve1->GetType() == WGNurbsCurve2d::Type::Instance()) {
                WGNurbsCurve2d* nurbs1 = (WGNurbsCurve2d*)curve1;
                if (nurbs1->GetWeights()) {
                    std::vector<WGIntersectHelper2d::CurveCurveIntersection> base_intersections;
                    WGNurbsCurve2d::BezierIterator bezier_iterator0(nurbs0, piece_index0, piece_count0);
                    while (!bezier_iterator0.Eof()) {
                        WGNurbsCurve2d::BezierIterator bezier_iterator1(nurbs1, piece_index1, piece_count1);
                        while (!bezier_iterator1.Eof()) {
                            WGIntersectHelper2d::RationalBezierCurveRationalBezierCurveIntersect(
                                nurbs0->GetDegree(), bezier_iterator0.GetCurrentColtrolPoints(), bezier_iterator0.GetCurrentWeights(),
                                nurbs1->GetDegree(), bezier_iterator1.GetCurrentColtrolPoints(), bezier_iterator1.GetCurrentWeights(),
                                distance_epsilon, intersect_cache, is_singularity0, is_singularity1, base_intersections);
                            if (is_singularity0 || is_singularity1) {
                                return;
                            }
                            if (base_intersections.size() > 0) {
                                int current_piece_index0 = bezier_iterator0.GetCurrentPieceIndex();
                                WSInterval current_domain0 = bezier_iterator0.GetCurrentDomain();
                                int current_piece_index1 = bezier_iterator1.GetCurrentPieceIndex();
                                WSInterval current_domain1 = bezier_iterator1.GetCurrentDomain();
                                AddIntersections(intersections, base_intersections, &current_domain0, &current_domain1, 
                                    current_piece_index0, current_piece_index1, false);
                                base_intersections.clear();
                            }
                            bezier_iterator1.Next();
                        }
                        bezier_iterator0.Next();
                    }
                }
                else {
                    std::vector<WGIntersectHelper2d::CurveCurveIntersection> base_intersections;
                    WGNurbsCurve2d::BezierIterator bezier_iterator0(nurbs0, piece_index0, piece_count0);
                    while (!bezier_iterator0.Eof()) {
                        WGNurbsCurve2d::BezierIterator bezier_iterator1(nurbs1, piece_index1, piece_count1);
                        while (!bezier_iterator1.Eof()) {
                            WGIntersectHelper2d::BezierCurveRationalBezierCurveIntersect(
                                nurbs1->GetDegree(), bezier_iterator1.GetCurrentColtrolPoints(),
                                nurbs0->GetDegree(), bezier_iterator0.GetCurrentColtrolPoints(), bezier_iterator0.GetCurrentWeights(),
                                distance_epsilon, intersect_cache, is_singularity1, is_singularity0, base_intersections);
                            if (is_singularity0 || is_singularity1) {
                                return;
                            }
                            if (base_intersections.size() > 0) {
                                int current_piece_index0 = bezier_iterator0.GetCurrentPieceIndex();
                                WSInterval current_domain0 = bezier_iterator0.GetCurrentDomain();
                                int current_piece_index1 = bezier_iterator1.GetCurrentPieceIndex();
                                WSInterval current_domain1 = bezier_iterator1.GetCurrentDomain();
                                AddIntersections(intersections, base_intersections, &current_domain1, &current_domain0,
                                    current_piece_index1, current_piece_index0, true);
                                base_intersections.clear();
                            }
                            bezier_iterator1.Next();
                        }
                        bezier_iterator0.Next();
                    }
                }
            }
        }
        else {
            if (curve1->GetType() == WGLine2d::Type::Instance()) {
                assert(piece_index1 == 0 && piece_count1 == 1);
                WGLine2d* line1 = (WGLine2d*)curve1;
                std::vector<WGIntersectHelper2d::CurveCurveIntersection> base_intersections;
                WGNurbsCurve2d::BezierIterator bezier_iterator0(nurbs0, piece_index0, piece_count0);
                while (!bezier_iterator0.Eof()) {
                    WGIntersectHelper2d::LineBezierCurveIntersect(line1->GetStartPoint(), line1->GetEndPoint(),
                        nurbs0->GetDegree(), bezier_iterator0.GetCurrentColtrolPoints(),
                        distance_epsilon, intersect_cache, is_singularity1, is_singularity0, base_intersections);
                    if (is_singularity0 || is_singularity1) {
                        return;
                    }
                    if (base_intersections.size() > 0) {
                        int current_piece_index0 = bezier_iterator0.GetCurrentPieceIndex();
                        WSInterval current_domain0 = bezier_iterator0.GetCurrentDomain();
                        AddIntersections(intersections, base_intersections, nullptr, &current_domain0, piece_index1, current_piece_index0, true);
                        base_intersections.clear();
                    }
                    bezier_iterator0.Next();
                }
            }
            else if (curve1->GetType() == WGArc2d::Type::Instance()) {
                assert(piece_index1 == 0 && piece_count1 == 1);
                WGArc2d* arc1 = (WGArc2d*)curve1;
                std::vector<WGIntersectHelper2d::CurveCurveIntersection> base_intersections;
                WGNurbsCurve2d::BezierIterator bezier_iterator0(nurbs0, piece_index0, piece_count0);
                while (!bezier_iterator0.Eof()) {
                    WGIntersectHelper2d::ArcBezierCurveIntersect(arc1->GetCenter(), arc1->GetRadius(), arc1->GetStartAngle(), arc1->GetDeltaAngle(),
                        nurbs0->GetDegree(), bezier_iterator0.GetCurrentColtrolPoints(), 
                        distance_epsilon, intersect_cache, is_singularity1, is_singularity0, base_intersections);
                    if (is_singularity0 || is_singularity1) {
                        return;
                    }
                    if (base_intersections.size() > 0) {
                        int current_piece_index0 = bezier_iterator0.GetCurrentPieceIndex();
                        WSInterval current_domain0 = bezier_iterator0.GetCurrentDomain();
                        AddIntersections(intersections, base_intersections, nullptr, &current_domain0, piece_index1, current_piece_index0, true);
                        base_intersections.clear();
                    }
                    bezier_iterator0.Next();
                }
            }
            else if (curve1->GetType() == WGNurbsCurve2d::Type::Instance()) {
                WGNurbsCurve2d* nurbs1 = (WGNurbsCurve2d*)curve1;
                if (nurbs1->GetWeights()) {
                    std::vector<WGIntersectHelper2d::CurveCurveIntersection> base_intersections;
                    WGNurbsCurve2d::BezierIterator bezier_iterator0(nurbs0, piece_index0, piece_count0);
                    while (!bezier_iterator0.Eof()) {
                        WGNurbsCurve2d::BezierIterator bezier_iterator1(nurbs1, piece_index1, piece_count1);
                        while (!bezier_iterator1.Eof()) {
                            WGIntersectHelper2d::BezierCurveRationalBezierCurveIntersect(
                                nurbs0->GetDegree(), bezier_iterator0.GetCurrentColtrolPoints(),
                                nurbs1->GetDegree(), bezier_iterator1.GetCurrentColtrolPoints(), bezier_iterator1.GetCurrentWeights(),
                                distance_epsilon, intersect_cache, is_singularity0, is_singularity1, base_intersections);
                            if (is_singularity0 || is_singularity1) {
                                return;
                            }
                            if (base_intersections.size() > 0) {
                                int current_piece_index0 = bezier_iterator0.GetCurrentPieceIndex();
                                WSInterval current_domain0 = bezier_iterator0.GetCurrentDomain();
                                int current_piece_index1 = bezier_iterator1.GetCurrentPieceIndex();
                                WSInterval current_domain1 = bezier_iterator1.GetCurrentDomain();
                                AddIntersections(intersections, base_intersections, &current_domain0, &current_domain1,
                                    current_piece_index0, current_piece_index1, false);
                                base_intersections.clear();
                            }
                            bezier_iterator1.Next();
                        }
                        bezier_iterator0.Next();
                    }
                }
                else {
                    std::vector<WGIntersectHelper2d::CurveCurveIntersection> base_intersections;
                    WGNurbsCurve2d::BezierIterator bezier_iterator0(nurbs0, piece_index0, piece_count0);
                    while (!bezier_iterator0.Eof()) {
                        WGNurbsCurve2d::BezierIterator bezier_iterator1(nurbs1, piece_index1, piece_count1);
                        while (!bezier_iterator1.Eof()) {
                            WGIntersectHelper2d::BezierCurveBezierCurveIntersect(
                                nurbs0->GetDegree(), bezier_iterator0.GetCurrentColtrolPoints(),
                                nurbs1->GetDegree(), bezier_iterator1.GetCurrentColtrolPoints(),
                                distance_epsilon, intersect_cache, is_singularity0, is_singularity1, base_intersections);
                            if (is_singularity0 || is_singularity1) {
                                return;
                            }
                            if (base_intersections.size() > 0) {
                                int current_piece_index0 = bezier_iterator0.GetCurrentPieceIndex();
                                WSInterval current_domain0 = bezier_iterator0.GetCurrentDomain();
                                int current_piece_index1 = bezier_iterator1.GetCurrentPieceIndex();
                                WSInterval current_domain1 = bezier_iterator1.GetCurrentDomain();
                                AddIntersections(intersections, base_intersections, &current_domain0, &current_domain1,
                                    current_piece_index0, current_piece_index1, false);
                                base_intersections.clear();
                            }
                            bezier_iterator1.Next();
                        }
                        bezier_iterator0.Next();
                    }
                }
            }
        }
    }
}

void WGCurveIntersecter2d::RemoveSingularity(const WGCurve2d* curve, double distance_epsilon, std::vector<WGCurve2d*>& curves) {
    if (curve->GetType() == WGLine2d::Type::Instance()) {
        WGLine2d* line = (WGLine2d*)curve;
        if ((line->GetEndPoint() - line->GetStartPoint()).SqrLength() > distance_epsilon * distance_epsilon) {
            curves.push_back(line);
        }
    }
    else if (curve->GetType() == WGArc2d::Type::Instance()) {
        WGArc2d* arc = (WGArc2d*)curve;
        if (abs(arc->GetRadius() * arc->GetDeltaAngle()) > distance_epsilon) {
            curves.push_back(arc);
        }
    }
    else if (curve->GetType() == WGNurbsCurve2d::Type::Instance()) {
        WGNurbsCurve2d* nurbs = (WGNurbsCurve2d*)curve;
        std::vector<WSInterval> domains;
        if (nurbs->GetWeights()) {
            std::vector<WSInterval> base_domains;
            WGNurbsCurve2d::BezierIterator bezier_iterator(nurbs, 0, nurbs->GetPieceCount());
            while (!bezier_iterator.Eof()) {
                WGIntersectHelper2d::RemoveRationalBezierCurveSingularity(
                    nurbs->GetDegree(), bezier_iterator.GetCurrentColtrolPoints(), bezier_iterator.GetCurrentWeights(),
                    distance_epsilon, base_domains);
                MergeDomains(base_domains);
                WSInterval current_domain = bezier_iterator.GetCurrentDomain();
                for (auto base_domain_itr = base_domains.begin(); base_domain_itr != base_domains.end(); ++base_domain_itr) {
                    WSInterval base_domain = *base_domain_itr;
                    domains.push_back(WSInterval((1 - base_domain.Min) * current_domain.Min + base_domain.Min * current_domain.Max, 
                        (1 - base_domain.Max) * current_domain.Min + base_domain.Max * current_domain.Max));
                }
                base_domains.clear();
                bezier_iterator.Next();
            }
        }
        else {
            std::vector<WSInterval> base_domains;
            WGNurbsCurve2d::BezierIterator bezier_iterator(nurbs, 0, nurbs->GetPieceCount());
            while (!bezier_iterator.Eof()) {
                WGIntersectHelper2d::RemoveBezierCurveSingularity(nurbs->GetDegree(), bezier_iterator.GetCurrentColtrolPoints(),
                    distance_epsilon, base_domains);
                MergeDomains(base_domains);
                WSInterval current_domain = bezier_iterator.GetCurrentDomain();
                for (auto base_domain_itr = base_domains.begin(); base_domain_itr != base_domains.end(); ++base_domain_itr) {
                    WSInterval base_domain = *base_domain_itr;
                    domains.push_back(WSInterval((1 - base_domain.Min) * current_domain.Min + base_domain.Min * current_domain.Max,
                        (1 - base_domain.Max) * current_domain.Min + base_domain.Max * current_domain.Max));
                }
                base_domains.clear();
                bezier_iterator.Next();
            }
        }
        MergeDomains(domains);
        curves.reserve(domains.size());
        if (domains.size() == 1) {
            const WSInterval& domain = domains.at(0);
            if (domain.Min == nurbs->GetPieceDomain(0).Min && domain.Max == nurbs->GetPieceDomain(nurbs->GetPieceCount() - 1).Max) {
                curves.push_back((WGCurve2d*)curve);
                return;
            }
        }
        int piece_index = 0;
        for (auto domain_itr = domains.begin(); domain_itr != domains.end(); ++domain_itr) {
            const WSInterval& domain = *domain_itr;
            int start_piece_index = -1;
            while (piece_index < nurbs->GetPieceCount()) {
                WSInterval domain2 = nurbs->GetPieceDomain(piece_index);
                if (domain.Min < domain2.Max) {
                    start_piece_index = piece_index;
                    break;
                }
                ++piece_index;
            }
            if (start_piece_index == -1) {
                break;
            }
            int end_piece_index = -1;
            while (piece_index < nurbs->GetPieceCount()) {
                WSInterval domain2 = nurbs->GetPieceDomain(piece_index);
                if (domain.Max <= domain2.Max) {
                    end_piece_index = piece_index;
                    break;
                }
                ++piece_index;
            }
            if (end_piece_index == -1) {
                break;
            }
            curves.push_back(curve->CreateSubCurve(start_piece_index, domain.Min, end_piece_index, domain.Max));
        }
    }
}

void WGCurveIntersecter2d::AddIntersections(std::vector<WGPointCurveIntersection2d>& intersections,
    WGIntersectHelper2d::PointCurveIntersection* base_intersections, int base_intersection_count, int piece_index) {
    if (intersections.size() == 0) {
        intersections.reserve(base_intersection_count);
    }
    for (int i = 0; i < base_intersection_count; ++i) {
        WGPointCurveIntersection2d intersection;
        intersection.PieceIndex = piece_index;
        intersection.T = base_intersections[i].T;
        intersections.push_back(intersection);
    }
}

void WGCurveIntersecter2d::AddIntersections(std::vector<WGPointCurveIntersection2d>& intersections,
    const std::vector<WGIntersectHelper2d::PointCurveIntersection>& base_intersections, const WSInterval& domain, int piece_index) {
    if (intersections.size() == 0) {
        intersections.reserve(base_intersections.size());
    }
    for (auto itr = base_intersections.begin(); itr != base_intersections.end(); ++itr) {
        const WGIntersectHelper2d::PointCurveIntersection& base_intersection = *itr;
        WGPointCurveIntersection2d intersection;
        intersection.PieceIndex = piece_index;
        intersection.T = (1 - base_intersection.T) * domain.Min + base_intersection.T * domain.Max;
        intersections.push_back(intersection);
    }
}

void WGCurveIntersecter2d::AddIntersections(std::vector<WGCurveCurveIntersection2d>& intersections,
    WGIntersectHelper2d::CurveCurveIntersection* base_intersections, int base_intersection_count,
    int piece_index0, int piece_index1, bool exchanged) {
    if (exchanged) {
        if (intersections.size() == 0) {
            intersections.reserve(base_intersection_count);
        }
        for (int i = 0; i < base_intersection_count; ++i) {
            WGCurveCurveIntersection2d intersection;
            intersection.PieceIndex0 = piece_index1;
            intersection.PieceIndex1 = piece_index0;
            intersection.IsOverlap = base_intersections[i].PointCount == 2;
            if (intersection.IsOverlap) {
                intersection.Ts0[0] = base_intersections[i].Ts1[0];
                intersection.Ts1[0] = base_intersections[i].Ts0[0];
                intersection.IsSamples[0] = base_intersections[i].IsSamples[0];
                intersection.Ts0[1] = base_intersections[i].Ts1[1];
                intersection.Ts1[1] = base_intersections[i].Ts0[1];
                intersection.IsSamples[1] = base_intersections[i].IsSamples[1];
            }
            else {
                intersection.Ts0[0] = base_intersections[i].Ts1[0];
                intersection.Ts1[0] = base_intersections[i].Ts0[0];
                intersection.IsSamples[0] = base_intersections[i].IsSamples[0];
            }
            intersections.push_back(intersection);
        }
    }
    else {
        if (intersections.size() == 0) {
            intersections.reserve(base_intersection_count);
        }
        for (int i = 0; i < base_intersection_count; ++i) {
            WGCurveCurveIntersection2d intersection;
            intersection.PieceIndex0 = piece_index0;
            intersection.PieceIndex1 = piece_index1;
            intersection.IsOverlap = base_intersections[i].PointCount == 2;
            if (intersection.IsOverlap) {
                intersection.Ts0[0] = base_intersections[i].Ts0[0];
                intersection.Ts1[0] = base_intersections[i].Ts1[0];
                intersection.IsSamples[0] = base_intersections[i].IsSamples[0];
                intersection.Ts0[1] = base_intersections[i].Ts0[1];
                intersection.Ts1[1] = base_intersections[i].Ts1[1];
                intersection.IsSamples[1] = base_intersections[i].IsSamples[1];
            }
            else {
                intersection.Ts0[0] = base_intersections[i].Ts0[0];
                intersection.Ts1[0] = base_intersections[i].Ts1[0];
                intersection.IsSamples[0] = base_intersections[i].IsSamples[0];
            }
            intersections.push_back(intersection);
        }
    }
}

void WGCurveIntersecter2d::AddIntersections(std::vector<WGCurveCurveIntersection2d>& intersections,
    const std::vector<WGIntersectHelper2d::CurveCurveIntersection>& base_intersections,
    int piece_index0, int piece_index1, bool exchanged) {
    if (exchanged) {
        if (intersections.size() == 0) {
            intersections.reserve(base_intersections.size());
        }
        for (auto itr = base_intersections.begin(); itr != base_intersections.end(); ++itr) {
            const WGIntersectHelper2d::CurveCurveIntersection& base_intersection = *itr;
            WGCurveCurveIntersection2d intersection;
            intersection.PieceIndex0 = piece_index1;
            intersection.PieceIndex1 = piece_index0;
            intersection.IsOverlap = base_intersection.PointCount == 2;
            if (intersection.IsOverlap) {
                intersection.Ts0[0] = base_intersection.Ts1[0];
                intersection.Ts1[0] = base_intersection.Ts0[0];
                intersection.IsSamples[0] = base_intersection.IsSamples[0];
                intersection.Ts0[1] = base_intersection.Ts1[1];
                intersection.Ts1[1] = base_intersection.Ts0[1];
                intersection.IsSamples[1] = base_intersection.IsSamples[1];
            }
            else {
                intersection.Ts0[0] = base_intersection.Ts1[0];
                intersection.Ts1[0] = base_intersection.Ts0[0];
                intersection.IsSamples[0] = base_intersection.IsSamples[0];
            }
            intersections.push_back(intersection);
        }
    }
    else {
        if (intersections.size() == 0) {
            intersections.reserve(base_intersections.size());
        }
        for (auto itr = base_intersections.begin(); itr != base_intersections.end(); ++itr) {
            const WGIntersectHelper2d::CurveCurveIntersection& base_intersection = *itr;
            WGCurveCurveIntersection2d intersection;
            intersection.PieceIndex0 = piece_index0;
            intersection.PieceIndex1 = piece_index1;
            intersection.IsOverlap = base_intersection.PointCount == 2;
            if (intersection.IsOverlap) {
                intersection.Ts0[0] = base_intersection.Ts0[0];
                intersection.Ts1[0] = base_intersection.Ts1[0];
                intersection.IsSamples[0] = base_intersection.IsSamples[0];
                intersection.Ts0[1] = base_intersection.Ts0[1];
                intersection.Ts1[1] = base_intersection.Ts1[1];
                intersection.IsSamples[1] = base_intersection.IsSamples[1];
            }
            else {
                intersection.Ts0[0] = base_intersection.Ts0[0];
                intersection.Ts1[0] = base_intersection.Ts1[0];
                intersection.IsSamples[0] = base_intersection.IsSamples[0];
            }
            intersections.push_back(intersection);
        }
    }
}

void WGCurveIntersecter2d::AddIntersections(std::vector<WGCurveCurveIntersection2d>& intersections,
    const std::vector<WGIntersectHelper2d::CurveCurveIntersection>& base_intersections,
    const WSInterval* domain0, const WSInterval* domain1, int piece_index0, int piece_index1, bool exchanged) {
    if (exchanged) {
        if (intersections.size() == 0) {
            intersections.reserve(base_intersections.size());
        }
        if (domain0) {
            if (domain1) {
                for (auto itr = base_intersections.begin(); itr != base_intersections.end(); ++itr) {
                    const WGIntersectHelper2d::CurveCurveIntersection& base_intersection = *itr;
                    WGCurveCurveIntersection2d intersection;
                    intersection.PieceIndex0 = piece_index1;
                    intersection.PieceIndex1 = piece_index0;
                    intersection.IsOverlap = base_intersection.PointCount == 2;
                    if (intersection.IsOverlap) {
                        intersection.Ts0[0] = (1 - base_intersection.Ts1[0]) * domain1->Min + base_intersection.Ts1[0] * domain1->Max;
                        intersection.Ts1[0] = (1 - base_intersection.Ts0[0]) * domain0->Min + base_intersection.Ts0[0] * domain0->Max;
                        intersection.IsSamples[0] = base_intersection.IsSamples[0];
                        intersection.Ts0[1] = (1 - base_intersection.Ts1[1]) * domain1->Min + base_intersection.Ts1[1] * domain1->Max;
                        intersection.Ts1[1] = (1 - base_intersection.Ts0[1]) * domain0->Min + base_intersection.Ts0[1] * domain0->Max;
                        intersection.IsSamples[1] = base_intersection.IsSamples[1];
                    }
                    else {
                        intersection.Ts0[0] = (1 - base_intersection.Ts1[0]) * domain1->Min + base_intersection.Ts1[0] * domain1->Max;
                        intersection.Ts1[0] = (1 - base_intersection.Ts0[0]) * domain0->Min + base_intersection.Ts0[0] * domain0->Max;
                        intersection.IsSamples[0] = base_intersection.IsSamples[0];
                    }
                    intersections.push_back(intersection);
                }
            }
            else {
                for (auto itr = base_intersections.begin(); itr != base_intersections.end(); ++itr) {
                    const WGIntersectHelper2d::CurveCurveIntersection& base_intersection = *itr;
                    WGCurveCurveIntersection2d intersection;
                    intersection.PieceIndex0 = piece_index1;
                    intersection.PieceIndex1 = piece_index0;
                    intersection.IsOverlap = base_intersection.PointCount == 2;
                    if (intersection.IsOverlap) {
                        intersection.Ts0[0] = base_intersection.Ts1[0];
                        intersection.Ts1[0] = (1 - base_intersection.Ts0[0]) * domain0->Min + base_intersection.Ts0[0] * domain0->Max;
                        intersection.IsSamples[0] = base_intersection.IsSamples[0];
                        intersection.Ts0[1] = base_intersection.Ts1[1];
                        intersection.Ts1[1] = (1 - base_intersection.Ts0[1]) * domain0->Min + base_intersection.Ts0[1] * domain0->Max;
                        intersection.IsSamples[1] = base_intersection.IsSamples[1];
                    }
                    else {
                        intersection.Ts0[0] = base_intersection.Ts1[0];
                        intersection.Ts1[0] = (1 - base_intersection.Ts0[0]) * domain0->Min + base_intersection.Ts0[0] * domain0->Max;
                        intersection.IsSamples[0] = base_intersection.IsSamples[0];
                    }
                    intersections.push_back(intersection);
                }
            }
        }
        else {
            if (domain1) {
                for (auto itr = base_intersections.begin(); itr != base_intersections.end(); ++itr) {
                    const WGIntersectHelper2d::CurveCurveIntersection& base_intersection = *itr;
                    WGCurveCurveIntersection2d intersection;
                    intersection.PieceIndex0 = piece_index1;
                    intersection.PieceIndex1 = piece_index0;
                    intersection.IsOverlap = base_intersection.PointCount == 2;
                    if (intersection.IsOverlap) {
                        intersection.Ts0[0] = (1 - base_intersection.Ts1[0]) * domain1->Min + base_intersection.Ts1[0] * domain1->Max;
                        intersection.Ts1[0] = base_intersection.Ts0[0];
                        intersection.IsSamples[0] = base_intersection.IsSamples[0];
                        intersection.Ts0[1] = (1 - base_intersection.Ts1[1]) * domain1->Min + base_intersection.Ts1[1] * domain1->Max;
                        intersection.Ts1[1] = base_intersection.Ts0[1];
                        intersection.IsSamples[1] = base_intersection.IsSamples[1];
                    }
                    else {
                        intersection.Ts0[0] = (1 - base_intersection.Ts1[0]) * domain1->Min + base_intersection.Ts1[0] * domain1->Max;
                        intersection.Ts1[0] = base_intersection.Ts0[0];
                        intersection.IsSamples[0] = base_intersection.IsSamples[0];
                    }
                    intersections.push_back(intersection);
                }
            }
            else {
                for (auto itr = base_intersections.begin(); itr != base_intersections.end(); ++itr) {
                    const WGIntersectHelper2d::CurveCurveIntersection& base_intersection = *itr;
                    WGCurveCurveIntersection2d intersection;
                    intersection.PieceIndex0 = piece_index1;
                    intersection.PieceIndex1 = piece_index0;
                    intersection.IsOverlap = base_intersection.PointCount == 2;
                    if (intersection.IsOverlap) {
                        intersection.Ts0[0] = base_intersection.Ts1[0];
                        intersection.Ts1[0] = base_intersection.Ts0[0];
                        intersection.IsSamples[0] = base_intersection.IsSamples[0];
                        intersection.Ts0[1] = base_intersection.Ts1[1];
                        intersection.Ts1[1] = base_intersection.Ts0[1];
                        intersection.IsSamples[1] = base_intersection.IsSamples[1];
                    }
                    else {
                        intersection.Ts0[0] = base_intersection.Ts1[0];
                        intersection.Ts1[0] = base_intersection.Ts0[0];
                        intersection.IsSamples[0] = base_intersection.IsSamples[0];
                    }
                    intersections.push_back(intersection);
                }
            }
        }        
    }
    else {
        if (intersections.size() == 0) {
            intersections.reserve(base_intersections.size());
        }
        if (domain0) {
            if (domain1) {
                for (auto itr = base_intersections.begin(); itr != base_intersections.end(); ++itr) {
                    const WGIntersectHelper2d::CurveCurveIntersection& base_intersection = *itr;
                    WGCurveCurveIntersection2d intersection;
                    intersection.PieceIndex0 = piece_index0;
                    intersection.PieceIndex1 = piece_index1;
                    intersection.IsOverlap = base_intersection.PointCount == 2;
                    if (intersection.IsOverlap) {
                        intersection.Ts0[0] = (1 - base_intersection.Ts0[0]) * domain0->Min + base_intersection.Ts0[0] * domain0->Max;
                        intersection.Ts1[0] = (1 - base_intersection.Ts1[0]) * domain1->Min + base_intersection.Ts1[0] * domain1->Max;
                        intersection.IsSamples[0] = base_intersection.IsSamples[0];
                        intersection.Ts0[1] = (1 - base_intersection.Ts0[1]) * domain0->Min + base_intersection.Ts0[1] * domain0->Max;
                        intersection.Ts1[1] = (1 - base_intersection.Ts1[1]) * domain1->Min + base_intersection.Ts1[1] * domain1->Max;
                        intersection.IsSamples[1] = base_intersection.IsSamples[1];
                    }
                    else {
                        intersection.Ts0[0] = (1 - base_intersection.Ts0[0]) * domain0->Min + base_intersection.Ts0[0] * domain0->Max;
                        intersection.Ts1[0] = (1 - base_intersection.Ts1[0]) * domain1->Min + base_intersection.Ts1[0] * domain1->Max;
                        intersection.IsSamples[0] = base_intersection.IsSamples[0];
                    }
                    intersections.push_back(intersection);
                }
            }
            else {
                for (auto itr = base_intersections.begin(); itr != base_intersections.end(); ++itr) {
                    const WGIntersectHelper2d::CurveCurveIntersection& base_intersection = *itr;
                    WGCurveCurveIntersection2d intersection;
                    intersection.PieceIndex0 = piece_index0;
                    intersection.PieceIndex1 = piece_index1;
                    intersection.IsOverlap = base_intersection.PointCount == 2;
                    if (intersection.IsOverlap) {
                        intersection.Ts0[0] = (1 - base_intersection.Ts0[0]) * domain0->Min + base_intersection.Ts0[0] * domain0->Max;
                        intersection.Ts1[0] = base_intersection.Ts1[0];
                        intersection.IsSamples[0] = base_intersection.IsSamples[0];
                        intersection.Ts0[1] = (1 - base_intersection.Ts0[1]) * domain0->Min + base_intersection.Ts0[1] * domain0->Max;
                        intersection.Ts1[1] = base_intersection.Ts1[1];
                        intersection.IsSamples[1] = base_intersection.IsSamples[1];
                    }
                    else {
                        intersection.Ts0[0] = (1 - base_intersection.Ts0[0]) * domain0->Min + base_intersection.Ts0[0] * domain0->Max;
                        intersection.Ts1[0] = base_intersection.Ts1[0];
                        intersection.IsSamples[0] = base_intersection.IsSamples[0];
                    }
                    intersections.push_back(intersection);
                }
            }
        }
        else {
            if (domain1) {
                for (auto itr = base_intersections.begin(); itr != base_intersections.end(); ++itr) {
                    const WGIntersectHelper2d::CurveCurveIntersection& base_intersection = *itr;
                    WGCurveCurveIntersection2d intersection;
                    intersection.PieceIndex0 = piece_index0;
                    intersection.PieceIndex1 = piece_index1;
                    intersection.IsOverlap = base_intersection.PointCount == 2;
                    if (intersection.IsOverlap) {
                        intersection.Ts0[0] = base_intersection.Ts0[0];
                        intersection.Ts1[0] = (1 - base_intersection.Ts1[0]) * domain1->Min + base_intersection.Ts1[0] * domain1->Max;
                        intersection.IsSamples[0] = base_intersection.IsSamples[0];
                        intersection.Ts0[1] = base_intersection.Ts0[1];
                        intersection.Ts1[1] = (1 - base_intersection.Ts1[1]) * domain1->Min + base_intersection.Ts1[1] * domain1->Max;
                        intersection.IsSamples[1] = base_intersection.IsSamples[1];
                    }
                    else {
                        intersection.Ts0[0] = base_intersection.Ts0[0];
                        intersection.Ts1[0] = (1 - base_intersection.Ts1[0]) * domain1->Min + base_intersection.Ts1[0] * domain1->Max;
                        intersection.IsSamples[0] = base_intersection.IsSamples[0];
                    }
                    intersections.push_back(intersection);
                }
            }
            else {
                for (auto itr = base_intersections.begin(); itr != base_intersections.end(); ++itr) {
                    const WGIntersectHelper2d::CurveCurveIntersection& base_intersection = *itr;
                    WGCurveCurveIntersection2d intersection;
                    intersection.PieceIndex0 = piece_index0;
                    intersection.PieceIndex1 = piece_index1;
                    intersection.IsOverlap = base_intersection.PointCount == 2;
                    if (intersection.IsOverlap) {
                        intersection.Ts0[0] = base_intersection.Ts0[0];
                        intersection.Ts1[0] = base_intersection.Ts1[0];
                        intersection.IsSamples[0] = base_intersection.IsSamples[0];
                        intersection.Ts0[1] = base_intersection.Ts0[1];
                        intersection.Ts1[1] = base_intersection.Ts1[1];
                        intersection.IsSamples[1] = base_intersection.IsSamples[1];
                    }
                    else {
                        intersection.Ts0[0] = base_intersection.Ts0[0];
                        intersection.Ts1[0] = base_intersection.Ts1[0];
                        intersection.IsSamples[0] = base_intersection.IsSamples[0];
                    }
                    intersections.push_back(intersection);
                }
            }
        }
    }
}

void WGCurveIntersecter2d::MergeDomains(std::vector<WSInterval>& domains) {
    if (domains.size() == 0) {
        return;
    }
    struct LessDomain {
        bool operator()(const WSInterval& domain0, const WSInterval& domain1) {
            return domain0.Middle() < domain1.Middle();
        }
    };
    std::sort(domains.begin(), domains.end(), LessDomain());
    int n = 1;
    WSInterval* domain0 = &domains.at(0);
    for (int i = 1; i < (int)domains.size(); ++i) {
        const WSInterval& domain1 = domains.at(i);
        if (domain0->Max >= domain1.Min) {
            domain0->Max = domain1.Max;
        }
        else {
            domain0 = &domains.at(n);
            if (n != i) {
                *domain0 = domain1;
            }
            ++n;
        }
    }
    domains.resize(n);
}