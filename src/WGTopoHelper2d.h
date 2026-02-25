#ifndef _WG_TOPO_HELPER_2D_
#define _WG_TOPO_HELPER_2D_

#include "WGPolygon.h"
#include <map>
#include <set>
#include <algorithm>
#include <assert.h>
#include "WGCurveIntersecter2d.h"

class WGTopoHelper2d {
public:
    struct SplitPointKey {
        SplitPointKey() {}
        SplitPointKey(const WGWire2d* wire, int edge_index, int piece_index, double t);
        SplitPointKey& Normalize();
        bool operator<(const SplitPointKey& other) const;
        bool operator==(const SplitPointKey& other) const;
        WGVector2d CalculatePoint() const;

        const WGWire2d* Wire;
        int EdgeIndex;
        int PieceIndex;
        double T;
    };
public:
    enum class SplitPointType : uint8_t {
        Inner,
        Joint,
        StartHeader,
        EndHeader
    };
public:
    struct SplitPoint {
        SplitPoint() {}
        SplitPoint(const SplitPointKey& key);
        bool operator<(const SplitPoint& other) const;
        bool PrevIsLine() const;
        bool NextIsLine() const;
        WGVector2d CalculatePrevD1() const;
        WGVector2d CalculateNextD1() const;

        SplitPointKey Key;
        SplitPointType Type;
        WGVector2d Point;
    };
public:
    struct SplitPointPairKey {
        SplitPointPairKey() {}
        SplitPointPairKey(int point_index0, int point_index1);
        bool operator<(const SplitPointPairKey& other) const;

        int PointIndices[2];
    };
public:
    struct SplitPointPair {
        SplitPointPair() {}
        SplitPointPair(const SplitPointPairKey& key);

        SplitPointPairKey Key;
        bool IsSample;
        int RefCount;
    };
public:
    struct Splitter {
        Splitter() {}
        Splitter(int pair_index0, int pair_index1, bool same_direction, bool is_point);

        int PairIndices[2];
        bool SameDirection;
        bool IsPoint;
    };
public:
    enum class SegmentPosition : uint8_t {
        Undefine,
        Unknown,
        End,
        Inner,
        Outter,
        Coincide
    };
public:
    struct MergedSplitter {
        int MinPointIndices[2];
        int MaxPointIndices[2];
        bool Singular;
        std::vector<int> Splitters;
    };
public:
    struct MergedSplitterPosition {
        MergedSplitterPosition() {}
        MergedSplitterPosition(SegmentPosition prev_position, int prev_confidence, SegmentPosition next_position, int next_confidence);

        SegmentPosition PrevPosition;
        int PrevConfidence;
        SegmentPosition NextPosition;
        int NextConfidence;
    };
public:
    struct IntersectionPart {
        IntersectionPart() {}

        IntersectionPart(int min_point_index, int max_point_index, const MergedSplitterPosition& position);

        int MinPointIndex;
        int MaxPointIndex;
        MergedSplitterPosition Position;
        int PrevSetgmentIndex;
        int NextSetgmentIndex;
    };
public:
    struct Intersection {
        IntersectionPart Part[2];
        std::vector<IntersectionPart> Parts[2];
    };
public:
    struct IntersectionPartIndex {
        IntersectionPartIndex() {}
        IntersectionPartIndex(int intersection_index, int part_index);

        int IntersectionIndex;
        int PartIndex;
    };
public:
    struct Segment {
        Segment() {}
        Segment(const IntersectionPartIndex& start, const IntersectionPartIndex& end);

        IntersectionPartIndex Start;
        IntersectionPartIndex End;
        SegmentPosition Position;
        int Confidence;
        int CoincideSegmentIndex;
        int SelectedFlag;
        SplitPointKey StartKey;
        SplitPointKey EndKey;
        bool IsStartExtend;
        WGVector2d StartExtendPoint;
        bool IsEndExtend;
        WGVector2d EndExtendPoint;
    };
public:
    struct CoincideSegment {
        CoincideSegment() {}
        CoincideSegment(int splitter_index);

        int SplitterIndex;
        std::vector<int> SplitterIndices[2];
        std::vector<int> DirectionFlags[2];
        int Flag; //0: Unknown  1: Same direction or inner  2: Reverse direction or outter  3: Ignore
    };
public:
    struct SplitterContainer {
        std::vector<SplitPoint> SplitPoints[2];
        std::map<SplitPointKey, int> SplitPointMap[2];
        std::vector<SplitPointPair> SplitPointPairs;
        std::map<SplitPointPairKey, int> SplitPointPairMap;
        std::vector<Splitter> Splitters;
        std::vector<MergedSplitter> MergedSplitters;
        std::vector<MergedSplitterPosition> MergedSplitterPositions[2];
        std::vector<Intersection> Intersections;
        std::vector<IntersectionPartIndex> SortedIntersectionPartIndices[2];
        std::vector<Segment> Segments[2];
        std::vector<CoincideSegment> CoincideSegments;
    };
public:
    class SegmentSelector {
    public:
        virtual ~SegmentSelector() {}
        //Set SelectedFlag for Segments
        //0: Not selected  1: Selected  -1 : Selected but needs to be reversed
        virtual void SelectSegments(SplitterContainer* container) = 0;
        //0: Not selected  1: Selected  -1 : Selected but needs to be reversed
        virtual int SelectWire(int index, bool is_inner) = 0;
    };
public:
    struct SplitCoincideSplitterError {
        SplitCoincideSplitterError() {}
        SplitCoincideSplitterError(int split_point_index, int splitter_index);

        int SplitPointIndex;
        int SplitterIndex;
    };
public:
    struct SegmentIndex {
        SegmentIndex() {}
        SegmentIndex(int index, int segment_index);

        int Index;
        int SegIndex;
    };
public:
    struct Loop {
        std::vector<SegmentIndex> SegmentIndices;
        bool Closed;
    };
public:
    class SampleIterator {
    public:
        SampleIterator(const SplitPointKey& start_key, const SplitPointKey& end_key, int sample_count_per_piece);
        void First();
        bool Eof();
        void Next();
        WGVector2d GetCurrentPoint();
    private:
        SplitPointKey m_start_key;
        SplitPointKey m_end_key;
        int m_sample_count_per_piece;
    private:
        int m_state;
        int m_edge_index;
        int m_piece_index;
        double m_t;
        double m_ti;
        double m_dt;
        int m_i;
    };
public:
    static const int HighestConfidence = 100;
    static const int HighConfidence = 70;
    static const int MediumConfidence = 30;
    static const int LowConfidence = 20;
    static const int LowestConfidence = 10;
public:
    static Splitter NewSplitter(SplitterContainer* container, const WGWire2d* wire0, int edge_index0,
        const WGWire2d* wire1, int edge_index1, WGCurveCurveIntersection2d* intersection);
    static void AddSplitters(SplitterContainer* container, std::vector<Splitter>& splitters);
    static void SetSplittersIsPoint(SplitterContainer* container, double distance_epsilon);
    static void AddCoincideSideSplitter(SplitterContainer* container);
    static SplitCoincideSplitterError* SplitCoincideSplitters(SplitterContainer* container, int index, double distance_epsilon,
        WGIntersectHelper2d::IntersectCache* intersect_cache, bool& is_singularity, bool check_only);
    static WGVector2d GetSlightMoveDirectionByError(SplitterContainer* container, int index, SplitCoincideSplitterError* error);
    static void BuildMergedSplitters(SplitterContainer* container);
    static void BuildMergedSplitterPositions(SplitterContainer* container, bool build0, bool build1);
    static void InitializeIntersections(SplitterContainer* container);
    static void MergeIntersections(SplitterContainer* container, int index);
    static void MergeIntersectionParts(SplitterContainer* container, int index);
    static void BuildSortedIntersectionPartIndices(SplitterContainer* container, int index);
    static void BuildSegments(SplitterContainer* container, int index);
    static void CalculateSegmentsPosition(SplitterContainer* container, int index, const WGPolygon* other_polygon,
        double distance_epsilon, WGIntersectHelper2d::IntersectCache* intersect_cache, bool& is_singularity);
    static bool CheckSegment(SplitterContainer* container, int index, const Segment& segment);
    static WGVector2d GetSlightMoveDirectionByErrorSegmentIndex(SplitterContainer* container, int index, int error_segment_index);
    static void BuildCoincideSegments(SplitterContainer* container);
    static bool CheckCoincideSegmentDirectionFlags(const CoincideSegment& coincide_segment);
    static WGVector2d GetSlightMoveDirectionByErrorCoincideSegmentIndex(SplitterContainer* container, int error_coincide_segment_index);
    static void BuildCoincideSegmentFlags(SplitterContainer* container, const WGPolygon* polygon0, const WGPolygon* polygon1,
        double distance_epsilon, WGIntersectHelper2d::IntersectCache* intersect_cache, bool& is_singularity0, bool& is_singularity1);
    static void BuildSegmentsCoincideSegmentIndex(SplitterContainer* container, int index);
    static bool CheckCoincideSegmentIndex(const Segment& segment);
    static void SearchLoops(SplitterContainer* container, std::vector<Loop>& loops);
    static bool CheckLoop(const Loop& loop);
    static void RepairEdgeStartEnds(SplitterContainer* container);
    static void AddPolygonLoops(SplitterContainer* container, std::vector<Loop>& loops, WGPolygon* result_polygon);
    static bool AddPolygonSeparatedLoops(SplitterContainer* container, const WGPolygon* polygon0, const WGPolygon* polygon1,
        bool is_box_separated, SegmentSelector* selector, double distance_epsilon,
        bool& is_singularity0, bool& is_singularity1, WGPolygon* result_polygon);
    static bool CheckCoincideSegmentFlag(const CoincideSegment& coincide_segment);
    static Intersection BuildMinIntersection(SplitterContainer* container, int index);
    static void MergeIntersectionPart(SplitterContainer* container, int index, Intersection& intersection);
private:
    static int AddSplitPoint(SplitterContainer* container, int index, const WGWire2d* wire, int edge_index, int piece_index, double t);
    static int AddSplitPointPair(SplitterContainer* container, const SplitPointPairKey& key, bool is_sample);
    static void AddSplitter(SplitterContainer* container, const Splitter& splitter);
    static int FindSplitter(SplitterContainer* container, int pair_index);
    static bool CheckSplitPointPair(SplitterContainer* container, const SplitPointPair& pair);
    static void SetSplitterIsPoint(SplitterContainer* container, Splitter& splitter, double distance_epsilon);
    static bool IsInner(const SplitPointKey& min_point_key, const SplitPointKey& max_point_key, const SplitPointKey& point_key);
    static bool IsIntersect(const SplitPoint& min_point1, const SplitPoint& max_point1, const SplitPoint& min_point2, const SplitPoint& max_point2);
    static bool MergeSplitPointInterval(SplitterContainer* container, int index,
        int& dst_min_index, int& dst_max_index, int src_min_index, int src_max_index);
    static bool SplitCurve(const WGVector2d& point, const WGCurve2d* curve, int piece_index, const WSInterval& t_domain,
        double distance_epsilon, WGIntersectHelper2d::IntersectCache* intersect_cache, bool& is_singularity, double& t);
    static int CalculateConfidence(double delta, bool is_line, bool downgrade);
    static SegmentPosition CalculatePosition(double c01, double c02, double d01, double d02);
    static Intersection InitializeIntersection(SplitterContainer* container, int merged_splitter_index);
    static const SplitPoint& GetMinSplitPoint(SplitterContainer* container, int index, const IntersectionPartIndex& part_index);
    static const SplitPoint& GetMaxSplitPoint(SplitterContainer* container, int index, const IntersectionPartIndex& part_index);
    static bool MergeIntersection(SplitterContainer* container, int index, Intersection& dst_intersection, const Intersection& src_intersection);
    static int FindMinIntersection(SplitterContainer* container, int index);
    static bool GetNearestOn(const WGVector2d& point, const WGPolygon* polygon, double distance_epsilon,
        WGIntersectHelper2d::IntersectCache* intersect_cache, bool& is_singularity, SplitPointKey& result_key, WGVector2d& result_point);
    static bool SameCoincideSplitter(SplitterContainer* container, int index, int splitter_index0, int splitter_index1);
    static bool MergeCoincideSegment(SplitterContainer* container, CoincideSegment& dst_coincide_segment, CoincideSegment& src_coincide_segment);
    //0: Unknown  1: Inner  2: Outter
    static int CheckCoincideSegmentInner(SplitterContainer* container, int index, const WGPolygon* polygon,
        double distance_epsilon, WGIntersectHelper2d::IntersectCache* intersect_cache, bool& is_singularity, const CoincideSegment& coincide_segment);
    static bool IsInnerSplitter(SplitterContainer* container, int index, const SplitPointKey& point_key, const Splitter& splitter);
    static bool IsInnerCoincideSegment(SplitterContainer* container, int index, const SplitPointKey& point_key, const CoincideSegment& coincide_segment);
    static void SearchLoop(SplitterContainer* container, int index, int segment_index, std::vector<Loop>& loops);
    static void AddEdges(WGWire2d* wire, const SplitPointKey& start_key, const SplitPointKey& end_key, bool reverse);
};

class PolygonPolygonSubtractionSegmentSelector : public WGTopoHelper2d::SegmentSelector {
public:
    virtual void SelectSegments(WGTopoHelper2d::SplitterContainer* container);
    virtual int SelectWire(int index, bool is_inner);
};

class PolygonPolygonIntersectionSegmentSelector : public WGTopoHelper2d::SegmentSelector {
public:
    virtual void SelectSegments(WGTopoHelper2d::SplitterContainer* container);
    virtual int SelectWire(int index, bool is_inner);
};

class PolygonPolygonUnionSegmentSelector : public WGTopoHelper2d::SegmentSelector {
public:
    virtual void SelectSegments(WGTopoHelper2d::SplitterContainer* container);
    virtual int SelectWire(int index, bool is_inner);
};

#endif
