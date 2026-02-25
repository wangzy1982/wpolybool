#include "WGTopoHelper2d.h"
#include "WGPointPolygonPosition2d.h"
#include "WGLine2d.h"

WGTopoHelper2d::SplitPointKey::SplitPointKey(const WGWire2d* wire, int edge_index, int piece_index, double t) :
    Wire(wire),
    EdgeIndex(edge_index),
    PieceIndex(piece_index),
    T(t) {
}

WGTopoHelper2d::SplitPointKey& WGTopoHelper2d::SplitPointKey::Normalize() {
    const WGCurve2d* curve = Wire->GetEdge(EdgeIndex)->GetCurve();
    WSInterval piece_domain = curve->GetPieceDomain(PieceIndex);
    if (T == piece_domain.Max) {
        int next_piece_index = PieceIndex + 1;
        if (next_piece_index < curve->GetPieceCount()) {
            PieceIndex = next_piece_index;
            T = curve->GetPieceDomain(PieceIndex).Min;
        }
        else {
            int next_edge_index = EdgeIndex + 1;
            if (next_edge_index == Wire->GetEdgeCount() && Wire->GetClosed()) {
                next_edge_index = 0;
            }
            if (next_edge_index < Wire->GetEdgeCount()) {
                EdgeIndex = next_edge_index;
                PieceIndex = 0;
                T = Wire->GetEdge(EdgeIndex)->GetCurve()->GetPieceDomain(PieceIndex).Min;
            }
        }
    }
    return *this;
}

bool WGTopoHelper2d::SplitPointKey::operator<(const SplitPointKey& other) const {
    if (Wire < other.Wire) {
        return true;
    }
    if (Wire > other.Wire) {
        return false;
    }
    if (EdgeIndex < other.EdgeIndex) {
        return true;
    }
    if (EdgeIndex > other.EdgeIndex) {
        return false;
    }
    if (PieceIndex < other.PieceIndex) {
        return true;
    }
    if (PieceIndex > other.PieceIndex) {
        return false;
    }
    return T < other.T;
}

bool WGTopoHelper2d::SplitPointKey::operator==(const SplitPointKey& other) const {
    return Wire == other.Wire && EdgeIndex == other.EdgeIndex && PieceIndex == other.PieceIndex && T == other.T;
}

WGVector2d WGTopoHelper2d::SplitPointKey::CalculatePoint() const {
    return Wire->GetEdge(EdgeIndex)->GetCurve()->CalculateD0(PieceIndex, T);
}

WGTopoHelper2d::SplitPoint::SplitPoint(const SplitPointKey& key) :
    Key(key) {
    const WGCurve2d* curve = key.Wire->GetEdge(key.EdgeIndex)->GetCurve();
    WSInterval piece_domain = curve->GetPieceDomain(key.PieceIndex);
    if (key.T == piece_domain.Min) {
        if (key.PieceIndex == 0 && key.EdgeIndex == 0 && !key.Wire->GetClosed()) {
            Type = SplitPointType::StartHeader;
        }
        else {
            Type = SplitPointType::Joint;
        }
    }
    else if (key.T == piece_domain.Max) {
        assert(key.PieceIndex == curve->GetPieceCount() - 1 && key.EdgeIndex == key.Wire->GetEdgeCount() - 1 && !key.Wire->GetClosed());
        Type = SplitPointType::EndHeader;
    }
    else {
        Type = SplitPointType::Inner;
    }
    Point = curve->CalculateD0(key.PieceIndex, key.T);
    if (Type == SplitPointType::Joint) {
        if (key.PieceIndex == 0) {
            int edge_index1 = key.EdgeIndex == 0 ? key.Wire->GetEdgeCount() - 1 : key.EdgeIndex - 1;
            WGVector2d point1 = key.Wire->GetEdge(edge_index1)->GetCurve()->GetEndPoint();
            Point = (Point + point1) * 0.5;
        }
        else {
            int piece_index1 = key.PieceIndex - 1;
            double t1 = curve->GetPieceDomain(piece_index1).Max;
            WGVector2d point1 = curve->CalculateD0(piece_index1, t1);
            Point = (Point + point1) * 0.5;
        }
    }
}

bool WGTopoHelper2d::SplitPoint::operator<(const SplitPoint& other) const {
    return Key < other.Key;
}

bool WGTopoHelper2d::SplitPoint::PrevIsLine() const {
    if (Type == SplitPointType::Joint) {
        int edge_index = Key.EdgeIndex - 1;
        if (edge_index < 0) {
            edge_index = Key.Wire->GetEdgeCount() - 1;
        }
        const WGCurve2d* curve = Key.Wire->GetEdge(edge_index)->GetCurve();
        int piece_index = curve->GetPieceCount() - 1;
        return curve->GetType() == WGLine2d::Type::Instance();
    }
    return Key.Wire->GetEdge(Key.EdgeIndex)->GetCurve()->GetType() == WGLine2d::Type::Instance();
}

bool WGTopoHelper2d::SplitPoint::NextIsLine() const {
    return Key.Wire->GetEdge(Key.EdgeIndex)->GetCurve()->GetType() == WGLine2d::Type::Instance();
}

WGVector2d WGTopoHelper2d::SplitPoint::CalculatePrevD1() const {
    if (Type == SplitPointType::Joint) {
        int edge_index = Key.EdgeIndex - 1;
        if (edge_index < 0) {
            edge_index = Key.Wire->GetEdgeCount() - 1;
        }
        const WGCurve2d* curve = Key.Wire->GetEdge(edge_index)->GetCurve();
        int piece_index = curve->GetPieceCount() - 1;
        return curve->CalculateD1(piece_index, curve->GetPieceDomain(piece_index).Max);
    }
    return Key.Wire->GetEdge(Key.EdgeIndex)->GetCurve()->CalculateD1(Key.PieceIndex, Key.T);
}

WGVector2d WGTopoHelper2d::SplitPoint::CalculateNextD1() const {
    return Key.Wire->GetEdge(Key.EdgeIndex)->GetCurve()->CalculateD1(Key.PieceIndex, Key.T);
}

WGTopoHelper2d::SplitPointPairKey::SplitPointPairKey(int point_index0, int point_index1) {
    PointIndices[0] = point_index0;
    PointIndices[1] = point_index1;
}

bool WGTopoHelper2d::SplitPointPairKey::operator<(const SplitPointPairKey& other) const {
    if (PointIndices[0] < other.PointIndices[0]) {
        return true;
    }
    if (PointIndices[0] > other.PointIndices[0]) {
        return false;
    }
    return PointIndices[1] < other.PointIndices[1];
}

WGTopoHelper2d::SplitPointPair::SplitPointPair(const SplitPointPairKey& key) :
    Key(key),
    IsSample(true),
    RefCount(0) {
}

WGTopoHelper2d::Splitter::Splitter(int pair_index0, int pair_index1, bool same_direction, bool is_point) {
    PairIndices[0] = pair_index0;
    PairIndices[1] = pair_index1;
    SameDirection = same_direction;
    IsPoint = is_point;
}

WGTopoHelper2d::MergedSplitterPosition::MergedSplitterPosition(SegmentPosition prev_position, int prev_confidence, 
    SegmentPosition next_position, int next_confidence) :
    PrevPosition(prev_position),
    PrevConfidence(prev_confidence),
    NextPosition(next_position),
    NextConfidence(next_confidence) {
}

WGTopoHelper2d::IntersectionPart::IntersectionPart(int min_point_index, int max_point_index, const MergedSplitterPosition& position) :
    MinPointIndex(min_point_index),
    MaxPointIndex(max_point_index),
    Position(position) {
}

WGTopoHelper2d::IntersectionPartIndex::IntersectionPartIndex(int intersection_index, int part_index) :
    IntersectionIndex(intersection_index),
    PartIndex(part_index) {
}

WGTopoHelper2d::Segment::Segment(const IntersectionPartIndex& start, const IntersectionPartIndex& end) :
    Start(start),
    End(end) {
}

WGTopoHelper2d::CoincideSegment::CoincideSegment(int splitter_index) :
    SplitterIndex(splitter_index) {
}

int WGTopoHelper2d::AddSplitPoint(SplitterContainer* container, int index, const WGWire2d* wire, int edge_index, int piece_index, double t) {
    SplitPointKey key = SplitPointKey(wire, edge_index, piece_index, t);
    key.Normalize();
    auto split_point_itr = container->SplitPointMap[index].find(key);
    if (split_point_itr != container->SplitPointMap[index].end()) {
        return split_point_itr->second;
    }
    int point_index = (int)container->SplitPoints[index].size();
    container->SplitPoints[index].push_back(SplitPoint(key));
    container->SplitPointMap[index].insert(std::pair<SplitPointKey, int>(key, point_index));
    return point_index;
}

int WGTopoHelper2d::AddSplitPointPair(SplitterContainer* container, const SplitPointPairKey& key, bool is_sample) {
    int pair_index;
    auto pair_itr = container->SplitPointPairMap.find(key);
    if (pair_itr != container->SplitPointPairMap.end()) {
        pair_index = pair_itr->second;
    }
    else {
        pair_index = (int)container->SplitPointPairs.size();
        container->SplitPointPairs.push_back(SplitPointPair(key));
        container->SplitPointPairMap.insert(std::pair<SplitPointPairKey, int>(key, pair_index));
    }
    if (!is_sample) {
        SplitPointPair& pair = container->SplitPointPairs.at(pair_index);
        pair.IsSample = false;
    }
    return pair_index;
}

WGTopoHelper2d::Splitter WGTopoHelper2d::NewSplitter(SplitterContainer* container, const WGWire2d* wire0, int edge_index0,
    const WGWire2d* wire1, int edge_index1, WGCurveCurveIntersection2d* intersection) {
    Splitter splitter;
    if (intersection->IsOverlap) {
        if (intersection->Ts0[0] == intersection->Ts0[1]) {
            if (intersection->Ts1[0] < intersection->Ts1[1]) {
                SplitPointPairKey pair_key0(AddSplitPoint(container, 0, wire0, edge_index0, intersection->PieceIndex0, intersection->Ts0[0]),
                    AddSplitPoint(container, 1, wire1, edge_index1, intersection->PieceIndex1, intersection->Ts1[0]));
                SplitPointPairKey pair_key1(pair_key0.PointIndices[0],
                    AddSplitPoint(container, 1, wire1, edge_index1, intersection->PieceIndex1, intersection->Ts1[1]));
                int pair_index0 = AddSplitPointPair(container, pair_key0, intersection->IsSamples[0]);
                int pair_index1 = AddSplitPointPair(container, pair_key1, intersection->IsSamples[1]);
                splitter = Splitter(pair_index0, pair_index1, true, true);
            }
            else {
                SplitPointPairKey pair_key0(AddSplitPoint(container, 0, wire0, edge_index0, intersection->PieceIndex0, intersection->Ts0[0]),
                    AddSplitPoint(container, 1, wire1, edge_index1, intersection->PieceIndex1, intersection->Ts1[1]));
                SplitPointPairKey pair_key1(pair_key0.PointIndices[0],
                    AddSplitPoint(container, 1, wire1, edge_index1, intersection->PieceIndex1, intersection->Ts1[0]));
                int pair_index0 = AddSplitPointPair(container, pair_key0, intersection->IsSamples[1]);
                int pair_index1 = AddSplitPointPair(container, pair_key1, intersection->IsSamples[0]);
                splitter = Splitter(pair_index0, pair_index1, true, true);
            }
        }
        else if (intersection->Ts1[0] == intersection->Ts1[1]) {
            SplitPointPairKey pair_key0(AddSplitPoint(container, 0, wire0, edge_index0, intersection->PieceIndex0, intersection->Ts0[0]),
                AddSplitPoint(container, 1, wire1, edge_index1, intersection->PieceIndex1, intersection->Ts1[0]));
            SplitPointPairKey pair_key1(AddSplitPoint(container, 0, wire0, edge_index0, intersection->PieceIndex0, intersection->Ts0[1]),
                pair_key0.PointIndices[1]);
            int pair_index0 = AddSplitPointPair(container, pair_key0, intersection->IsSamples[0]);
            int pair_index1 = AddSplitPointPair(container, pair_key1, intersection->IsSamples[1]);
            splitter = Splitter(pair_index0, pair_index1, true, true);
        }
        else {
            SplitPointPairKey pair_key0(AddSplitPoint(container, 0, wire0, edge_index0, intersection->PieceIndex0, intersection->Ts0[0]),
                AddSplitPoint(container, 1, wire1, edge_index1, intersection->PieceIndex1, intersection->Ts1[0]));
            SplitPointPairKey pair_key1(AddSplitPoint(container, 0, wire0, edge_index0, intersection->PieceIndex0, intersection->Ts0[1]),
                AddSplitPoint(container, 1, wire1, edge_index1, intersection->PieceIndex1, intersection->Ts1[1]));
            int pair_index0 = AddSplitPointPair(container, pair_key0, intersection->IsSamples[0]);
            int pair_index1 = AddSplitPointPair(container, pair_key1, intersection->IsSamples[1]);
            splitter = Splitter(pair_index0, pair_index1, intersection->Ts1[0] < intersection->Ts1[1], false);
        }
    }
    else {
        SplitPointPairKey pair_key0(AddSplitPoint(container, 0, wire0, edge_index0, intersection->PieceIndex0, intersection->Ts0[0]),
            AddSplitPoint(container, 1, wire1, edge_index1, intersection->PieceIndex1, intersection->Ts1[0]));
        int pair_index0 = AddSplitPointPair(container, pair_key0, intersection->IsSamples[0]);
        splitter = Splitter(pair_index0, pair_index0, true, true);
    }
    return splitter;
}

void WGTopoHelper2d::AddSplitter(SplitterContainer* container, const Splitter& splitter) {
    if (splitter.PairIndices[0] == splitter.PairIndices[1]) {
        SplitPointPair& pair = container->SplitPointPairs.at(splitter.PairIndices[0]);
        ++pair.RefCount;
        if (pair.RefCount == 1) {
            container->Splitters.push_back(splitter);
        }
    }
    else {
        SplitPointPair& pair0 = container->SplitPointPairs.at(splitter.PairIndices[0]);
        ++pair0.RefCount;
        bool b = true;
        if (pair0.RefCount > 1) {
            for (auto splitter_itr = container->Splitters.begin(); splitter_itr != container->Splitters.end(); ++splitter_itr) {
                if (splitter_itr->PairIndices[0] == splitter_itr->PairIndices[1]) {
                    if (splitter_itr->PairIndices[0] == splitter.PairIndices[0]) {
                        *splitter_itr = splitter;
                        b = false;
                        break;
                    }
                }
            }
        }
        SplitPointPair& pair1 = container->SplitPointPairs.at(splitter.PairIndices[1]);
        ++pair1.RefCount;
        if (pair1.RefCount > 1) {
            for (auto splitter_itr = container->Splitters.begin(); splitter_itr != container->Splitters.end(); ++splitter_itr) {
                if (splitter_itr->PairIndices[0] == splitter_itr->PairIndices[1]) {
                    if (splitter_itr->PairIndices[0] == splitter.PairIndices[1]) {
                        if (b) {
                            *splitter_itr = splitter;
                            b = false;
                        }
                        else {
                            container->Splitters.erase(splitter_itr);
                        }
                        break;
                    }
                }
            }
        }
        if (b) {
            container->Splitters.push_back(splitter);
        }
    }
}

void WGTopoHelper2d::AddSplitters(SplitterContainer* container, std::vector<Splitter>& splitters) {
    if (container->Splitters.capacity() < container->Splitters.size() + splitters.size()) {
        container->Splitters.reserve(container->Splitters.size() + splitters.size());
    }
    if (container->Splitters.size() > 0) {
        for (auto splitter_itr = splitters.begin(); splitter_itr != splitters.end(); ++splitter_itr) {
            AddSplitter(container, *splitter_itr);
        }
    }
    else {
        auto splitter_itr1 = splitters.begin();
        auto splitter_itr2 = splitter_itr1;
        while (splitter_itr2 != splitters.end()) {
            if (splitter_itr2->PairIndices[0] != splitter_itr2->PairIndices[1]) {
                if (splitter_itr1 != splitter_itr2) {
                    Splitter t = *splitter_itr1;
                    *splitter_itr1 = *splitter_itr2;
                    *splitter_itr2 = t;
                }
                ++splitter_itr1;
                ++splitter_itr2;
            }
            else {
                ++splitter_itr2;
            }
        }
        for (auto splitter_itr = splitters.begin(); splitter_itr != splitters.end(); ++splitter_itr) {
            if (splitter_itr->PairIndices[0] == splitter_itr->PairIndices[1]) {
                SplitPointPair& pair = container->SplitPointPairs.at(splitter_itr->PairIndices[0]);
                ++pair.RefCount;
                if (pair.RefCount == 1) {
                    container->Splitters.push_back(*splitter_itr);
                }
            }
            else {
                SplitPointPair& pair0 = container->SplitPointPairs.at(splitter_itr->PairIndices[0]);
                ++pair0.RefCount;
                SplitPointPair& pair1 = container->SplitPointPairs.at(splitter_itr->PairIndices[1]);
                ++pair1.RefCount;
                container->Splitters.push_back(*splitter_itr);
            }
        }
    }
}

int WGTopoHelper2d::FindSplitter(SplitterContainer* container, int pair_index) {
    for (int i = (int)container->Splitters.size() - 1; i >= 0; --i) {
        Splitter& splitter = container->Splitters.at(i);
        if (splitter.PairIndices[0] == pair_index && splitter.PairIndices[1] == pair_index) {
            return i;
        }
    }
    return -1;
}

bool WGTopoHelper2d::CheckSplitPointPair(SplitterContainer* container, const SplitPointPair& pair) {
    if (pair.RefCount == -1) {
        return true;
    }
    const SplitPoint& point0 = container->SplitPoints[0].at(pair.Key.PointIndices[0]);
    const SplitPoint& point1 = container->SplitPoints[1].at(pair.Key.PointIndices[1]);
    switch (point0.Type) {
    case SplitPointType::Joint: {
            switch (point1.Type) {
            case SplitPointType::Joint: {
                    return pair.RefCount == 4;
                }
            case SplitPointType::StartHeader:
            case SplitPointType::EndHeader: {
                    return pair.RefCount == 2;
                }
            default: {
                    assert(point1.Type == SplitPointType::Inner);
                    return pair.RefCount == 2;
                }
            }
        }
    case SplitPointType::StartHeader:
    case SplitPointType::EndHeader: {
            switch (point1.Type) {
            case SplitPointType::Joint: {
                    return pair.RefCount == 2;
                }
            case SplitPointType::StartHeader:
            case SplitPointType::EndHeader: {
                    return pair.RefCount == 1;
                }
            default: {
                    assert(point1.Type == SplitPointType::Inner);
                    return pair.RefCount == 1;
                }
            }
        }
    default: {
            assert(point0.Type == SplitPointType::Inner);
            switch (point1.Type) {
            case SplitPointType::Joint: {
                    return pair.RefCount == 2;
                }
            case SplitPointType::StartHeader:
            case SplitPointType::EndHeader: {
                    return pair.RefCount == 1;
                }
            default: {
                    assert(point1.Type == SplitPointType::Inner);
                    return pair.RefCount == 1;
                }
            }
        }
    }
}

void WGTopoHelper2d::SetSplitterIsPoint(SplitterContainer* container, Splitter& splitter, double distance_epsilon) {
    assert(!splitter.IsPoint);
    const SplitPointPair& pair0 = container->SplitPointPairs.at(splitter.PairIndices[0]);
    const SplitPointPair& pair1 = container->SplitPointPairs.at(splitter.PairIndices[1]);
    const SplitPoint& point00 = container->SplitPoints[0].at(pair0.Key.PointIndices[0]);
    const SplitPoint& point10 = container->SplitPoints[0].at(pair1.Key.PointIndices[0]);
    if ((point00.Point - point10.Point).SqrLength() <= distance_epsilon * distance_epsilon) {
        splitter.IsPoint = true;
        return;
    }
    const SplitPoint& point01 = container->SplitPoints[1].at(pair0.Key.PointIndices[1]);
    const SplitPoint& point11 = container->SplitPoints[1].at(pair1.Key.PointIndices[1]);
    if ((point01.Point - point11.Point).SqrLength() <= distance_epsilon * distance_epsilon) {
        splitter.IsPoint = true;
        return;
    }
}

void WGTopoHelper2d::SetSplittersIsPoint(SplitterContainer* container, double distance_epsilon) {
    for (auto splitter_itr = container->Splitters.begin(); splitter_itr != container->Splitters.end(); ++splitter_itr) {
        if (!splitter_itr->IsPoint) {
            SetSplitterIsPoint(container, *splitter_itr, distance_epsilon);
        }
    }
}

bool WGTopoHelper2d::IsInner(const SplitPointKey& min_point_key, const SplitPointKey& max_point_key, const SplitPointKey& point_key) {
    if (min_point_key.Wire != point_key.Wire) {
        return false;
    }
    if (min_point_key.EdgeIndex == max_point_key.EdgeIndex) {
        bool b = true;
        if (min_point_key.PieceIndex > max_point_key.PieceIndex) {
            b = false;
        }
        else if (min_point_key.PieceIndex == max_point_key.PieceIndex) {
            if (min_point_key.T > max_point_key.T) {
                b = false;
            }
        }
        if (b) {
            if (point_key.EdgeIndex != min_point_key.EdgeIndex) {
                return false;
            }
            if (point_key.PieceIndex < min_point_key.PieceIndex) {
                return false;
            }
            if (point_key.PieceIndex == min_point_key.PieceIndex) {
                if (point_key.T < min_point_key.T) {
                    return false;
                }
            }
            if (point_key.PieceIndex > max_point_key.PieceIndex) {
                return false;
            }
            if (point_key.PieceIndex == max_point_key.PieceIndex) {
                if (point_key.T > max_point_key.T) {
                    return false;
                }
            }
            return true;
        }
        else {
            if (point_key.EdgeIndex != min_point_key.EdgeIndex) {
                return true;
            }
            if (point_key.PieceIndex < min_point_key.PieceIndex) {
                return true;
            }
            if (point_key.PieceIndex == min_point_key.PieceIndex) {
                if (point_key.T < min_point_key.T) {
                    return true;
                }
            }
            if (point_key.PieceIndex > max_point_key.PieceIndex) {
                return true;
            }
            if (point_key.PieceIndex == max_point_key.PieceIndex) {
                if (point_key.T > max_point_key.T) {
                    return true;
                }
            }
            return false;
        }
    }
    if (point_key.EdgeIndex == min_point_key.EdgeIndex) {
        if (point_key.PieceIndex < min_point_key.PieceIndex) {
            return false;
        }
        if (point_key.PieceIndex == min_point_key.PieceIndex) {
            if (point_key.T < min_point_key.T) {
                return false;
            }
        }
        return true;
    }
    if (point_key.EdgeIndex == max_point_key.EdgeIndex) {
        if (point_key.PieceIndex > max_point_key.PieceIndex) {
            return false;
        }
        if (point_key.PieceIndex == max_point_key.PieceIndex) {
            if (point_key.T > max_point_key.T) {
                return false;
            }
        }
        return true;
    }
    if (min_point_key.EdgeIndex < max_point_key.EdgeIndex) {
        return point_key.EdgeIndex > min_point_key.EdgeIndex && point_key.EdgeIndex < max_point_key.EdgeIndex;
    }
    else {
        return point_key.EdgeIndex > min_point_key.EdgeIndex || point_key.EdgeIndex < max_point_key.EdgeIndex;
    }
}

bool WGTopoHelper2d::IsIntersect(const SplitPoint& min_point1, const SplitPoint& max_point1, const SplitPoint& min_point2, const SplitPoint& max_point2) {
    return IsInner(min_point1.Key, max_point1.Key, min_point2.Key) ||
        IsInner(min_point1.Key, max_point1.Key, max_point2.Key) ||
        IsInner(min_point2.Key, max_point2.Key, min_point1.Key) ||
        IsInner(min_point2.Key, max_point2.Key, max_point1.Key);
}

bool WGTopoHelper2d::MergeSplitPointInterval(SplitterContainer* container, int index,
    int& dst_min_index, int& dst_max_index, int src_min_index, int src_max_index) {
    bool b = false;
    const SplitPoint& dst_min_point = container->SplitPoints[index].at(dst_min_index);
    const SplitPoint& dst_max_point = container->SplitPoints[index].at(dst_max_index);
    const SplitPoint& src_min_point = container->SplitPoints[index].at(src_min_index);
    const SplitPoint& src_max_point = container->SplitPoints[index].at(src_max_index);
    if (IsInner(src_min_point.Key, src_max_point.Key, dst_min_point.Key)) {
        dst_min_index = src_min_index;
        b = true;
    }
    if (IsInner(src_min_point.Key, src_max_point.Key, dst_max_point.Key)) {
        dst_max_index = src_max_index;
        b = true;
    }
    if (!b) {
        return IsInner(dst_min_point.Key, dst_max_point.Key, src_min_point.Key);
    }
    return true;
}

void WGTopoHelper2d::AddCoincideSideSplitter(SplitterContainer* container) {
    for (int i = (int)container->Splitters.size() - 1; i >= 0; --i) {
        Splitter& splitter = container->Splitters.at(i);
        if (splitter.IsPoint) {
            continue;
        }
        int pair_index0 = splitter.PairIndices[0];
        int pair_index1 = splitter.PairIndices[1];
        if (FindSplitter(container, pair_index0) == -1) {
            container->Splitters.push_back(Splitter(pair_index0, pair_index0, true, true));
        }
        if (FindSplitter(container, pair_index1) == -1) {
            container->Splitters.push_back(Splitter(pair_index1, pair_index1, true, true));
        }
    }
}

bool WGTopoHelper2d::SplitCurve(const WGVector2d& point, const WGCurve2d* curve, int piece_index, const WSInterval& t_domain,
    double distance_epsilon, WGIntersectHelper2d::IntersectCache* intersect_cache, bool& is_singularity, double& t) {
    is_singularity = false;
    std::vector<WGPointCurveIntersection2d> intersections;
    WGCurveIntersecter2d::Intersect(point, curve, piece_index, 1, distance_epsilon, intersect_cache, is_singularity, intersections);
    if (is_singularity) {
        return false;
    }
    double d = distance_epsilon * 2 + 100;
    for (auto intersection_itr = intersections.begin(); intersection_itr != intersections.end(); ++intersection_itr) {
        WGPointCurveIntersection2d intersection = *intersection_itr;
        WSInterval t2 = intersection.T;
        if (t2.Min < t_domain.Min) {
            t2.Min = t_domain.Min;
        }
        if (t2.Max > t_domain.Max) {
            t2.Max = t_domain.Max;
        }
        if (t2.Min <= t2.Max) {
            double t1 = t2.Middle();
            WGVector2d point1 = curve->CalculateD0(piece_index, t1);
            double d1 = (point1 - point).Length();
            if (d1 < d) {
                d = d1;
                t = t1;
            }
        }
    }
    return d <= distance_epsilon;
}

WGTopoHelper2d::SplitCoincideSplitterError::SplitCoincideSplitterError(int split_point_index, int splitter_index) :
    SplitPointIndex(split_point_index),
    SplitterIndex(splitter_index) {
}

WGTopoHelper2d::SplitCoincideSplitterError* WGTopoHelper2d::SplitCoincideSplitters(SplitterContainer* container, int index, double distance_epsilon,
    WGIntersectHelper2d::IntersectCache* intersect_cache, bool& is_singularity, bool check_only) {
    is_singularity = false;
    for (int i = 0; i < (int)container->SplitPoints[index].size(); ++i) {
        const SplitPoint& split_point = container->SplitPoints[index].at(i);
        for (int j = (int)container->Splitters.size() - 1; j >= 0; --j) {
            Splitter& splitter = container->Splitters.at(j);
            if (splitter.IsPoint) {
                continue;
            }
            const SplitPointPair& pair0 = container->SplitPointPairs.at(splitter.PairIndices[0]);
            const SplitPoint& split_point0 = container->SplitPoints[index].at(pair0.Key.PointIndices[index]);
            if (split_point.Key == split_point0.Key) {
                continue;
            }
            const SplitPointPair& pair1 = container->SplitPointPairs.at(splitter.PairIndices[1]);
            const SplitPoint& split_point1 = container->SplitPoints[index].at(pair1.Key.PointIndices[index]);
            if (split_point.Key == split_point1.Key) {
                continue;
            }
            if (index == 0 || splitter.SameDirection) {
                if (!IsInner(split_point0.Key, split_point1.Key, split_point.Key)) {
                    continue;
                }
            }
            else {
                if (!IsInner(split_point1.Key, split_point0.Key, split_point.Key)) {
                    continue;
                }
            }
            int other_index = index ^ 1;
            const SplitPoint& other_split_point0 = container->SplitPoints[other_index].at(pair0.Key.PointIndices[other_index]);
            const SplitPoint& other_split_point1 = container->SplitPoints[other_index].at(pair1.Key.PointIndices[other_index]);
            const SplitPointKey& min_key = (other_index == 0 || splitter.SameDirection) ? other_split_point0.Key : other_split_point1.Key;
            double t = 0;
            const WGCurve2d* curve = min_key.Wire->GetEdge(min_key.EdgeIndex)->GetCurve();
            WSInterval t_domain;
            t_domain.Min = min_key.T;
            if (other_split_point0.Key.EdgeIndex != other_split_point1.Key.EdgeIndex || other_split_point0.Key.PieceIndex != split_point1.Key.PieceIndex) {
                t_domain.Max = curve->GetPieceDomain(min_key.PieceIndex).Max;
            }
            else {
                t_domain.Max = (other_index == 0 || splitter.SameDirection) ? other_split_point1.Key.T : other_split_point0.Key.T;
            }
            if (check_only || !SplitCurve(split_point.Point, curve, min_key.PieceIndex, t_domain,
                distance_epsilon, intersect_cache, is_singularity, t)) {
                if (is_singularity) {
                    return nullptr;
                }
                return new SplitCoincideSplitterError(i, j);
            }
            int i1 = AddSplitPoint(container, other_index, min_key.Wire, min_key.EdgeIndex, min_key.PieceIndex, t);
            int pair_index = AddSplitPointPair(container, index == 0 ? SplitPointPairKey(i, i1) : SplitPointPairKey(i1, i), true);
            container->SplitPointPairs.at(pair_index).RefCount = -1;
            Splitter splitter1 = Splitter(pair_index, splitter.PairIndices[1], splitter.SameDirection, false);
            splitter.PairIndices[1] = pair_index;
            SetSplitterIsPoint(container, splitter, distance_epsilon);
            SetSplitterIsPoint(container, splitter1, distance_epsilon);
            container->Splitters.push_back(splitter1);
            if (FindSplitter(container, pair_index) == -1) {
                container->Splitters.push_back(Splitter(pair_index, pair_index, true, true));
            }
        }
    }
    return nullptr;
}

WGVector2d WGTopoHelper2d::GetSlightMoveDirectionByError(SplitterContainer* container, int index, SplitCoincideSplitterError* error) {
    const SplitPoint& split_point = container->SplitPoints[index].at(error->SplitPointIndex);
    WGVector2d vt = split_point.CalculatePrevD1() + split_point.CalculateNextD1();
    if (abs(vt.X) < abs(vt.Y)) {
        return WGVector2d(1, 0);
    }
    else {
        return WGVector2d(0, 1);
    }
}

void WGTopoHelper2d::BuildMergedSplitters(SplitterContainer* container) {
    container->MergedSplitters.clear();
    container->MergedSplitters.reserve(container->Splitters.size());
    for (int i = 0; i < container->Splitters.size(); ++i) {
        const Splitter& splitter = container->Splitters.at(i);
        const SplitPointPair& pair0 = container->SplitPointPairs.at(splitter.PairIndices[0]);
        const SplitPointPair& pair1 = container->SplitPointPairs.at(splitter.PairIndices[1]);
        MergedSplitter merged_splitter;
        if (splitter.SameDirection) {
            merged_splitter.MinPointIndices[0] = pair0.Key.PointIndices[0];
            merged_splitter.MaxPointIndices[0] = pair1.Key.PointIndices[0];
            merged_splitter.MinPointIndices[1] = pair0.Key.PointIndices[1];
            merged_splitter.MaxPointIndices[1] = pair1.Key.PointIndices[1];
        }
        else {
            merged_splitter.MinPointIndices[0] = pair0.Key.PointIndices[0];
            merged_splitter.MaxPointIndices[0] = pair1.Key.PointIndices[0];
            merged_splitter.MinPointIndices[1] = pair1.Key.PointIndices[1];
            merged_splitter.MaxPointIndices[1] = pair0.Key.PointIndices[1];
        }
        merged_splitter.Singular = pair0.IsSample || pair1.IsSample ||
            !CheckSplitPointPair(container, pair0) || !CheckSplitPointPair(container, pair1);
        container->MergedSplitters.push_back(merged_splitter);
    }
    int i = 0;
    int j = 0;
    while (j < container->MergedSplitters.size()) {
        bool b = true;
        MergedSplitter& merged_splitter0 = container->MergedSplitters.at(j);
        if (container->Splitters.at(j).IsPoint) {
            for (int k = j + 1; k < container->MergedSplitters.size(); ++k) {
                MergedSplitter& merged_splitter1 = container->MergedSplitters.at(k);
                if (container->Splitters.at(k).IsPoint) {
                    int min_index0 = merged_splitter1.MinPointIndices[0];
                    int max_index0 = merged_splitter1.MaxPointIndices[0];
                    int min_index1 = merged_splitter1.MinPointIndices[1];
                    int max_index1 = merged_splitter1.MaxPointIndices[1];
                    if (MergeSplitPointInterval(container, 0, min_index0, max_index0,
                        merged_splitter0.MinPointIndices[0], merged_splitter0.MaxPointIndices[0]) &&
                        MergeSplitPointInterval(container, 1, min_index1, max_index1,
                            merged_splitter0.MinPointIndices[1], merged_splitter0.MaxPointIndices[1])) {
                        merged_splitter1.MinPointIndices[0] = min_index0;
                        merged_splitter1.MaxPointIndices[0] = max_index0;
                        merged_splitter1.MinPointIndices[1] = min_index1;
                        merged_splitter1.MaxPointIndices[1] = max_index1;
                        if (merged_splitter1.Splitters.size() == 0) {
                            merged_splitter1.Splitters.reserve(10);
                            merged_splitter1.Splitters.push_back(j);
                            merged_splitter1.Splitters.push_back(k);
                        }
                        else {
                            merged_splitter1.Splitters.push_back(j);
                        }
                        merged_splitter1.Singular |= merged_splitter0.Singular;
                        b = false;
                        break;
                    }
                }
            }
        }
        else {
            b = false;
        }
        if (b) {
            if (i != j) {
                container->MergedSplitters.at(i) = std::move(merged_splitter0);
            }
            ++i;
        }
        ++j;
    }
    container->MergedSplitters.resize(i);
}

int WGTopoHelper2d::CalculateConfidence(double delta, bool is_line, bool downgrade) {
    const double delta0 = g_unit_epsilon;
    const double delta1 = 0.001;
    const double delta2 = 0.01;
    const double delta3 = 0.1;
    if (is_zero(delta, delta0)) {
        return 0;
    }
    if (is_zero(delta, delta1)) {
        if (is_line) {
            return downgrade ? LowestConfidence : LowConfidence;
        }
        return downgrade ? 0 : LowestConfidence;
    }
    if (is_zero(delta, delta2)) {
        if (is_line) {
            return downgrade ? MediumConfidence : HighConfidence;
        }
        return downgrade ? LowConfidence : MediumConfidence;
    }
    if (is_zero(delta, delta3)) {
        if (is_line) {
            return downgrade ? HighConfidence : HighestConfidence;
        }
        return downgrade ? MediumConfidence : HighConfidence;
    }
    return downgrade ? HighConfidence : HighestConfidence;
}

WGTopoHelper2d::SegmentPosition WGTopoHelper2d::CalculatePosition(double c01, double c02, double d01, double d02) {
    if (c01 * c02 < 0) {
        if (c01 > 0) {
            return SegmentPosition::Inner;
        }
        else {
            return SegmentPosition::Outter;
        }
    }
    else {
        if (c01 + c02 > 0) {
            if (d02 < d01) {
                return SegmentPosition::Inner;
            }
            else {
                return SegmentPosition::Outter;
            }
        }
        else {
            if (d02 < d01) {
                return SegmentPosition::Outter;
            }
            else {
                return SegmentPosition::Inner;
            }
        }
    }
}

void WGTopoHelper2d::BuildMergedSplitterPositions(SplitterContainer* container, bool build0, bool build1) {
    assert(build0 || build1);
    if (build0) {
        container->MergedSplitterPositions[0].resize(container->MergedSplitters.size());
        for (int i = 0; i < (int)container->MergedSplitters.size(); ++i) {
            const MergedSplitter& merged_splitter = container->MergedSplitters.at(i);
            MergedSplitterPosition& merged_splitter_position = container->MergedSplitterPositions[0].at(i);
            if (merged_splitter.Singular) {
                merged_splitter_position.PrevPosition = SegmentPosition::Unknown;
                merged_splitter_position.PrevConfidence = 0;
                merged_splitter_position.NextPosition = SegmentPosition::Unknown;
                merged_splitter_position.NextConfidence = 0;
            }
            else {
                merged_splitter_position.PrevPosition = SegmentPosition::Undefine;
                merged_splitter_position.PrevConfidence = 0;
                merged_splitter_position.NextPosition = SegmentPosition::Undefine;
                merged_splitter_position.NextConfidence = 0;
            }
        }
    }
    if (build1) {
        container->MergedSplitterPositions[1].resize(container->MergedSplitters.size());
        for (int i = 0; i < (int)container->MergedSplitters.size(); ++i) {
            const MergedSplitter& merged_splitter = container->MergedSplitters.at(i);
            MergedSplitterPosition& merged_splitter_position = container->MergedSplitterPositions[1].at(i);
            if (merged_splitter.Singular) {
                merged_splitter_position.PrevPosition = SegmentPosition::Unknown;
                merged_splitter_position.PrevConfidence = 0;
                merged_splitter_position.NextPosition = SegmentPosition::Unknown;
                merged_splitter_position.NextConfidence = 0;
            }
            else {
                merged_splitter_position.PrevPosition = SegmentPosition::Undefine;
                merged_splitter_position.PrevConfidence = 0;
                merged_splitter_position.NextPosition = SegmentPosition::Undefine;
                merged_splitter_position.NextConfidence = 0;
            }
        }
    }
    for (int i = 0; i < (int)container->MergedSplitters.size(); ++i) {
        const MergedSplitter& merged_splitter = container->MergedSplitters.at(i);
        MergedSplitterPosition* merged_splitter_position0 = build0 ? &container->MergedSplitterPositions[0].at(i) : nullptr;
        MergedSplitterPosition* merged_splitter_position1 = build1 ? &container->MergedSplitterPositions[1].at(i) : nullptr;
        const SplitPoint& min_point00 = container->SplitPoints[0].at(merged_splitter.MinPointIndices[0]);
        const SplitPoint& max_point00 = container->SplitPoints[0].at(merged_splitter.MaxPointIndices[0]);
        const SplitPoint& min_point01 = container->SplitPoints[1].at(merged_splitter.MinPointIndices[1]);
        const SplitPoint& max_point01 = container->SplitPoints[1].at(merged_splitter.MaxPointIndices[1]);
        for (auto splitter_itr = container->Splitters.begin(); splitter_itr != container->Splitters.end(); ++splitter_itr) {
            Splitter& splitter = *splitter_itr;
            if (splitter.IsPoint) {
                continue;
            }
            SplitPointPair pair0 = container->SplitPointPairs.at(splitter.PairIndices[0]);
            SplitPointPair pair1 = container->SplitPointPairs.at(splitter.PairIndices[1]);
            const SplitPoint& min_point10 = container->SplitPoints[0].at(pair0.Key.PointIndices[0]);
            const SplitPoint& max_point10 = container->SplitPoints[0].at(pair1.Key.PointIndices[0]);
            const SplitPoint& min_point11 = container->SplitPoints[1].at(pair0.Key.PointIndices[1]);
            const SplitPoint& max_point11 = container->SplitPoints[1].at(pair1.Key.PointIndices[1]);
            if (splitter.SameDirection) {
                if (build0) {
                    if (IsInner(min_point00.Key, max_point00.Key, min_point10.Key)) {
                        merged_splitter_position0->NextPosition = SegmentPosition::Coincide;
                        merged_splitter_position0->NextConfidence = HighestConfidence;
                    }
                    else if (IsInner(min_point00.Key, max_point00.Key, max_point10.Key)) {
                        merged_splitter_position0->PrevPosition = SegmentPosition::Coincide;
                        merged_splitter_position0->PrevConfidence = HighestConfidence;
                    }
                }
                if (build1) {
                    if (IsInner(min_point01.Key, max_point01.Key, min_point11.Key)) {
                        merged_splitter_position1->NextPosition = SegmentPosition::Coincide;
                        merged_splitter_position1->NextConfidence = HighestConfidence;
                    }
                    else if (IsInner(min_point01.Key, max_point01.Key, max_point11.Key)) {
                        merged_splitter_position1->PrevPosition = SegmentPosition::Coincide;
                        merged_splitter_position1->PrevConfidence = HighestConfidence;
                    }
                }
                /*
                if (IsInner(min_point00.Key, max_point00.Key, min_point10.Key) && IsInner(min_point01.Key, max_point01.Key, min_point11.Key)) {
                    if (build0) {
                        merged_splitter_position0->NextPosition = SegmentPosition::Coincide;
                        merged_splitter_position0->NextConfidence = HighestConfidence;
                    }
                    if (build1) {
                        merged_splitter_position1->NextPosition = SegmentPosition::Coincide;
                        merged_splitter_position1->NextConfidence = HighestConfidence;
                    }
                }
                else if (IsInner(min_point00.Key, max_point00.Key, max_point10.Key) && IsInner(min_point01.Key, max_point01.Key, max_point11.Key)) {
                    if (build0) {
                        merged_splitter_position0->PrevPosition = SegmentPosition::Coincide;
                        merged_splitter_position0->PrevConfidence = HighestConfidence;
                    }
                    if (build1) {
                        merged_splitter_position1->PrevPosition = SegmentPosition::Coincide;
                        merged_splitter_position1->PrevConfidence = HighestConfidence;
                    }
                }
                */
            }
            else {
                if (build0) {
                    if (IsInner(min_point00.Key, max_point00.Key, min_point10.Key)) {
                        merged_splitter_position0->NextPosition = SegmentPosition::Coincide;
                        merged_splitter_position0->NextConfidence = HighestConfidence;
                    }
                    else if (IsInner(min_point00.Key, max_point00.Key, max_point10.Key)) {
                        merged_splitter_position0->PrevPosition = SegmentPosition::Coincide;
                        merged_splitter_position0->PrevConfidence = HighestConfidence;
                    }
                }
                if (build1) {
                    if (IsInner(min_point01.Key, max_point01.Key, min_point11.Key)) {
                        merged_splitter_position1->PrevPosition = SegmentPosition::Coincide;
                        merged_splitter_position1->PrevConfidence = HighestConfidence;
                    }
                    else if (IsInner(min_point01.Key, max_point01.Key, max_point11.Key)) {
                        merged_splitter_position1->NextPosition = SegmentPosition::Coincide;
                        merged_splitter_position1->NextConfidence = HighestConfidence;
                    }
                }
                /*
                if (IsInner(min_point00.Key, max_point00.Key, min_point10.Key) && IsInner(min_point01.Key, max_point01.Key, min_point11.Key)) {
                    if (build0) {
                        merged_splitter_position0->NextPosition = SegmentPosition::Coincide;
                        merged_splitter_position0->NextConfidence = HighestConfidence;
                    }
                    if (build1) {
                        merged_splitter_position1->PrevPosition = SegmentPosition::Coincide;
                        merged_splitter_position1->PrevConfidence = HighestConfidence;
                    }
                }
                else if (IsInner(min_point00.Key, max_point00.Key, max_point10.Key) && IsInner(min_point01.Key, max_point01.Key, max_point11.Key)) {
                    if (build0) {
                        merged_splitter_position0->PrevPosition = SegmentPosition::Coincide;
                        merged_splitter_position0->PrevConfidence = HighestConfidence;
                    }
                    if (build1) {
                        merged_splitter_position1->NextPosition = SegmentPosition::Coincide;
                        merged_splitter_position1->NextConfidence = HighestConfidence;
                    }
                }
                */
            }
        }
        if (build0) {
            if (min_point00.Type == SplitPointType::StartHeader) {
                merged_splitter_position0->PrevPosition = SegmentPosition::End;
                merged_splitter_position0->PrevConfidence = HighestConfidence;
            }
            if (max_point00.Type == SplitPointType::EndHeader) {
                merged_splitter_position0->NextPosition = SegmentPosition::End;
                merged_splitter_position0->NextConfidence = HighestConfidence;
            }
        }
        if (build1) {
            if (min_point01.Type == SplitPointType::StartHeader) {
                merged_splitter_position1->PrevPosition = SegmentPosition::End;
                merged_splitter_position1->PrevConfidence = HighestConfidence;
            }
            if (max_point01.Type == SplitPointType::EndHeader) {
                merged_splitter_position1->NextPosition = SegmentPosition::End;
                merged_splitter_position1->NextConfidence = HighestConfidence;
            }
        }
    }
    for (int i = 0; i < (int)container->MergedSplitters.size(); ++i) {
        const MergedSplitter& merged_splitter = container->MergedSplitters.at(i);
        MergedSplitterPosition* merged_splitter_position0 = build0 ? &container->MergedSplitterPositions[0].at(i) : nullptr;
        MergedSplitterPosition* merged_splitter_position1 = build1 ? &container->MergedSplitterPositions[1].at(i) : nullptr;
        int flag0 = 0;
        int flag1 = 0;
        if (build0) {
            if (merged_splitter_position0->PrevPosition != SegmentPosition::Undefine) {
                flag0 |= 1;
            }
            if (merged_splitter_position0->NextPosition != SegmentPosition::Undefine) {
                flag0 |= 2;
            }
        }
        else {
            flag0 = 3;
        }
        if (build1) {
            if (merged_splitter_position1->PrevPosition != SegmentPosition::Undefine) {
                flag1 |= 1;
            }
            if (merged_splitter_position1->NextPosition != SegmentPosition::Undefine) {
                flag1 |= 2;
            }
        }
        else {
            flag1 = 3;
        }
        if (flag0 == 3 && flag1 == 3) {
            continue;
        }
        const SplitPoint& min_point0 = container->SplitPoints[0].at(merged_splitter.MinPointIndices[0]);
        const SplitPoint& max_point0 = container->SplitPoints[0].at(merged_splitter.MaxPointIndices[0]);
        const SplitPoint& min_point1 = container->SplitPoints[1].at(merged_splitter.MinPointIndices[1]);
        const SplitPoint& max_point1 = container->SplitPoints[1].at(merged_splitter.MaxPointIndices[1]);
        if (merged_splitter.MinPointIndices[0] == merged_splitter.MaxPointIndices[0] &&
            merged_splitter.MinPointIndices[1] == merged_splitter.MaxPointIndices[1]) {
            if (min_point0.Type == SplitPointType::Inner && min_point1.Type == SplitPointType::Inner) {
                WGVector2d vt0 = min_point0.CalculateNextD1();
                WGVector2d vt1 = min_point1.CalculateNextD1();
                if (vt0.Normalize(g_double_epsilon) <= g_double_epsilon || vt1.Normalize(g_double_epsilon) <= g_double_epsilon) {
                    if ((flag0 & 1) == 0) {
                        merged_splitter_position0->PrevPosition = SegmentPosition::Unknown;
                        merged_splitter_position0->PrevConfidence = 0;
                        flag0 |= 1;
                    }
                    if ((flag0 & 2) == 0) {
                        merged_splitter_position0->NextPosition = SegmentPosition::Unknown;
                        merged_splitter_position0->NextConfidence = 0;
                        flag0 |= 2;
                    }
                    if ((flag1 & 1) == 0) {
                        merged_splitter_position1->PrevPosition = SegmentPosition::Unknown;
                        merged_splitter_position1->PrevConfidence = 0;
                        flag1 |= 1;
                    }
                    if ((flag1 & 2) == 0) {
                        merged_splitter_position1->NextPosition = SegmentPosition::Unknown;
                        merged_splitter_position1->NextConfidence = 0;
                        flag1 |= 2;
                    }
                    continue;
                }
                double c = vt0.Cross(vt1);
                int confidence = CalculateConfidence(c, min_point0.NextIsLine() && min_point1.NextIsLine(), false);
                if (confidence == 0) {
                    if ((flag0 & 1) == 0) {
                        merged_splitter_position0->PrevPosition = SegmentPosition::Unknown;
                        merged_splitter_position0->PrevConfidence = confidence;
                        flag0 |= 1;
                    }
                    if ((flag0 & 2) == 0) {
                        merged_splitter_position0->NextPosition = SegmentPosition::Unknown;
                        merged_splitter_position0->NextConfidence = confidence;
                        flag0 |= 2;
                    }
                    if ((flag1 & 1) == 0) {
                        merged_splitter_position1->PrevPosition = SegmentPosition::Unknown;
                        merged_splitter_position1->PrevConfidence = confidence;
                        flag1 |= 1;
                    }
                    if ((flag1 & 2) == 0) {
                        merged_splitter_position1->NextPosition = SegmentPosition::Unknown;
                        merged_splitter_position1->NextConfidence = confidence;
                        flag1 |= 2;
                    }
                }
                else if (c > 0) {
                    if ((flag0 & 1) == 0) {
                        merged_splitter_position0->PrevPosition = SegmentPosition::Inner;
                        merged_splitter_position0->PrevConfidence = confidence;
                        flag0 |= 1;
                    }
                    if ((flag0 & 2) == 0) {
                        merged_splitter_position0->NextPosition = SegmentPosition::Outter;
                        merged_splitter_position0->NextConfidence = confidence;
                        flag0 |= 2;
                    }
                    if ((flag1 & 1) == 0) {
                        merged_splitter_position1->PrevPosition = SegmentPosition::Outter;
                        merged_splitter_position1->PrevConfidence = confidence;
                        flag1 |= 1;
                    }
                    if ((flag1 & 2) == 0) {
                        merged_splitter_position1->NextPosition = SegmentPosition::Inner;
                        merged_splitter_position1->NextConfidence = confidence;
                        flag1 |= 2;
                    }
                }
                else {
                    if ((flag0 & 1) == 0) {
                        merged_splitter_position0->PrevPosition = SegmentPosition::Outter;
                        merged_splitter_position0->PrevConfidence = confidence;
                        flag0 |= 1;
                    }
                    if ((flag0 & 2) == 0) {
                        merged_splitter_position0->NextPosition = SegmentPosition::Inner;
                        merged_splitter_position0->NextConfidence = confidence;
                        flag0 |= 2;
                    }
                    if ((flag1 & 1) == 0) {
                        merged_splitter_position1->PrevPosition = SegmentPosition::Inner;
                        merged_splitter_position1->PrevConfidence = confidence;
                        flag1 |= 1;
                    }
                    if ((flag1 & 2) == 0) {
                        merged_splitter_position1->NextPosition = SegmentPosition::Outter;
                        merged_splitter_position1->NextConfidence = confidence;
                        flag1 |= 2;
                    }
                }
                continue;
            }
        }
        if (flag0 != 3) {
            if (min_point1.Type == SplitPointType::StartHeader || min_point1.Type == SplitPointType::EndHeader) {
                if ((flag0 & 1) == 0) {
                    merged_splitter_position0->PrevPosition = SegmentPosition::Unknown;
                    merged_splitter_position0->PrevConfidence = 0;
                    flag0 |= 1;
                }
                if ((flag0 & 2) == 0) {
                    merged_splitter_position0->NextPosition = SegmentPosition::Unknown;
                    merged_splitter_position0->NextConfidence = 0;
                    flag0 |= 2;
                }
            }
        }
        if (flag1 != 3) {
            if (min_point0.Type == SplitPointType::StartHeader || min_point0.Type == SplitPointType::EndHeader) {
                if ((flag1 & 1) == 0) {
                    merged_splitter_position1->PrevPosition = SegmentPosition::Unknown;
                    merged_splitter_position1->PrevConfidence = 0;
                    flag1 |= 1;
                }
                if ((flag1 & 2) == 0) {
                    merged_splitter_position1->NextPosition = SegmentPosition::Unknown;
                    merged_splitter_position1->NextConfidence = 0;
                    flag1 |= 2;
                }
            }
        }
        if (flag0 == 3 && flag1 == 3) {
            continue;
        }
        WGVector2d vt00 = min_point0.CalculatePrevD1();
        WGVector2d vt01 = max_point0.CalculateNextD1();
        WGVector2d vt10 = min_point1.CalculatePrevD1();
        WGVector2d vt11 = max_point1.CalculateNextD1();
        if (vt00.Normalize(g_double_epsilon) <= g_double_epsilon) {
            if ((flag0 & 1) == 0) {
                merged_splitter_position0->PrevPosition = SegmentPosition::Unknown;
                merged_splitter_position0->PrevConfidence = 0;
                flag0 |= 1;
            }
            if ((flag1 & 1) == 0) {
                merged_splitter_position1->PrevPosition = SegmentPosition::Unknown;
                merged_splitter_position1->PrevConfidence = 0;
                flag1 |= 1;
            }
            if ((flag1 & 2) == 0) {
                merged_splitter_position1->NextPosition = SegmentPosition::Unknown;
                merged_splitter_position1->NextConfidence = 0;
                flag1 |= 2;
            }
        }
        if (vt01.Normalize(g_double_epsilon) <= g_double_epsilon) {
            if ((flag0 & 2) == 0) {
                merged_splitter_position0->NextPosition = SegmentPosition::Unknown;
                merged_splitter_position0->NextConfidence = 0;
                flag0 |= 2;
            }
            if ((flag1 & 1) == 0) {
                merged_splitter_position1->PrevPosition = SegmentPosition::Unknown;
                merged_splitter_position1->PrevConfidence = 0;
                flag1 |= 1;
            }
            if ((flag1 & 2) == 0) {
                merged_splitter_position1->NextPosition = SegmentPosition::Unknown;
                merged_splitter_position1->NextConfidence = 0;
                flag1 |= 2;
            }
        }
        if (vt10.Normalize(g_double_epsilon) <= g_double_epsilon) {
            if ((flag1 & 1) == 0) {
                merged_splitter_position1->PrevPosition = SegmentPosition::Unknown;
                merged_splitter_position1->PrevConfidence = 0;
                flag1 |= 1;
            }
            if ((flag0 & 1) == 0) {
                merged_splitter_position0->PrevPosition = SegmentPosition::Unknown;
                merged_splitter_position0->PrevConfidence = 0;
                flag0 |= 1;
            }
            if ((flag0 & 2) == 0) {
                merged_splitter_position0->NextPosition = SegmentPosition::Unknown;
                merged_splitter_position0->NextConfidence = 0;
                flag0 |= 2;
            }
        }
        if (vt11.Normalize(g_double_epsilon) <= g_double_epsilon) {
            if ((flag1 & 2) == 0) {
                merged_splitter_position1->NextPosition = SegmentPosition::Unknown;
                merged_splitter_position1->NextConfidence = 0;
                flag1 |= 2;
            }
            if ((flag0 & 1) == 0) {
                merged_splitter_position0->PrevPosition = SegmentPosition::Unknown;
                merged_splitter_position0->PrevConfidence = 0;
                flag0 |= 1;
            }
            if ((flag0 & 2) == 0) {
                merged_splitter_position0->NextPosition = SegmentPosition::Unknown;
                merged_splitter_position0->NextConfidence = 0;
                flag0 |= 2;
            }
        }
        if (flag0 == 3 && flag1 == 3) {
            continue;
        }
        bool is_line_00 = min_point0.PrevIsLine();
        bool is_line_01 = max_point0.NextIsLine();
        bool is_line_10 = min_point1.PrevIsLine();
        bool is_line_11 = max_point1.NextIsLine();
        bool downgrade = merged_splitter.Splitters.size() > 1;
        if (flag0 != 3) {
            double d01 = -vt10.Dot(vt11);
            int confidence = CalculateConfidence(1 - d01, is_line_10 && is_line_11, downgrade);
            if (confidence == 0) {
                if ((flag0 & 1) == 0) {
                    merged_splitter_position0->PrevPosition = SegmentPosition::Unknown;
                    merged_splitter_position0->PrevConfidence = 0;
                    flag0 |= 1;
                }
                if ((flag0 & 2) == 0) {
                    merged_splitter_position0->NextPosition = SegmentPosition::Unknown;
                    merged_splitter_position0->NextConfidence = 0;
                    flag0 |= 2;
                }
            }
            else {
                double c01 = -vt10.Cross(vt11);
                if ((flag0 & 1) == 0) {
                    double d02 = vt10.Dot(vt00);
                    int confidence1 = CalculateConfidence(1 - d02, is_line_00 && is_line_10, downgrade);
                    if (confidence1 < confidence) {
                        confidence = confidence1;
                    }
                    if (confidence == 0) {
                        merged_splitter_position0->PrevPosition = SegmentPosition::Unknown;
                        merged_splitter_position0->PrevConfidence = 0;
                        flag0 |= 1;
                    }
                    else {
                        double d12 = -vt11.Dot(vt00);
                        confidence1 = CalculateConfidence(1 - d12, is_line_00 && is_line_11, downgrade);
                        if (confidence1 < confidence) {
                            confidence = confidence1;
                        }
                        if (confidence == 0) {
                            merged_splitter_position0->PrevPosition = SegmentPosition::Unknown;
                            merged_splitter_position0->PrevConfidence = 0;
                            flag0 |= 1;
                        }
                        else {
                            double c02 = vt10.Cross(vt00);
                            merged_splitter_position0->PrevPosition = CalculatePosition(c01, c02, d01, d02);
                            merged_splitter_position0->PrevConfidence = confidence;
                            flag0 |= 1;
                        }
                    }
                }
                if ((flag0 & 2) == 0) {
                    double d02 = -vt10.Dot(vt01);
                    int confidence1 = CalculateConfidence(1 - d02, is_line_01 && is_line_10, downgrade);
                    if (confidence1 < confidence) {
                        confidence = confidence1;
                    }
                    if (confidence == 0) {
                        merged_splitter_position0->NextPosition = SegmentPosition::Unknown;
                        merged_splitter_position0->NextConfidence = 0;
                        flag0 |= 2;
                    }
                    else {
                        double d12 = vt11.Dot(vt01);
                        confidence1 = CalculateConfidence(1 - d12, is_line_01 && is_line_11, downgrade);
                        if (confidence1 < confidence) {
                            confidence = confidence1;
                        }
                        if (confidence == 0) {
                            merged_splitter_position0->NextPosition = SegmentPosition::Unknown;
                            merged_splitter_position0->NextConfidence = 0;
                            flag0 |= 2;
                        }
                        else {
                            double c02 = -vt10.Cross(vt01);
                            merged_splitter_position0->NextPosition = CalculatePosition(c01, c02, d01, d02);
                            merged_splitter_position0->NextConfidence = confidence;
                            flag0 |= 2;
                        }
                    }
                }
            }
        }
        if (flag1 != 3) {
            double d01 = -vt00.Dot(vt01);
            int confidence = CalculateConfidence(1 - d01, is_line_00 && is_line_01, downgrade);
            if (confidence == 0) {
                if ((flag1 & 1) == 0) {
                    merged_splitter_position1->PrevPosition = SegmentPosition::Unknown;
                    merged_splitter_position1->PrevConfidence = 0;
                    flag1 |= 1;
                }
                if ((flag1 & 2) == 0) {
                    merged_splitter_position1->NextPosition = SegmentPosition::Unknown;
                    merged_splitter_position1->NextConfidence = 0;
                    flag1 |= 2;
                }
            }
            else {
                double c01 = -vt00.Cross(vt01);
                if ((flag1 & 1) == 0) {
                    double d02 = vt00.Dot(vt10);
                    int confidence1 = CalculateConfidence(1 - d02, is_line_10 && is_line_00, downgrade);
                    if (confidence1 < confidence) {
                        confidence = confidence1;
                    }
                    if (confidence == 0) {
                        merged_splitter_position1->PrevPosition = SegmentPosition::Unknown;
                        merged_splitter_position1->PrevConfidence = 0;
                        flag1 |= 1;
                    }
                    else {
                        double d12 = -vt01.Dot(vt10);
                        confidence1 = CalculateConfidence(1 - d12, is_line_10 && is_line_01, downgrade);
                        if (confidence1 < confidence) {
                            confidence = confidence1;
                        }
                        if (confidence == 0) {
                            merged_splitter_position1->PrevPosition = SegmentPosition::Unknown;
                            merged_splitter_position1->PrevConfidence = 0;
                            flag1 |= 1;
                        }
                        else {
                            double c02 = vt00.Cross(vt10);
                            merged_splitter_position1->PrevPosition = CalculatePosition(c01, c02, d01, d02);
                            merged_splitter_position1->PrevConfidence = confidence;
                            flag1 |= 1;
                        }
                    }
                }
                if ((flag1 & 2) == 0) {
                    double d02 = -vt00.Dot(vt11);
                    int confidence1 = CalculateConfidence(1 - d02, is_line_11 && is_line_00, downgrade);
                    if (confidence1 < confidence) {
                        confidence = confidence1;
                    }
                    if (confidence == 0) {
                        merged_splitter_position1->NextPosition = SegmentPosition::Unknown;
                        merged_splitter_position1->NextConfidence = 0;
                        flag1 |= 2;
                    }
                    else {
                        double d12 = vt01.Dot(vt11);
                        confidence1 = CalculateConfidence(1 - d12, is_line_11 && is_line_01, downgrade);
                        if (confidence1 < confidence) {
                            confidence = confidence1;
                        }
                        if (confidence == 0) {
                            merged_splitter_position1->NextPosition = SegmentPosition::Unknown;
                            merged_splitter_position1->NextConfidence = 0;
                            flag1 |= 2;
                        }
                        else {
                            double c02 = -vt00.Cross(vt11);
                            merged_splitter_position1->NextPosition = CalculatePosition(c01, c02, d01, d02);
                            merged_splitter_position1->NextConfidence = confidence;
                            flag1 |= 2;
                        }
                    }
                }
            }
        }
    }
}

WGTopoHelper2d::Intersection WGTopoHelper2d::InitializeIntersection(SplitterContainer* container, int merged_splitter_index) {
    const MergedSplitter& merged_splitter = container->MergedSplitters.at(merged_splitter_index);
    Intersection intersection;
    intersection.Part[0] = IntersectionPart(merged_splitter.MinPointIndices[0], merged_splitter.MaxPointIndices[0],
        container->MergedSplitterPositions[0].at(merged_splitter_index));
    intersection.Part[1] = IntersectionPart(merged_splitter.MinPointIndices[1], merged_splitter.MaxPointIndices[1],
        container->MergedSplitterPositions[1].at(merged_splitter_index));
    return intersection;
}

void WGTopoHelper2d::InitializeIntersections(SplitterContainer* container) {
    container->Intersections.reserve(container->MergedSplitters.size());
    for (int i = 0; i < (int)container->MergedSplitters.size(); ++i) {
        container->Intersections.push_back(InitializeIntersection(container, i));
    }
}

const WGTopoHelper2d::SplitPoint& WGTopoHelper2d::GetMinSplitPoint(SplitterContainer* container, int index, const IntersectionPartIndex& part_index) {
    const Intersection& intersection = container->Intersections.at(part_index.IntersectionIndex);
    int point_index;
    if (intersection.Parts[index].size() > 0) {
        point_index = intersection.Parts[index].at(part_index.PartIndex).MinPointIndex;
    }
    else {
        point_index = intersection.Part[index].MinPointIndex;
    }
    return container->SplitPoints[index].at(point_index);
}

const WGTopoHelper2d::SplitPoint& WGTopoHelper2d::GetMaxSplitPoint(SplitterContainer* container, int index, const IntersectionPartIndex& part_index) {
    const Intersection& intersection = container->Intersections.at(part_index.IntersectionIndex);
    int point_index;
    if (intersection.Parts[index].size() > 0) {
        point_index = intersection.Parts[index].at(part_index.PartIndex).MaxPointIndex;
    }
    else {
        point_index = intersection.Part[index].MaxPointIndex;
    }
    return container->SplitPoints[index].at(point_index);
}

struct LessIntersectionPartIndex {
    LessIntersectionPartIndex(WGTopoHelper2d::SplitterContainer* container, int index) : m_container(container), m_index(index) {}

    bool operator()(const WGTopoHelper2d::IntersectionPartIndex& part_index0, const WGTopoHelper2d::IntersectionPartIndex& part_index1) const {
        const WGTopoHelper2d::Intersection& intersection0 = m_container->Intersections.at(part_index0.IntersectionIndex);
        const WGTopoHelper2d::Intersection& intersection1 = m_container->Intersections.at(part_index1.IntersectionIndex);
        int point_index0;
        int point_index1;
        if (intersection0.Parts[m_index].size() > 0) {
            point_index0 = intersection0.Parts[m_index].at(part_index0.PartIndex).MinPointIndex;
        }
        else {
            point_index0 = intersection0.Part[m_index].MinPointIndex;
        }
        if (intersection1.Parts[m_index].size() > 0) {
            point_index1 = intersection1.Parts[m_index].at(part_index1.PartIndex).MinPointIndex;
        }
        else {
            point_index1 = intersection1.Part[m_index].MinPointIndex;
        }
        const WGTopoHelper2d::SplitPoint& point0 = m_container->SplitPoints[m_index].at(point_index0);
        const WGTopoHelper2d::SplitPoint& point1 = m_container->SplitPoints[m_index].at(point_index1);
        return point0.Key < point1.Key;
    }
private:
    WGTopoHelper2d::SplitterContainer* m_container;
    int m_index;
};

WGTopoHelper2d::Intersection WGTopoHelper2d::BuildMinIntersection(SplitterContainer* container, int index) {
    assert(container->Intersections.size() > 0);
    int intersection_count = (int)container->Intersections.size();
    bool* flags = new bool[intersection_count];
    for (int i = 0; i < intersection_count; ++i) {
        flags[i] = true;
    }
    LessIntersectionPartIndex less(container, index);
    IntersectionPartIndex min_part_index(0, 0);
    assert(container->Intersections.at(0).Parts[0].size() == 0 && container->Intersections.at(0).Parts[1].size() == 0);
    for (int i = 1; i < intersection_count; ++i) {
        assert(container->Intersections.at(i).Parts[0].size() == 0 && container->Intersections.at(i).Parts[1].size() == 0);
        IntersectionPartIndex part_index(i, 0);
        if (less(part_index, min_part_index)) {
            min_part_index = part_index;
        }
    }
    flags[min_part_index.IntersectionIndex] = false;
    Intersection intersection = InitializeIntersection(container, min_part_index.IntersectionIndex);
    bool b = true;
    while (b) {
        b = false;
        for (int i = 0; i < intersection_count; ++i) {
            if (flags[i]) {
                if (MergeIntersection(container, index, intersection, container->Intersections.at(i))) {
                    flags[i] = false;
                    b = true;
                    break;
                }
            }
        }
    }
    delete[] flags;
    return intersection;
}

bool WGTopoHelper2d::MergeIntersection(SplitterContainer* container, int index, Intersection& dst_intersection, const Intersection& src_intersection) {
    if (src_intersection.Parts[index].size() > 0) {
        if (dst_intersection.Parts[index].size() > 0) {
            bool b2 = false;
            for (int i0 = 0; i0 < (int)src_intersection.Parts[index].size(); ++i0) {
                const IntersectionPart& part0 = src_intersection.Parts[index].at(i0);
                const SplitPoint& min_point0 = container->SplitPoints[index].at(part0.MinPointIndex);
                const SplitPoint& max_point0 = container->SplitPoints[index].at(part0.MaxPointIndex);
                for (int i1 = 0; i1 < (int)dst_intersection.Parts[index].size(); ++i1) {
                    const IntersectionPart& part1 = dst_intersection.Parts[index].at(i1);
                    const SplitPoint& min_point1 = container->SplitPoints[index].at(part1.MinPointIndex);
                    const SplitPoint& max_point1 = container->SplitPoints[index].at(part1.MaxPointIndex);
                    if (IsIntersect(min_point0, max_point0, min_point1, max_point1)) {
                        b2 = true;
                        break;
                    }
                }
                if (b2) {
                    break;
                }
            }
            if (b2) {
                dst_intersection.Parts[0].insert(dst_intersection.Parts[0].begin(), src_intersection.Parts[0].begin(), src_intersection.Parts[0].end());
                dst_intersection.Parts[1].insert(dst_intersection.Parts[1].begin(), src_intersection.Parts[1].begin(), src_intersection.Parts[1].end());
                return true;
            }
        }
        else {
            const SplitPoint& min_point1 = container->SplitPoints[index].at(dst_intersection.Part[index].MinPointIndex);
            const SplitPoint& max_point1 = container->SplitPoints[index].at(dst_intersection.Part[index].MaxPointIndex);
            bool b2 = false;
            for (int i0 = 0; i0 < (int)src_intersection.Parts[index].size(); ++i0) {
                const IntersectionPart& part0 = src_intersection.Parts[index].at(i0);
                const SplitPoint& min_point0 = container->SplitPoints[index].at(part0.MinPointIndex);
                const SplitPoint& max_point0 = container->SplitPoints[index].at(part0.MaxPointIndex);
                if (IsIntersect(min_point0, max_point0, min_point1, max_point1)) {
                    b2 = true;
                    break;
                }
            }
            if (b2) {
                dst_intersection.Parts[0].reserve(16);
                dst_intersection.Parts[0].insert(dst_intersection.Parts[0].begin(), src_intersection.Parts[0].begin(), src_intersection.Parts[0].end());
                dst_intersection.Parts[0].push_back(dst_intersection.Part[0]);
                dst_intersection.Parts[1].reserve(16);
                dst_intersection.Parts[1].insert(dst_intersection.Parts[1].begin(), src_intersection.Parts[1].begin(), src_intersection.Parts[1].end());
                dst_intersection.Parts[1].push_back(dst_intersection.Part[1]);
                return true;
            }
        }
    }
    else {
        if (dst_intersection.Parts[index].size() > 0) {
            const SplitPoint& min_point0 = container->SplitPoints[index].at(src_intersection.Part[index].MinPointIndex);
            const SplitPoint& max_point0 = container->SplitPoints[index].at(src_intersection.Part[index].MaxPointIndex);
            bool b2 = false;
            for (int i1 = 0; i1 < (int)dst_intersection.Parts[index].size(); ++i1) {
                const IntersectionPart& part1 = dst_intersection.Parts[index].at(i1);
                const SplitPoint& min_point1 = container->SplitPoints[index].at(part1.MinPointIndex);
                const SplitPoint& max_point1 = container->SplitPoints[index].at(part1.MaxPointIndex);
                if (IsIntersect(min_point0, max_point0, min_point1, max_point1)) {
                    b2 = true;
                    break;
                }
            }
            if (b2) {
                dst_intersection.Parts[0].push_back(src_intersection.Part[0]);
                dst_intersection.Parts[1].push_back(src_intersection.Part[1]);
                return true;
            }
        }
        else {
            const SplitPoint& min_point0 = container->SplitPoints[index].at(src_intersection.Part[index].MinPointIndex);
            const SplitPoint& max_point0 = container->SplitPoints[index].at(src_intersection.Part[index].MaxPointIndex);
            const SplitPoint& min_point1 = container->SplitPoints[index].at(dst_intersection.Part[index].MinPointIndex);
            const SplitPoint& max_point1 = container->SplitPoints[index].at(dst_intersection.Part[index].MaxPointIndex);
            if (IsIntersect(min_point0, max_point0, min_point1, max_point1)) {
                dst_intersection.Parts[0].reserve(16);
                dst_intersection.Parts[0].push_back(src_intersection.Part[0]);
                dst_intersection.Parts[0].push_back(dst_intersection.Part[0]);
                dst_intersection.Parts[1].reserve(16);
                dst_intersection.Parts[1].push_back(src_intersection.Part[1]);
                dst_intersection.Parts[1].push_back(dst_intersection.Part[1]);
                return true;
            }
        }
    }
    return false;
}

void WGTopoHelper2d::MergeIntersections(SplitterContainer* container, int index) {
    int i = 0;
    for (int j = 0; j < (int)container->Intersections.size(); ++j) {
        bool b = true;
        Intersection& intersection0 = container->Intersections.at(j);
        for (int k = j + 1; k < (int)container->Intersections.size(); ++k) {
            Intersection& intersection1 = container->Intersections.at(k);
            if (MergeIntersection(container, index, intersection1, intersection0)) {
                b = false;
                break;
            }
        }
        if (b) {
            if (i != j) {
                container->Intersections.at(i) = std::move(intersection0);
            }
            ++i;
        }
    }
    container->Intersections.resize(i);
}

void WGTopoHelper2d::MergeIntersectionPart(SplitterContainer* container, int index, Intersection& intersection) {
    if (intersection.Parts[index].size() > 1) {
        int i = 0;
        for (int j = 0; j < (int)intersection.Parts[index].size(); ++j) {
            bool b = true;
            IntersectionPart& part0 = intersection.Parts[index].at(j);
            for (int k = j + 1; k < (int)intersection.Parts[index].size(); ++k) {
                IntersectionPart& part1 = intersection.Parts[index].at(k);
                int min_index1 = part1.MinPointIndex;
                int max_index1 = part1.MaxPointIndex;
                if (MergeSplitPointInterval(container, index, min_index1, max_index1, part0.MinPointIndex, part0.MaxPointIndex)) {
                    if (min_index1 == part0.MinPointIndex) {
                        part1.MinPointIndex = min_index1;
                        part1.Position.PrevPosition = part0.Position.PrevPosition;
                        part1.Position.PrevConfidence = part0.Position.PrevConfidence;
                    }
                    if (max_index1 == part0.MaxPointIndex) {
                        part1.MaxPointIndex = max_index1;
                        part1.Position.NextPosition = part0.Position.NextPosition;
                        part1.Position.NextConfidence = part0.Position.NextConfidence;
                    }
                    b = false;
                    break;
                }
            }
            if (b) {
                if (j != i) {
                    intersection.Parts[index].at(i) = std::move(part0);
                }
                ++i;
            }
        }
        intersection.Parts[index].resize(i);
    }
}

void WGTopoHelper2d::MergeIntersectionParts(SplitterContainer* container, int index) {
    for (auto intersection_itr = container->Intersections.begin(); intersection_itr != container->Intersections.end(); ++intersection_itr) {
        Intersection& intersection = *intersection_itr;
        MergeIntersectionPart(container, index, intersection);
    }
}

int WGTopoHelper2d::FindMinIntersection(SplitterContainer* container, int index) {
    LessIntersectionPartIndex less(container, index);
    IntersectionPartIndex min_part_index(0, 0);
    for (int i = 0; i < (int)container->Intersections.size(); ++i) {
        const Intersection& intersection = container->Intersections.at(i);
        if (intersection.Parts[index].size() > 0) {
            for (int j = 0; j < (int)intersection.Parts[index].size(); ++j) {
                IntersectionPartIndex part_index(i, j);
                if (less(part_index, min_part_index)) {
                    min_part_index = part_index;
                }
            }
        }
        else {
            IntersectionPartIndex part_index(i, 0);
            if (less(part_index, min_part_index)) {
                min_part_index = part_index;
            }
        }
    }
    return min_part_index.IntersectionIndex;
}

void WGTopoHelper2d::BuildSortedIntersectionPartIndices(SplitterContainer* container, int index) {
    int count = 0;
    for (auto intersection_itr = container->Intersections.begin(); intersection_itr != container->Intersections.end(); ++intersection_itr) {
        const Intersection& intersection = *intersection_itr;
        if (intersection.Parts[index].size() > 0) {
            count += (int)intersection.Parts[index].size();
        }
        else {
            ++count;
        }
    }
    container->SortedIntersectionPartIndices[index].reserve(count);
    for (int i = 0; i < (int)container->Intersections.size(); ++i) {
        const Intersection& intersection = container->Intersections.at(i);
        if (intersection.Parts[index].size() > 0) {
            for (int j = 0; j < (int)intersection.Parts[index].size(); ++j) {
                container->SortedIntersectionPartIndices[index].push_back(IntersectionPartIndex(i, j));
            }
        }
        else {
            container->SortedIntersectionPartIndices[index].push_back(IntersectionPartIndex(i, 0));
        }
    }
    std::sort(container->SortedIntersectionPartIndices[index].begin(), container->SortedIntersectionPartIndices[index].end(),
        LessIntersectionPartIndex(container, index));
}

WGTopoHelper2d::SampleIterator::SampleIterator(const SplitPointKey& start_key, const SplitPointKey& end_key, int sample_count_per_piece) :
    m_start_key(start_key),
    m_end_key(end_key),
    m_sample_count_per_piece(sample_count_per_piece) {
    First();
}

void WGTopoHelper2d::SampleIterator::First() {
    m_state = 0;
    m_edge_index = m_start_key.EdgeIndex;
    m_piece_index = m_start_key.PieceIndex;
    m_t = m_start_key.T;
    Next();
}

bool WGTopoHelper2d::SampleIterator::Eof() {
    return m_state == 10;
}

void WGTopoHelper2d::SampleIterator::Next() {
    assert(!Eof());
    if (m_state == 1) {
        ++m_i;
        if (m_i == m_sample_count_per_piece) {
            m_state = 10;
            return;
        }
        m_ti += m_dt;
    }
    else if (m_state == 2) {
        ++m_i;
        if (m_i == m_sample_count_per_piece) {
            const WGWire2d* wire = m_start_key.Wire;
            const WGCurve2d* curve = wire->GetEdge(m_edge_index)->GetCurve();
            ++m_piece_index;
            if (m_piece_index == curve->GetPieceCount()) {
                ++m_edge_index;
                if (m_edge_index == wire->GetEdgeCount()) {
                    m_edge_index = 0;
                }
                m_piece_index = 0;
            }
            m_t = wire->GetEdge(m_edge_index)->GetCurve()->GetPieceDomain(m_piece_index).Min;
            m_state = 0;
        }
        else {
            m_ti += m_dt;
        }
    }
    while (m_state == 0) {
        m_state = 2;
        if (m_edge_index == m_end_key.EdgeIndex && m_piece_index == m_end_key.PieceIndex) {
            if (m_t < m_end_key.T) {
                double dt = (m_end_key.T - m_t) / (m_sample_count_per_piece + 1);
                if (dt > 0) {
                    const WGWire2d* wire = m_start_key.Wire;
                    const WGCurve2d* curve = wire->GetEdge(m_edge_index)->GetCurve();
                    m_ti = m_t + dt;
                    m_dt = dt;
                    m_i = 0;
                    m_state = 1;
                }
                else {
                    m_state = 10;
                }
            }
        }
        if (m_state == 2) {
            const WGWire2d* wire = m_start_key.Wire;
            const WGCurve2d* curve = wire->GetEdge(m_edge_index)->GetCurve();
            WSInterval domain = curve->GetPieceDomain(m_piece_index);
            double dt = (domain.Max - m_t) / m_sample_count_per_piece;
            if (dt > 0) {
                m_ti = m_t + dt;
                m_dt = dt;
                m_i = 0;
            }
            else {
                m_state = 0;
            }
        }
    }
}

WGVector2d WGTopoHelper2d::SampleIterator::GetCurrentPoint() {
    assert(!Eof());
    const WGWire2d* wire = m_start_key.Wire;
    const WGCurve2d* curve = wire->GetEdge(m_edge_index)->GetCurve();
    return curve->CalculateD0(m_piece_index, m_ti);
}

bool WGTopoHelper2d::GetNearestOn(const WGVector2d& point, const WGPolygon* polygon, double distance_epsilon,
    WGIntersectHelper2d::IntersectCache* intersect_cache, bool& is_singularity, SplitPointKey& result_key, WGVector2d& result_point) {
    is_singularity = false;
    bool exist = false;
    double min_sqr = distance_epsilon * distance_epsilon * 100 + 100;
    std::vector<WGPointCurveIntersection2d> intersections;
    intersections.reserve(16);
    for (int i = 0; i < polygon->GetLoopCount(); ++i) {
        const WGLoop2d* loop = polygon->GetLoop(i);
        const WGWire2d* wire = loop->GetWire();
        for (int j = 0; j < wire->GetEdgeCount(); ++j) {
            const WGCurve2d* curve = wire->GetEdge(j)->GetCurve();
            WGCurveIntersecter2d::Intersect(point, curve, 0, curve->GetPieceCount(), distance_epsilon,
                intersect_cache, is_singularity, intersections);
            if (is_singularity) {
                return false;
            }
            for (auto intersection_itr = intersections.begin(); intersection_itr != intersections.end(); ++intersection_itr) {
                const WGPointCurveIntersection2d& intersection = *intersection_itr;
                double t = intersection.T;
                WGVector2d pt = curve->CalculateD0(intersection.PieceIndex, t);
                double d = (pt - point).SqrLength();
                if (d < min_sqr) {
                    min_sqr = d;
                    result_key.Wire = wire;
                    result_key.EdgeIndex = j;
                    result_key.PieceIndex = intersection.PieceIndex;
                    result_key.T = t;
                    result_point = pt;
                    exist = true;
                }
            }
            intersections.clear();
        }
    }
    return exist;
}

void WGTopoHelper2d::BuildSegments(SplitterContainer* container, int index) {
    container->Segments[index].reserve(container->SortedIntersectionPartIndices[index].size() * 2);
    int i = 0;
    int n = (int)container->SortedIntersectionPartIndices[index].size();
    while (i < n) {
        IntersectionPartIndex* part_index = &container->SortedIntersectionPartIndices[index].at(i);
        Intersection* intersection = &container->Intersections.at(part_index->IntersectionIndex);
        IntersectionPart* part;
        if (intersection->Parts[index].size() > 0) {
            part = &intersection->Parts[index].at(part_index->PartIndex);
        }
        else {
            part = &intersection->Part[index];
        }
        const SplitPoint& min_point = container->SplitPoints[index].at(part->MinPointIndex);
        const WGWire2d* wire = min_point.Key.Wire;
        int m = i + 1;
        while (m < n) {
            IntersectionPartIndex* part_index1 = &container->SortedIntersectionPartIndices[index].at(m);
            Intersection* intersection1 = &container->Intersections.at(part_index1->IntersectionIndex);
            IntersectionPart* part1;
            if (intersection1->Parts[index].size() > 0) {
                part1 = &intersection1->Parts[index].at(part_index1->PartIndex);
            }
            else {
                part1 = &intersection1->Part[index];
            }
            const SplitPoint& min_point1 = container->SplitPoints[index].at(part1->MinPointIndex);
            if (min_point1.Key.Wire != wire) {
                break;
            }
            ++m;
        }
        --m;
        for (int j = i; j < m; ++j) {
            container->Segments[index].push_back(
                Segment(container->SortedIntersectionPartIndices[index].at(j), container->SortedIntersectionPartIndices[index].at(j + 1)));
        }
        if (wire->GetClosed()) {
            container->Segments[index].push_back(
                Segment(container->SortedIntersectionPartIndices[index].at(m), container->SortedIntersectionPartIndices[index].at(i)));
        }
        i = m + 1;
    }
}

void WGTopoHelper2d::CalculateSegmentsPosition(SplitterContainer* container, int index, const WGPolygon* other_polygon,
    double distance_epsilon, WGIntersectHelper2d::IntersectCache* intersect_cache, bool& is_singularity) {
    is_singularity = false;
    for (int i = 0; i < (int)container->Segments[index].size(); ++i) {
        Segment& segment = container->Segments[index].at(i);
        SegmentPosition prev_position;
        int prev_confidence;
        SegmentPosition next_position;
        int next_confidence;
        Intersection* intersection0 = &container->Intersections.at(segment.Start.IntersectionIndex);
        IntersectionPart* part0;
        if (intersection0->Parts[index].size() > 0) {
            part0 = &intersection0->Parts[index].at(segment.Start.PartIndex);
            if (intersection0->Parts[index ^ 1].size() > 1) {
                prev_position = SegmentPosition::Unknown;
                prev_confidence = 0;
            }
            else {
                prev_position = part0->Position.NextPosition;
                prev_confidence = part0->Position.NextConfidence;
            }
        }
        else {
            part0 = &intersection0->Part[index];
            prev_position = part0->Position.NextPosition;
            prev_confidence = part0->Position.NextConfidence;
        }
        Intersection* intersection1 = &container->Intersections.at(segment.End.IntersectionIndex);
        IntersectionPart* part1;
        if (intersection1->Parts[index].size() > 0) {
            part1 = &intersection1->Parts[index].at(segment.End.PartIndex);
            if (intersection1->Parts[index ^ 1].size() > 1) {
                next_position = SegmentPosition::Unknown;
                next_confidence = 0;
            }
            else {
                next_position = part1->Position.PrevPosition;
                next_confidence = part1->Position.PrevConfidence;
            }
        }
        else {
            part1 = &intersection1->Part[index];
            next_position = part1->Position.PrevPosition;
            next_confidence = part1->Position.PrevConfidence;
        }
        part0->NextSetgmentIndex = i;
        part1->PrevSetgmentIndex = i;
        if (prev_position == next_position) {
            segment.Position = prev_position;
            segment.Confidence = prev_confidence + next_confidence;
        }
        else {
            segment.Confidence = prev_confidence - next_confidence;
            if (segment.Confidence > 0) {
                segment.Position = prev_position;
            }
            else {
                segment.Position = next_position;
                segment.Confidence = -segment.Confidence;
            }
        }
        if (segment.Confidence < HighConfidence && other_polygon) {
            const SplitPoint& start_point = GetMaxSplitPoint(container, index, segment.Start);
            const SplitPoint& end_point = GetMinSplitPoint(container, index, segment.End);
            int sample_count_per_pieces[10] = {
                2, 3, 5, 7, 11, 13, 17, 19, 23, 29
            };
            for (int k = 0; k < 10; ++k) {
                SampleIterator iterator2(start_point.Key, end_point.Key, sample_count_per_pieces[k]);
                while (!iterator2.Eof()) {
                    WGPointPolygonPosition2d::Result position = WGPointPolygonPosition2d::Execute(
                        iterator2.GetCurrentPoint(), other_polygon, distance_epsilon);
                    if (position == WGPointPolygonPosition2d::Result::Singularity) {
                        is_singularity = true;
                        return;
                    }
                    if (position == WGPointPolygonPosition2d::Result::Inner) {
                        segment.Position = SegmentPosition::Inner;
                        segment.Confidence = HighestConfidence;
                        break;
                    }
                    if (position == WGPointPolygonPosition2d::Result::Outter) {
                        segment.Position = SegmentPosition::Outter;
                        segment.Confidence = HighestConfidence;
                        break;
                    }
                    iterator2.Next();
                }
                if (segment.Confidence >= HighConfidence) {
                    break;
                }
            }
            if (segment.Confidence < HighConfidence) {
                int total_count = 0;
                int inner_count = 0;
                int outter_count = 0;
                double distance_epsilon2 = distance_epsilon * 2;
                SampleIterator iterator2(start_point.Key, end_point.Key, 29);
                while (!iterator2.Eof()) {
                    WGVector2d point = iterator2.GetCurrentPoint();
                    SplitPointKey nearest_key;
                    WGVector2d nearest_point;
                    if (GetNearestOn(point, other_polygon, distance_epsilon2, intersect_cache, is_singularity, nearest_key, nearest_point)) {
                        if (is_singularity) {
                            return;
                        }
                        const WGCurve2d* curve = nearest_key.Wire->GetEdge(nearest_key.EdgeIndex)->GetCurve();
                        WGVector2d vt = curve->CalculateD1(nearest_key.PieceIndex, nearest_key.T);
                        double d = vt.Cross(point - nearest_point);
                        if (d > 0) {
                            ++inner_count;
                        }
                        else if (d < 0) {
                            ++outter_count;
                        }
                        ++total_count;
                    }
                    iterator2.Next();
                }
                if (total_count > 0) {
                    if (inner_count > outter_count) {
                        if ((double)outter_count / total_count <= 0.3) {
                            segment.Position = SegmentPosition::Inner;
                            segment.Confidence = HighestConfidence;
                        }
                    }
                    else if (inner_count < outter_count) {
                        if ((double)inner_count / total_count <= 0.3) {
                            segment.Position = SegmentPosition::Outter;
                            segment.Confidence = HighestConfidence;
                        }
                    }
                }
            }
        }
    }
}

bool WGTopoHelper2d::CheckSegment(SplitterContainer* container, int index, const Segment& segment) {
    if (segment.Confidence < HighConfidence) {
        return false;
    }
    return true;
}

WGVector2d WGTopoHelper2d::GetSlightMoveDirectionByErrorSegmentIndex(SplitterContainer* container, int index, int error_segment_index) {
    const Segment& segment = container->Segments[index].at(error_segment_index);
    const SplitPoint& split_point = GetMaxSplitPoint(container, index, segment.Start);
    WGVector2d vt = split_point.CalculateNextD1();
    if (abs(vt.X) < abs(vt.Y)) {
        return WGVector2d(1, 0);
    }
    else {
        return WGVector2d(0, 1);
    }
}

bool WGTopoHelper2d::SameCoincideSplitter(SplitterContainer* container, int index, int splitter_index0, int splitter_index1) {
    const Splitter& splitter0 = container->Splitters.at(splitter_index0);
    const SplitPointPair& pair00 = container->SplitPointPairs.at(splitter0.PairIndices[0]);
    const SplitPointPair& pair01 = container->SplitPointPairs.at(splitter0.PairIndices[1]);
    const Splitter& splitter1 = container->Splitters.at(splitter_index1);
    const SplitPointPair& pair10 = container->SplitPointPairs.at(splitter1.PairIndices[0]);
    const SplitPointPair& pair11 = container->SplitPointPairs.at(splitter1.PairIndices[1]);
    if (index == 0) {
        return pair00.Key.PointIndices[0] == pair10.Key.PointIndices[0] && pair01.Key.PointIndices[0] == pair11.Key.PointIndices[0];
    }
    else if (splitter0.SameDirection == splitter1.SameDirection) {
        return pair00.Key.PointIndices[1] == pair10.Key.PointIndices[1] && pair01.Key.PointIndices[1] == pair11.Key.PointIndices[1];
    }
    else {
        return pair00.Key.PointIndices[1] == pair11.Key.PointIndices[1] && pair01.Key.PointIndices[1] == pair10.Key.PointIndices[1];
    }
}

bool WGTopoHelper2d::MergeCoincideSegment(SplitterContainer* container, CoincideSegment& dst_coincide_segment, CoincideSegment& src_coincide_segment) {
    if (dst_coincide_segment.SplitterIndex != -1) {
        if (src_coincide_segment.SplitterIndex != -1) {
            if (SameCoincideSplitter(container, 0, dst_coincide_segment.SplitterIndex, src_coincide_segment.SplitterIndex)) {
                dst_coincide_segment.SplitterIndices[0].push_back(dst_coincide_segment.SplitterIndex);
                dst_coincide_segment.SplitterIndices[1].push_back(dst_coincide_segment.SplitterIndex);
                dst_coincide_segment.SplitterIndices[1].push_back(src_coincide_segment.SplitterIndex);
                dst_coincide_segment.SplitterIndex = -1;
                return true;
            }
            else if (SameCoincideSplitter(container, 1, dst_coincide_segment.SplitterIndex, src_coincide_segment.SplitterIndex)) {
                dst_coincide_segment.SplitterIndices[0].push_back(dst_coincide_segment.SplitterIndex);
                dst_coincide_segment.SplitterIndices[0].push_back(src_coincide_segment.SplitterIndex);
                dst_coincide_segment.SplitterIndices[1].push_back(dst_coincide_segment.SplitterIndex);
                dst_coincide_segment.SplitterIndex = -1;
                return true;
            }
            return false;
        }
        else {
            bool b0 = false;
            bool b1 = false;
            for (auto src_itr = src_coincide_segment.SplitterIndices[0].begin(); src_itr != src_coincide_segment.SplitterIndices[0].end(); ++src_itr) {
                if (SameCoincideSplitter(container, 0, dst_coincide_segment.SplitterIndex, *src_itr)) {
                    b0 = true;
                    break;
                }
            }
            for (auto src_itr = src_coincide_segment.SplitterIndices[1].begin(); src_itr != src_coincide_segment.SplitterIndices[1].end(); ++src_itr) {
                if (SameCoincideSplitter(container, 1, dst_coincide_segment.SplitterIndex, *src_itr)) {
                    b1 = true;
                    break;
                }
            }
            if (b0 || b1) {
                dst_coincide_segment.SplitterIndices[0] = std::move(src_coincide_segment.SplitterIndices[0]);
                dst_coincide_segment.SplitterIndices[1] = std::move(src_coincide_segment.SplitterIndices[1]);
                if (!b0) {
                    dst_coincide_segment.SplitterIndices[0].push_back(dst_coincide_segment.SplitterIndex);
                }
                if (!b1) {
                    dst_coincide_segment.SplitterIndices[1].push_back(dst_coincide_segment.SplitterIndex);
                }
                dst_coincide_segment.SplitterIndex = -1;
                return true;
            }
            return false;
        }
    }
    else {
        if (src_coincide_segment.SplitterIndex != -1) {
            bool b0 = false;
            bool b1 = false;
            for (auto dst_itr = dst_coincide_segment.SplitterIndices[0].begin(); dst_itr != dst_coincide_segment.SplitterIndices[0].end(); ++dst_itr) {
                if (SameCoincideSplitter(container, 0, src_coincide_segment.SplitterIndex, *dst_itr)) {
                    b0 = true;
                    break;
                }
            }
            for (auto dst_itr = dst_coincide_segment.SplitterIndices[1].begin(); dst_itr != dst_coincide_segment.SplitterIndices[1].end(); ++dst_itr) {
                if (SameCoincideSplitter(container, 1, src_coincide_segment.SplitterIndex, *dst_itr)) {
                    b1 = true;
                    break;
                }
            }
            if (b0 || b1) {
                if (!b0) {
                    dst_coincide_segment.SplitterIndices[0].push_back(src_coincide_segment.SplitterIndex);
                }
                if (!b1) {
                    dst_coincide_segment.SplitterIndices[1].push_back(src_coincide_segment.SplitterIndex);
                }
                dst_coincide_segment.SplitterIndex = -1;
                return true;
            }
            return false;
        }
        else {
            bool b = false;
            for (auto src_itr = src_coincide_segment.SplitterIndices[0].begin(); src_itr != src_coincide_segment.SplitterIndices[0].end(); ++src_itr) {
                for (auto dst_itr = dst_coincide_segment.SplitterIndices[0].begin(); dst_itr != dst_coincide_segment.SplitterIndices[0].end(); ++dst_itr) {
                    if (SameCoincideSplitter(container, 0, *src_itr, *dst_itr)) {
                        b = true;
                        break;
                    }
                }
            }
            if (!b) {
                for (auto src_itr = src_coincide_segment.SplitterIndices[1].begin(); src_itr != src_coincide_segment.SplitterIndices[1].end(); ++src_itr) {
                    for (auto dst_itr = dst_coincide_segment.SplitterIndices[1].begin(); dst_itr != dst_coincide_segment.SplitterIndices[1].end(); ++dst_itr) {
                        if (SameCoincideSplitter(container, 1, *src_itr, *dst_itr)) {
                            b = true;
                            break;
                        }
                    }
                }
            }
            if (b) {
                int n = (int)dst_coincide_segment.SplitterIndices[0].size();
                for (auto src_itr = src_coincide_segment.SplitterIndices[0].begin(); src_itr != src_coincide_segment.SplitterIndices[0].end(); ++src_itr) {
                    bool b2 = true;
                    for (int i = 0; i < n; ++i) {
                        if (SameCoincideSplitter(container, 0, *src_itr, dst_coincide_segment.SplitterIndices[0].at(i))) {
                            b2 = false;
                            break;
                        }
                    }
                    if (b2) {
                        dst_coincide_segment.SplitterIndices[0].push_back(*src_itr);
                    }
                }
                n = (int)dst_coincide_segment.SplitterIndices[1].size();
                for (auto src_itr = src_coincide_segment.SplitterIndices[1].begin(); src_itr != src_coincide_segment.SplitterIndices[1].end(); ++src_itr) {
                    bool b2 = true;
                    for (int i = 0; i < n; ++i) {
                        if (SameCoincideSplitter(container, 1, *src_itr, dst_coincide_segment.SplitterIndices[1].at(i))) {
                            b2 = false;
                            break;
                        }
                    }
                    if (b2) {
                        dst_coincide_segment.SplitterIndices[1].push_back(*src_itr);
                    }
                }
                return true;
            }
            return false;
        }
    }
}

void WGTopoHelper2d::BuildCoincideSegments(SplitterContainer* container) {
    container->CoincideSegments.reserve(container->Splitters.size());
    for (int i = 0; i < (int)container->Splitters.size(); ++i) {
        const Splitter& splitter = container->Splitters.at(i);
        if (splitter.IsPoint) {
            continue;
        }
        container->CoincideSegments.push_back(CoincideSegment(i));
    }
    int i = 0;
    int j = 0;
    while (j < (int)container->CoincideSegments.size()) {
        bool b = true;
        CoincideSegment& coincide_segment = container->CoincideSegments.at(j);
        for (int k = j + 1; k < (int)container->CoincideSegments.size(); ++k) {
            if (MergeCoincideSegment(container, container->CoincideSegments.at(k), coincide_segment)) {
                b = false;
                break;
            }
        }
        if (b) {
            if (i != j) {
                container->CoincideSegments.at(i) = std::move(coincide_segment);
            }
            ++i;
        }
        ++j;
    }
    container->CoincideSegments.resize(i);
    for (auto coincide_segment_itr = container->CoincideSegments.begin(); coincide_segment_itr != container->CoincideSegments.end(); ++coincide_segment_itr) {
        CoincideSegment& coincide_segment = *coincide_segment_itr;
        if (coincide_segment.SplitterIndex == -1) {
            int n0 = (int)coincide_segment.SplitterIndices[0].size();
            int n1 = (int)coincide_segment.SplitterIndices[1].size();
            coincide_segment.DirectionFlags[0].reserve(n0);
            for (int i = 0; i < n0; ++i) {
                coincide_segment.DirectionFlags[0].push_back(0);
            }
            coincide_segment.DirectionFlags[1].reserve(n1);
            for (int i = 0; i < n1; ++i) {
                coincide_segment.DirectionFlags[1].push_back(0);
            }
            coincide_segment.DirectionFlags[0].at(0) = 1;
            for (int i = 0; i < n0; ++i) {
                int splitter_index0 = coincide_segment.SplitterIndices[0].at(i);
                if (coincide_segment.DirectionFlags[0].at(i) == 0) {
                    for (int j = 0; j < n1; ++j) {
                        if (coincide_segment.DirectionFlags[1].at(j) != 0) {
                            int splitter_index1 = coincide_segment.SplitterIndices[1].at(j);
                            if (SameCoincideSplitter(container, 1, splitter_index0, splitter_index1)) {
                                if (container->Splitters.at(splitter_index0).SameDirection) {
                                    coincide_segment.DirectionFlags[0].at(i) = coincide_segment.DirectionFlags[1].at(j);
                                }
                                else {
                                    coincide_segment.DirectionFlags[0].at(i) = -coincide_segment.DirectionFlags[1].at(j);
                                }
                                break;
                            }
                        }
                    }
                }
                for (int j = 0; j < n1; ++j) {
                    if (coincide_segment.DirectionFlags[1].at(j) == 0) {
                        int splitter_index1 = coincide_segment.SplitterIndices[1].at(j);
                        if (SameCoincideSplitter(container, 0, splitter_index0, splitter_index1)) {
                            if (container->Splitters.at(splitter_index1).SameDirection) {
                                coincide_segment.DirectionFlags[1].at(j) = coincide_segment.DirectionFlags[0].at(i);
                            }
                            else {
                                coincide_segment.DirectionFlags[1].at(j) = -coincide_segment.DirectionFlags[0].at(i);
                            }
                        }
                    }
                }
            }
        }
    }
}

bool WGTopoHelper2d::CheckCoincideSegmentDirectionFlags(const CoincideSegment& coincide_segment) {
    for (int i = 0; i <= 1; ++i) {
        int a = 0;
        int b = 0;
        for (auto itr = coincide_segment.DirectionFlags[i].begin(); itr != coincide_segment.DirectionFlags[i].end(); ++itr) {
            int n = *itr;
            if (n == 0) {
                return false;
            }
            if (n > 0) {
                ++a;
            }
            else {
                ++b;
            }
        }
        if (abs(a - b) > 1) {
            return false;
        }
    }
    return true;
}

bool WGTopoHelper2d::CheckCoincideSegmentFlag(const CoincideSegment& coincide_segment) {
    return coincide_segment.Flag != 0;
}

WGVector2d WGTopoHelper2d::GetSlightMoveDirectionByErrorCoincideSegmentIndex(SplitterContainer* container, int error_coincide_segment_index) {
    const CoincideSegment& coincide_segment = container->CoincideSegments.at(error_coincide_segment_index);
    int splitter_index = coincide_segment.SplitterIndex;
    if (splitter_index == -1) {
        splitter_index = coincide_segment.SplitterIndices[0].at(0);
    }
    const Splitter& splitter = container->Splitters.at(splitter_index);
    const SplitPointPair& pair = container->SplitPointPairs.at(splitter.PairIndices[0]);
    const SplitPoint& split_point = container->SplitPoints[0].at(pair.Key.PointIndices[0]);
    WGVector2d vt = split_point.CalculateNextD1();
    if (abs(vt.X) < abs(vt.Y)) {
        return WGVector2d(1, 0);
    }
    else {
        return WGVector2d(0, 1);
    }
}

//0: Unknown  1: Inner  2: Outter
int WGTopoHelper2d::CheckCoincideSegmentInner(SplitterContainer* container, int index, const WGPolygon* polygon,
    double distance_epsilon, WGIntersectHelper2d::IntersectCache* intersect_cache, bool& is_singularity, const CoincideSegment& coincide_segment) {
    is_singularity = false;
    SplitterContainer container2;
    container2.SplitPoints[index] = container->SplitPoints[index];
    container2.Intersections.resize(2);
    Intersection& intersection0 = container2.Intersections.at(0);
    Intersection& intersection1 = container2.Intersections.at(1);
    int n = (int)coincide_segment.SplitterIndices[index].size();
    for (int i = 0; i < n; ++i) {
        int splitter_index = coincide_segment.SplitterIndices[index].at(i);
        for (int j = 0; j < i; ++j) {
            if (splitter_index == (int)coincide_segment.SplitterIndices[index].at(j)) {
                splitter_index = -1;
                break;
            }
        }
        if (splitter_index == -1) {
            continue;
        }
        const Splitter& splitter = container->Splitters.at(splitter_index);
        int split_point_index0;
        int split_point_index1;
        if (index == 0 || splitter.SameDirection) {
            split_point_index0 = container->SplitPointPairs.at(splitter.PairIndices[0]).Key.PointIndices[index];
            split_point_index1 = container->SplitPointPairs.at(splitter.PairIndices[1]).Key.PointIndices[index];
        }
        else {
            split_point_index0 = container->SplitPointPairs.at(splitter.PairIndices[1]).Key.PointIndices[index];
            split_point_index1 = container->SplitPointPairs.at(splitter.PairIndices[0]).Key.PointIndices[index];
        }
        if (coincide_segment.DirectionFlags[index].at(i) == 1) {
            bool b = true;
            for (int j = 0; j < (int)intersection0.Parts[index].size(); ++j) {
                IntersectionPart& intersection_part = intersection0.Parts[index].at(j);
                if (intersection_part.MinPointIndex == split_point_index0) {
                    intersection_part.Position.PrevPosition = SegmentPosition::Unknown;
                    intersection_part.Position.PrevConfidence = 0;
                    intersection_part.Position.NextPosition = SegmentPosition::Unknown;
                    intersection_part.Position.NextConfidence = 0;
                    b = false;
                }
            }
            if (b) {
                intersection0.Parts[index].push_back(IntersectionPart(split_point_index0, split_point_index0,
                    MergedSplitterPosition(SegmentPosition::Inner, HighestConfidence, SegmentPosition::Coincide, HighestConfidence)));
            }
            b = true;
            for (int j = 0; j < (int)intersection1.Parts[index].size(); ++j) {
                IntersectionPart& intersection_part = intersection1.Parts[index].at(j);
                if (intersection_part.MinPointIndex == split_point_index1) {
                    intersection_part.Position.PrevPosition = SegmentPosition::Unknown;
                    intersection_part.Position.PrevConfidence = 0;
                    intersection_part.Position.NextPosition = SegmentPosition::Unknown;
                    intersection_part.Position.NextConfidence = 0;
                    b = false;
                }
            }
            if (b) {
                intersection1.Parts[index].push_back(IntersectionPart(split_point_index1, split_point_index1,
                    MergedSplitterPosition(SegmentPosition::Coincide, HighestConfidence, SegmentPosition::Inner, HighestConfidence)));
            }
        }
        else {
            bool b = true;
            for (int j = 0; j < (int)intersection0.Parts[index].size(); ++j) {
                IntersectionPart& intersection_part = intersection0.Parts[index].at(j);
                if (intersection_part.MinPointIndex == split_point_index1) {
                    intersection_part.Position.PrevPosition = SegmentPosition::Unknown;
                    intersection_part.Position.PrevConfidence = 0;
                    intersection_part.Position.NextPosition = SegmentPosition::Unknown;
                    intersection_part.Position.NextConfidence = 0;
                    b = false;
                }
            }
            if (b) {
                intersection0.Parts[index].push_back(IntersectionPart(split_point_index1, split_point_index1,
                    MergedSplitterPosition(SegmentPosition::Coincide, HighestConfidence, SegmentPosition::Inner, HighestConfidence)));
            }
            b = true;
            for (int j = 0; j < (int)intersection1.Parts[index].size(); ++j) {
                IntersectionPart& intersection_part = intersection1.Parts[index].at(j);
                if (intersection_part.MinPointIndex == split_point_index0) {
                    intersection_part.Position.PrevPosition = SegmentPosition::Unknown;
                    intersection_part.Position.PrevConfidence = 0;
                    intersection_part.Position.NextPosition = SegmentPosition::Unknown;
                    intersection_part.Position.NextConfidence = 0;
                    b = false;
                }
            }
            if (b) {
                intersection1.Parts[index].push_back(IntersectionPart(split_point_index0, split_point_index0,
                    MergedSplitterPosition(SegmentPosition::Inner, HighestConfidence, SegmentPosition::Coincide, HighestConfidence)));
            }
        }
    }
    BuildSortedIntersectionPartIndices(&container2, index);
    BuildSegments(&container2, index);
    CalculateSegmentsPosition(&container2, index, nullptr, 0, intersect_cache, is_singularity);
    if (is_singularity) {
        return 0;
    }
    for (auto segment_itr = container2.Segments[index].begin(); segment_itr != container2.Segments[index].end(); ++segment_itr) {
        Segment& segment = *segment_itr;
        if (segment.Position == SegmentPosition::Unknown) {
            return 0;
        }
        segment.SelectedFlag = segment.Position == SegmentPosition::Inner ? 1 : 0;
    }
    std::vector<Loop> loops;
    SearchLoops(&container2, loops);
    int error_loop_index = -1;
    for (int i = 0; i < (int)loops.size(); ++i) {
        if (!WGTopoHelper2d::CheckLoop(loops.at(i))) {
            error_loop_index = i;
            break;
        }
    }
    if (error_loop_index != -1) {
        return 0;
    }
    WGTopoHelper2d::RepairEdgeStartEnds(&container2);
    WGPolygon* polygon2 = new WGPolygon();
    WGTopoHelper2d::AddPolygonLoops(&container2, loops, polygon2);
    for (int i = 0; i < polygon->GetLoopCount(); ++i) {
        const WGLoop2d* loop = polygon->GetLoop(i);
        const WGWire2d* wire = loop->GetWire();
        bool b = true;
        for (int j = 0; j < n; ++j) {
            int splitter_index = coincide_segment.SplitterIndices[index].at(i);
            const Splitter& splitter = container->Splitters.at(splitter_index);
            int split_point_index0 = container->SplitPointPairs.at(splitter.PairIndices[0]).Key.PointIndices[index];
            const SplitPoint& split_point = container->SplitPoints[index].at(split_point_index0);
            if (split_point.Key.Wire == wire) {
                b = false;
                break;
            }
        }
        if (b) {
            polygon2->AddLoop(loop->Clone());
        }
    }
    const Splitter& splitter = container->Splitters.at(0);
    int split_point_index0;
    int split_point_index1;
    if (index == 0 || splitter.SameDirection) {
        split_point_index0 = container->SplitPointPairs.at(splitter.PairIndices[0]).Key.PointIndices[index];
        split_point_index1 = container->SplitPointPairs.at(splitter.PairIndices[1]).Key.PointIndices[index];
    }
    else {
        split_point_index0 = container->SplitPointPairs.at(splitter.PairIndices[1]).Key.PointIndices[index];
        split_point_index1 = container->SplitPointPairs.at(splitter.PairIndices[0]).Key.PointIndices[index];
    }
    SampleIterator iterator(container->SplitPoints[index].at(split_point_index0).Key, container->SplitPoints[index].at(split_point_index1).Key, 10);
    iterator.First();
    while (!iterator.Eof()) {
        WGPointPolygonPosition2d::Result position = WGPointPolygonPosition2d::Execute(iterator.GetCurrentPoint(), polygon2, distance_epsilon);
        if (position == WGPointPolygonPosition2d::Result::Singularity) {
            is_singularity = true;
            return 0;
        }
        if (position == WGPointPolygonPosition2d::Result::Inner) {
            delete polygon2;
            return 1;
        }
        if (position == WGPointPolygonPosition2d::Result::Outter) {
            delete polygon2;
            return 2;
        }
        iterator.Next();
    }
    delete polygon2;
    return 0;
}

void WGTopoHelper2d::BuildCoincideSegmentFlags(SplitterContainer* container, const WGPolygon* polygon0, const WGPolygon* polygon1,
    double distance_epsilon, WGIntersectHelper2d::IntersectCache* intersect_cache, bool& is_singularity0, bool& is_singularity1) {
    is_singularity0 = false;
    is_singularity1 = false;
    for (auto coincide_segment_itr = container->CoincideSegments.begin(); coincide_segment_itr != container->CoincideSegments.end(); ++coincide_segment_itr) {
        CoincideSegment& coincide_segment = *coincide_segment_itr;
        if (coincide_segment.SplitterIndex != -1) {
            coincide_segment.Flag = container->Splitters.at(coincide_segment.SplitterIndex).SameDirection ? 1 : 2;
        }
        else {
            int n0 = (int)coincide_segment.SplitterIndices[0].size();
            int n1 = (int)coincide_segment.SplitterIndices[1].size();
            if ((n0 & 1) == 0) {
                if ((n1 & 1) == 0) {
                    coincide_segment.Flag = 3;
                }
                else {
                    coincide_segment.Flag = CheckCoincideSegmentInner(container, 0, polygon0, distance_epsilon,
                        intersect_cache, is_singularity0, coincide_segment);
                    if (is_singularity0) {
                        return;
                    }
                }
            }
            else {
                if ((n1 & 1) == 0) {
                    coincide_segment.Flag = CheckCoincideSegmentInner(container, 1, polygon1, distance_epsilon,
                        intersect_cache, is_singularity1, coincide_segment);
                    if (is_singularity1) {
                        return;
                    }
                }
                else {
                    int s0 = 0;
                    for (auto itr = coincide_segment.DirectionFlags[0].begin(); itr != coincide_segment.DirectionFlags[0].end(); ++itr) {
                        if (*itr > 0) {
                            ++s0;
                        }
                    }
                    int s1 = 0;
                    for (auto itr = coincide_segment.DirectionFlags[1].begin(); itr != coincide_segment.DirectionFlags[1].end(); ++itr) {
                        if (*itr > 0) {
                            ++s1;
                        }
                    }
                    if (s0 > n0 - s0) {
                        if (s1 > n1 - s1) {
                            coincide_segment.Flag = 1;
                        }
                        else {
                            coincide_segment.Flag = 2;
                        }
                    }
                    else {
                        if (s1 > n1 - s1) {
                            coincide_segment.Flag = 2;
                        }
                        else {
                            coincide_segment.Flag = 1;
                        }
                    }
                }
            }
        }
    }
}

bool WGTopoHelper2d::IsInnerSplitter(SplitterContainer* container, int index, const SplitPointKey& point_key, const Splitter& splitter) {
    const SplitPointPair& pair0 = container->SplitPointPairs.at(splitter.PairIndices[0]);
    const SplitPointPair& pair1 = container->SplitPointPairs.at(splitter.PairIndices[1]);
    if (index == 0) {
        return IsInner(container->SplitPoints[index].at(pair0.Key.PointIndices[0]).Key,
            container->SplitPoints[index].at(pair1.Key.PointIndices[0]).Key, point_key);
    }
    else if (splitter.SameDirection) {
        return IsInner(container->SplitPoints[index].at(pair0.Key.PointIndices[1]).Key,
            container->SplitPoints[index].at(pair1.Key.PointIndices[1]).Key, point_key);
    }
    else {
        return IsInner(container->SplitPoints[index].at(pair1.Key.PointIndices[1]).Key,
            container->SplitPoints[index].at(pair0.Key.PointIndices[1]).Key, point_key);
    }
}

bool WGTopoHelper2d::IsInnerCoincideSegment(SplitterContainer* container, int index, const SplitPointKey& point_key, const CoincideSegment& coincide_segment) {
    if (coincide_segment.SplitterIndex != -1) {
        return IsInnerSplitter(container, index, point_key, container->Splitters.at(coincide_segment.SplitterIndex));
    }
    else {
        for (auto itr = coincide_segment.SplitterIndices[index].begin(); itr != coincide_segment.SplitterIndices[index].end(); ++itr) {
            if (IsInnerSplitter(container, index, point_key, container->Splitters.at(*itr))) {
                return true;
            }
        }
        return false;
    }
}

void WGTopoHelper2d::BuildSegmentsCoincideSegmentIndex(SplitterContainer* container, int index) {
    for (auto segment_itr = container->Segments[index].begin(); segment_itr != container->Segments[index].end(); ++segment_itr) {
        Segment& segment = *segment_itr;
        segment.CoincideSegmentIndex = -1;
        if (segment.Position == SegmentPosition::Coincide) {
            const SplitPoint& start_point = GetMaxSplitPoint(container, index, segment.Start);
            const SplitPoint& end_point = GetMinSplitPoint(container, index, segment.End);
            SplitPointKey point_key = start_point.Key;
            if (end_point.Key.EdgeIndex == start_point.Key.EdgeIndex && end_point.Key.PieceIndex == start_point.Key.PieceIndex) {
                point_key.T = (start_point.Key.T + end_point.Key.T) * 0.5;
            }
            else {
                const WGCurve2d* curve = start_point.Key.Wire->GetEdge(start_point.Key.EdgeIndex)->GetCurve();
                point_key.T = (start_point.Key.T + curve->GetPieceDomain(start_point.Key.PieceIndex).Max) * 0.5;
            }
            for (int i = 0; i < (int)container->CoincideSegments.size(); ++i) {
                if (IsInnerCoincideSegment(container, index, point_key, container->CoincideSegments.at(i))) {
                    segment.CoincideSegmentIndex = i;
                    break;
                }
            }
        }
    }
}

bool WGTopoHelper2d::CheckCoincideSegmentIndex(const Segment& segment) {
    if (segment.Position == SegmentPosition::Coincide && segment.CoincideSegmentIndex == -1) {
        return false;
    }
    return true;
}

WGTopoHelper2d::SegmentIndex::SegmentIndex(int index, int segment_index) :
    Index(index),
    SegIndex(segment_index) {
}

void WGTopoHelper2d::SearchLoop(SplitterContainer* container, int index, int segment_index, std::vector<Loop>& loops) {
    Segment& segment = container->Segments[index].at(segment_index);
    if (segment.SelectedFlag != 1 && segment.SelectedFlag != -1) {
        return;
    }
    Loop loop;
    loop.SegmentIndices.reserve(16);
    loop.SegmentIndices.push_back(SegmentIndex(index, segment_index));
    segment.SelectedFlag *= 2;
    if (segment.End.IntersectionIndex == segment.Start.IntersectionIndex) {
        loop.Closed = true;
        loops.push_back(std::move(loop));
        return;
    }
    int intersection_index = segment.SelectedFlag > 0 ? segment.End.IntersectionIndex : segment.Start.IntersectionIndex;
    while (true) {
        bool b = false;
        for (int i = index; i <= 1; ++i) {
            for (int j = i == index ? segment_index + 1 : 0; j < (int)container->Segments[i].size(); ++j) {
                Segment& segment = container->Segments[i].at(j);
                if (segment.SelectedFlag == 1) {
                    if (segment.Start.IntersectionIndex == intersection_index) {
                        loop.SegmentIndices.push_back(SegmentIndex(i, j));
                        segment.SelectedFlag = 2;
                        intersection_index = segment.End.IntersectionIndex;
                        b = true;
                        break;
                    }
                }
                else if (segment.SelectedFlag == -1) {
                    if (segment.End.IntersectionIndex == intersection_index) {
                        loop.SegmentIndices.push_back(SegmentIndex(i, j));
                        segment.SelectedFlag = -2;
                        intersection_index = segment.Start.IntersectionIndex;
                        b = true;
                        break;
                    }
                }
            }
            if (b) {
                break;
            }
        }
        if (!b) {
            for (int i = index; i <= 1; ++i) {
                for (int j = i == index ? segment_index + 1 : 0; j < (int)container->Segments[i].size(); ++j) {
                    Segment& segment = container->Segments[i].at(j);
                    if (segment.SelectedFlag == 0 && segment.Position == SegmentPosition::Coincide) {
                        const CoincideSegment& coincide_segment = container->CoincideSegments.at(segment.CoincideSegmentIndex);
                        if (coincide_segment.Flag == 0) {
                            if (segment.Start.IntersectionIndex == intersection_index) {
                                loop.SegmentIndices.push_back(SegmentIndex(i, j));
                                segment.SelectedFlag = 2;
                                intersection_index = segment.End.IntersectionIndex;
                                b = true;
                                break;
                            }
                            if (segment.End.IntersectionIndex == intersection_index) {
                                loop.SegmentIndices.push_back(SegmentIndex(i, j));
                                segment.SelectedFlag = -2;
                                intersection_index = segment.Start.IntersectionIndex;
                                b = true;
                                break;
                            }
                        }
                    }
                }
                if (b) {
                    break;
                }
            }
        }
        if (!b) {
            loop.Closed = false;
            loops.push_back(std::move(loop));
            return;
        }
        int i = (int)loop.SegmentIndices.size() - 1;
        while (i >= 0) {
            const SegmentIndex& segment_index = loop.SegmentIndices.at(i);
            const Segment& segment = container->Segments[segment_index.Index].at(segment_index.SegIndex);
            if (segment.SelectedFlag > 0) {
                if (intersection_index == segment.Start.IntersectionIndex) {
                    break;
                }
            }
            else {
                if (intersection_index == segment.End.IntersectionIndex) {
                    break;
                }
            }
            --i;
        }
        if (i == 0) {
            loop.Closed = true;
            loops.push_back(std::move(loop));
            return;
        }
        else if (i > 0) {
            Loop loop1;
            loop1.SegmentIndices.insert(loop1.SegmentIndices.begin(), loop.SegmentIndices.begin() + i, loop.SegmentIndices.end());
            loop1.Closed = true;
            loops.push_back(std::move(loop1));
            loop.SegmentIndices.resize(i);
        }
    }
}

void WGTopoHelper2d::SearchLoops(SplitterContainer* container, std::vector<Loop>& loops) {
    for (int i = 0; i <= 1; ++i) {
        for (int j = 0; j < (int)container->Segments[i].size(); ++j) {
            SearchLoop(container, i, j, loops);
        }
    }
}

bool WGTopoHelper2d::CheckLoop(const Loop& loop) {
    assert(loop.SegmentIndices.size() > 0);
    if (!loop.Closed) {
        return false;
    }
    return true;
}

void WGTopoHelper2d::RepairEdgeStartEnds(SplitterContainer* container) {
    //不同的需求会产生不同的策略
    //如需定制可和作者联系
    //Different strategies can be developed here based on different needs
    //If customization is required, please contact the author
    for (int index = 0; index <= 1; ++index) {
        for (auto segment_itr = container->Segments[index].begin(); segment_itr != container->Segments[index].end(); ++segment_itr) {
            Segment& segment = *segment_itr;
            if (segment.SelectedFlag == 2 || segment.SelectedFlag == -2) {
                segment.StartKey = GetMaxSplitPoint(container, index, segment.Start).Key;
                segment.EndKey = GetMinSplitPoint(container, index, segment.End).Key;
                segment.IsStartExtend = false;
                segment.IsEndExtend = false;
            }
        }
    }
}

void WGTopoHelper2d::AddEdges(WGWire2d* wire, const SplitPointKey& start_key, const SplitPointKey& end_key, bool reverse) {
    SplitPointKey current_key = start_key;
    int start_index = wire->GetEdgeCount();
    while (true) {
        if (current_key.EdgeIndex == end_key.EdgeIndex) {
            if (current_key.PieceIndex < end_key.PieceIndex) {
                WGEdge2d* edge = new WGEdge2d(current_key.Wire->GetEdge(current_key.EdgeIndex)->GetCurve()->CreateSubCurve(
                    current_key.PieceIndex, current_key.T, end_key.PieceIndex, end_key.T));
                wire->AddEdge(edge);
                break;
            }
            if (current_key.PieceIndex == end_key.PieceIndex) {
                if (current_key.T < end_key.T) {
                    WGEdge2d* edge = new WGEdge2d(current_key.Wire->GetEdge(current_key.EdgeIndex)->GetCurve()->CreateSubCurve(
                        current_key.PieceIndex, current_key.T, end_key.PieceIndex, end_key.T));
                    wire->AddEdge(edge);
                    break;
                }
                if (current_key.T == end_key.T) {
                    if (wire->GetEdgeCount() > start_index) {
                        break;
                    }
                }
            }
        }
        const WGCurve2d* src_curve = current_key.Wire->GetEdge(current_key.EdgeIndex)->GetCurve();
        int next_piece_index = src_curve->GetPieceCount() - 1;
        double next_t = src_curve->GetPieceDomain(next_piece_index).Max;
        WGEdge2d* edge = new WGEdge2d(src_curve->CreateSubCurve(current_key.PieceIndex, current_key.T, next_piece_index, next_t));
        wire->AddEdge(edge);
        ++current_key.EdgeIndex;
        if (current_key.EdgeIndex == current_key.Wire->GetEdgeCount()) {
            current_key.EdgeIndex = 0;
        }
        current_key.PieceIndex = 0;
        current_key.T = current_key.Wire->GetEdge(current_key.EdgeIndex)->GetCurve()->GetPieceDomain(current_key.PieceIndex).Min;
    }
    if (reverse) {
        wire->Reverse(start_index);
    }
}

void WGTopoHelper2d::AddPolygonLoops(SplitterContainer* container, std::vector<Loop>& loops, WGPolygon* result_polygon) {
    for (auto loop_itr = loops.begin(); loop_itr != loops.end(); ++loop_itr) {
        const Loop& loop = *loop_itr;
        WGWire2d* wire = new WGWire2d();
        for (auto segment_index_itr = loop.SegmentIndices.begin(); segment_index_itr != loop.SegmentIndices.end(); ++segment_index_itr) {
            const SegmentIndex& segment_index = *segment_index_itr;
            const Segment& segment = container->Segments[segment_index.Index].at(segment_index.SegIndex);
            if (segment.SelectedFlag > 0) {
                if (segment.IsStartExtend) {
                    wire->AddEdge(new WGEdge2d(new WGLine2d(segment.StartExtendPoint, segment.StartKey.CalculatePoint())));
                }
                AddEdges(wire, segment.StartKey, segment.EndKey, false);
                if (segment.IsEndExtend) {
                    wire->AddEdge(new WGEdge2d(new WGLine2d(segment.EndKey.CalculatePoint(), segment.EndExtendPoint)));
                }
            }
            else {
                if (segment.IsEndExtend) {
                    wire->AddEdge(new WGEdge2d(new WGLine2d(segment.EndExtendPoint, segment.EndKey.CalculatePoint())));
                }
                AddEdges(wire, segment.StartKey, segment.EndKey, true);
                if (segment.IsStartExtend) {
                    wire->AddEdge(new WGEdge2d(new WGLine2d(segment.StartKey.CalculatePoint(), segment.StartExtendPoint)));
                }
            }
        }
        wire->SetClosed(true);
        WGLoop2d* result_loop = new WGLoop2d();
        result_loop->SetWire(wire, WGLoop2d::CalculateSignArea(wire));
        result_polygon->AddLoop(result_loop);
    }
}

bool WGTopoHelper2d::AddPolygonSeparatedLoops(SplitterContainer* container, const WGPolygon* polygon0, const WGPolygon* polygon1,
    bool is_box_separated, SegmentSelector* selector, double distance_epsilon, 
    bool& is_singularity0, bool& is_singularity1, WGPolygon* result_polygon) {
    is_singularity0 = false;
    is_singularity1 = false;
    const WGPolygon* polygons[2] = { polygon0, polygon1 };
    for (int index = 0; index <= 1; ++index) {
        for (int i = 0; i < polygons[index]->GetLoopCount(); ++i) {
            const WGLoop2d* loop = polygons[index]->GetLoop(i);
            const WGWire2d* wire = loop->GetWire();
            if (is_box_separated) {
                int selected_flag = selector->SelectWire(index, false);
                if (selected_flag == 1) {
                    WGWire2d* result_wire = wire->Clone();
                    WGLoop2d* result_loop = new WGLoop2d();
                    result_loop->SetWire(result_wire, loop->GetSignArea());
                    result_polygon->AddLoop(result_loop);
                }
                else if (selected_flag == -1) {
                    WGWire2d* result_wire = wire->Clone();
                    result_wire->Reverse();
                    WGLoop2d* result_loop = new WGLoop2d();
                    result_loop->SetWire(result_wire, -loop->GetSignArea());
                    result_polygon->AddLoop(result_loop);
                }
            }
            else {
                bool b = true;
                if (container) {
                    for (auto split_point_itr = container->SplitPoints[index].begin(); split_point_itr != container->SplitPoints[index].end(); ++split_point_itr) {
                        if (split_point_itr->Key.Wire == wire) {
                            b = false;
                            break;
                        }
                    }
                }
                if (b) {
                    bool b2 = true;
                    SplitPointKey start_key = SplitPointKey(wire, 0, 0, wire->GetEdge(0)->GetCurve()->GetPieceDomain(0).Min);
                    const WGCurve2d* curve = wire->GetEdge(wire->GetEdgeCount() - 1)->GetCurve();
                    SplitPointKey end_key = SplitPointKey(wire, wire->GetEdgeCount() - 1, curve->GetPieceCount() - 1,
                        curve->GetPieceDomain(curve->GetPieceCount() - 1).Max);
                    SampleIterator sample_iterator(start_key, end_key, 4);
                    sample_iterator.First();
                    while (!sample_iterator.Eof()) {
                        WGPointPolygonPosition2d::Result position = WGPointPolygonPosition2d::Execute(
                            sample_iterator.GetCurrentPoint(), polygons[index ^ 1], distance_epsilon);
                        if (position == WGPointPolygonPosition2d::Result::Singularity) {
                            if (index == 0) {
                                is_singularity0 = true;
                            }
                            else {
                                is_singularity1 = true;
                            }
                            return false;
                        }
                        int selected_flag = 0;
                        if (position == WGPointPolygonPosition2d::Result::Inner) {
                            selected_flag = selector->SelectWire(index, true);
                        }
                        else if (position == WGPointPolygonPosition2d::Result::Outter) {
                            selected_flag = selector->SelectWire(index, false);
                        }
                        if (selected_flag == 1) {
                            WGWire2d* result_wire = wire->Clone();
                            WGLoop2d* result_loop = new WGLoop2d();
                            result_loop->SetWire(result_wire, loop->GetSignArea());
                            result_polygon->AddLoop(result_loop);
                            b2 = false;
                            break;
                        }
                        if (selected_flag == -1) {
                            WGWire2d* result_wire = wire->Clone();
                            result_wire->Reverse();
                            WGLoop2d* result_loop = new WGLoop2d();
                            result_loop->SetWire(result_wire, -loop->GetSignArea());
                            result_polygon->AddLoop(result_loop);
                            b2 = false;
                            break;
                        }
                        sample_iterator.Next();
                    }
                    if (b2) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

void PolygonPolygonSubtractionSegmentSelector::SelectSegments(WGTopoHelper2d::SplitterContainer* container) {
    for (auto segment_itr = container->Segments[0].begin(); segment_itr != container->Segments[0].end(); ++segment_itr) {
        WGTopoHelper2d::Segment& segment = *segment_itr;
        if (segment.Position == WGTopoHelper2d::SegmentPosition::Outter) {
            segment.SelectedFlag = 1;
        }
        else if (segment.Position == WGTopoHelper2d::SegmentPosition::Coincide) {
            const WGTopoHelper2d::CoincideSegment& coincide_segment = container->CoincideSegments.at(segment.CoincideSegmentIndex);
            if (coincide_segment.Flag != 0) {
                if (coincide_segment.SplitterIndex != -1) {
                    if (coincide_segment.Flag == 1) {
                        segment.SelectedFlag = 0;
                    }
                    else {
                        segment.SelectedFlag = 1;
                    }
                }
                else {
                    int n0 = (int)coincide_segment.SplitterIndices[0].size();
                    int n1 = (int)coincide_segment.SplitterIndices[1].size();
                    if ((n0 & 1) == 0) {
                        if ((n1 & 1) == 0) {
                            segment.SelectedFlag = 0;
                        }
                        else {
                            segment.SelectedFlag = 0;
                        }
                    }
                    else {
                        if ((n1 & 1) == 0) {
                            if (coincide_segment.Flag == 2) {
                                segment.SelectedFlag = 1;
                            }
                            else {
                                segment.SelectedFlag = 0;
                            }
                        }
                        else {
                            if (coincide_segment.Flag == 1) {
                                segment.SelectedFlag = 0;
                            }
                            else {
                                segment.SelectedFlag = 1;
                            }
                        }
                    }
                }
            }
            else {
                segment.SelectedFlag = 0;
            }
        }
        else {
            segment.SelectedFlag = 0;
        }
    }
    for (auto segment_itr = container->Segments[1].begin(); segment_itr != container->Segments[1].end(); ++segment_itr) {
        WGTopoHelper2d::Segment& segment = *segment_itr;
        if (segment.Position == WGTopoHelper2d::SegmentPosition::Inner) {
            segment.SelectedFlag = -1;
        }
        else if (segment.Position == WGTopoHelper2d::SegmentPosition::Coincide) {
            const WGTopoHelper2d::CoincideSegment& coincide_segment = container->CoincideSegments.at(segment.CoincideSegmentIndex);
            if (coincide_segment.Flag != 0) {
                if (coincide_segment.SplitterIndex != -1) {
                    segment.SelectedFlag = 0;
                }
                else {
                    int n0 = (int)coincide_segment.SplitterIndices[0].size();
                    int n1 = (int)coincide_segment.SplitterIndices[1].size();
                    if ((n0 & 1) == 0) {
                        if ((n1 & 1) == 0) {
                            segment.SelectedFlag = 0;
                        }
                        else {
                            if (coincide_segment.Flag == 1) {
                                segment.SelectedFlag = -1;
                            }
                            else {
                                segment.SelectedFlag = 0;
                            }
                        }
                    }
                    else {
                        if ((n1 & 1) == 0) {
                            segment.SelectedFlag = 0;
                        }
                        else {
                            segment.SelectedFlag = 0;
                        }
                    }
                }
            }
            else {
                segment.SelectedFlag = 0;
            }
        }
        else {
            segment.SelectedFlag = 0;
        }
    }
}

int PolygonPolygonSubtractionSegmentSelector::SelectWire(int index, bool is_inner) {
    if (index == 0) {
        if (is_inner) {
            return 0;
        }
        return 1;
    }
    else {
        if (is_inner) {
            return -1;
        }
        return 0;
    }
}

void PolygonPolygonIntersectionSegmentSelector::SelectSegments(WGTopoHelper2d::SplitterContainer* container) {
    for (auto segment_itr = container->Segments[0].begin(); segment_itr != container->Segments[0].end(); ++segment_itr) {
        WGTopoHelper2d::Segment& segment = *segment_itr;
        if (segment.Position == WGTopoHelper2d::SegmentPosition::Inner) {
            segment.SelectedFlag = 1;
        }
        else if (segment.Position == WGTopoHelper2d::SegmentPosition::Coincide) {
            const WGTopoHelper2d::CoincideSegment& coincide_segment = container->CoincideSegments.at(segment.CoincideSegmentIndex);
            if (coincide_segment.Flag != 0) {
                if (coincide_segment.SplitterIndex != -1) {
                    if (coincide_segment.Flag == 1) {
                        segment.SelectedFlag = 1;
                    }
                    else {
                        segment.SelectedFlag = 0;
                    }
                }
                else {
                    int n0 = (int)coincide_segment.SplitterIndices[0].size();
                    int n1 = (int)coincide_segment.SplitterIndices[1].size();
                    if ((n0 & 1) == 0) {
                        if ((n1 & 1) == 0) {
                            segment.SelectedFlag = 0;
                        }
                        else {
                            segment.SelectedFlag = 0;
                        }
                    }
                    else {
                        if ((n1 & 1) == 0) {
                            if (coincide_segment.Flag == 1) {
                                segment.SelectedFlag = 1;
                            }
                            else {
                                segment.SelectedFlag = 0;
                            }
                        }
                        else {
                            if (coincide_segment.Flag == 1) {
                                segment.SelectedFlag = 1;
                            }
                            else {
                                segment.SelectedFlag = 0;
                            }
                        }
                    }
                }
            }
            else {
                segment.SelectedFlag = 0;
            }
        }
        else {
            segment.SelectedFlag = 0;
        }
    }
    for (auto segment_itr = container->Segments[1].begin(); segment_itr != container->Segments[1].end(); ++segment_itr) {
        WGTopoHelper2d::Segment& segment = *segment_itr;
        if (segment.Position == WGTopoHelper2d::SegmentPosition::Inner) {
            segment.SelectedFlag = 1;
        }
        else if (segment.Position == WGTopoHelper2d::SegmentPosition::Coincide) {
            const WGTopoHelper2d::CoincideSegment& coincide_segment = container->CoincideSegments.at(segment.CoincideSegmentIndex);
            if (coincide_segment.Flag != 0) {
                if (coincide_segment.SplitterIndex != -1) {
                    segment.SelectedFlag = 0;
                }
                else {
                    int n0 = (int)coincide_segment.SplitterIndices[0].size();
                    int n1 = (int)coincide_segment.SplitterIndices[1].size();
                    if ((n0 & 1) == 0) {
                        if ((n1 & 1) == 0) {
                            segment.SelectedFlag = 0;
                        }
                        else {
                            if (coincide_segment.Flag == 1) {
                                segment.SelectedFlag = 1;
                            }
                            else {
                                segment.SelectedFlag = 0;
                            }
                        }
                    }
                    else {
                        if ((n1 & 1) == 0) {
                            segment.SelectedFlag = 0;
                        }
                        else {
                            segment.SelectedFlag = 0;
                        }
                    }
                }
            }
            else {
                segment.SelectedFlag = 0;
            }
        }
        else {
            segment.SelectedFlag = 0;
        }
    }
}

int PolygonPolygonIntersectionSegmentSelector::SelectWire(int index, bool is_inner) {
    if (index == 0) {
        if (is_inner) {
            return 1;
        }
        return 0;
    }
    else {
        if (is_inner) {
            return 1;
        }
        return 0;
    }
}

void PolygonPolygonUnionSegmentSelector::SelectSegments(WGTopoHelper2d::SplitterContainer* container) {
    for (auto segment_itr = container->Segments[0].begin(); segment_itr != container->Segments[0].end(); ++segment_itr) {
        WGTopoHelper2d::Segment& segment = *segment_itr;
        if (segment.Position == WGTopoHelper2d::SegmentPosition::Outter) {
            segment.SelectedFlag = 1;
        }
        else if (segment.Position == WGTopoHelper2d::SegmentPosition::Coincide) {
            const WGTopoHelper2d::CoincideSegment& coincide_segment = container->CoincideSegments.at(segment.CoincideSegmentIndex);
            if (coincide_segment.Flag != 0) {
                if (coincide_segment.SplitterIndex != -1) {
                    if (coincide_segment.Flag == 1) {
                        segment.SelectedFlag = 1;
                    }
                    else {
                        segment.SelectedFlag = 0;
                    }
                }
                else {
                    int n0 = (int)coincide_segment.SplitterIndices[0].size();
                    int n1 = (int)coincide_segment.SplitterIndices[1].size();
                    if ((n0 & 1) == 0) {
                        if ((n1 & 1) == 0) {
                            segment.SelectedFlag = 0;
                        }
                        else {
                            segment.SelectedFlag = 0;
                        }
                    }
                    else {
                        if ((n1 & 1) == 0) {
                            if (coincide_segment.Flag == 2) {
                                segment.SelectedFlag = 1;
                            }
                            else {
                                segment.SelectedFlag = 0;
                            }
                        }
                        else {
                            if (coincide_segment.Flag == 1) {
                                segment.SelectedFlag = 1;
                            }
                            else {
                                segment.SelectedFlag = 0;
                            }
                        }
                    }
                }
            }
            else {
                segment.SelectedFlag = 0;
            }
        }
        else {
            segment.SelectedFlag = 0;
        }
    }
    for (auto segment_itr = container->Segments[1].begin(); segment_itr != container->Segments[1].end(); ++segment_itr) {
        WGTopoHelper2d::Segment& segment = *segment_itr;
        if (segment.Position == WGTopoHelper2d::SegmentPosition::Outter) {
            segment.SelectedFlag = 1;
        }
        else if (segment.Position == WGTopoHelper2d::SegmentPosition::Coincide) {
            const WGTopoHelper2d::CoincideSegment& coincide_segment = container->CoincideSegments.at(segment.CoincideSegmentIndex);
            if (coincide_segment.Flag != 0) {
                if (coincide_segment.SplitterIndex != -1) {
                    segment.SelectedFlag = 0;
                }
                else {
                    int n0 = (int)coincide_segment.SplitterIndices[0].size();
                    int n1 = (int)coincide_segment.SplitterIndices[1].size();
                    if ((n0 & 1) == 0) {
                        if ((n1 & 1) == 0) {
                            segment.SelectedFlag = 0;
                        }
                        else {
                            if (coincide_segment.Flag == 1) {
                                segment.SelectedFlag = 0;
                            }
                            else {
                                segment.SelectedFlag = 1;
                            }
                        }
                    }
                    else {
                        if ((n1 & 1) == 0) {
                            segment.SelectedFlag = 0;
                        }
                        else {
                            segment.SelectedFlag = 0;
                        }
                    }
                }
            }
            else {
                segment.SelectedFlag = 0;
            }
        }
        else {
            segment.SelectedFlag = 0;
        }
    }
}

int PolygonPolygonUnionSegmentSelector::SelectWire(int index, bool is_inner) {
    if (index == 0) {
        if (is_inner) {
            return 0;
        }
        return 1;
    }
    else {
        if (is_inner) {
            return 0;
        }
        return 1;
    }
}