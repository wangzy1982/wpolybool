#ifndef _WG_WIRE_2D_
#define _WG_WIRE_2D_

#include "WGEdge2d.h"
#include <vector>

class WGWire2d : public WGObject2d {
public:
    WG_OBJECT_TYPE_DEF(WGWire2d, WGObject2d);
public:
    WGWire2d();
    virtual ~WGWire2d();
    void Reserve(int edge_capacity);
    int AddEdge(WGEdge2d* edge);
    void RemoveEdge(int index);
    void Clear();
    void SetClosed(bool closed);
    bool GetClosed() const;
    int GetEdgeCount() const;
    const WGEdge2d* GetEdge(int index) const;
    bool CheckSelfIntersection() const;
    void Reverse(int start_index);
    void Reverse();
    WGBox2d CalculateBox() const;
    WGWire2d* Clone() const;
public:
    bool RepairSingularity(double distance_epsilon);
public:
    virtual void Move(const WGVector2d& vt);
private:
    bool m_closed;
    std::vector<WGEdge2d*> m_edges;
};

#endif