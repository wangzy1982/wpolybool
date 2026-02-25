#ifndef _WG_POLYGON_
#define _WG_POLYGON_

#include "WGLoop2d.h"

class WGPolygon : public WGObject2d {
public:
    WG_OBJECT_TYPE_DEF(WGPolygon, WGObject2d);
public:
    WGPolygon();
    virtual ~WGPolygon();
    void Reserve(int loop_capacity);
    void AddLoop(WGLoop2d* loop);
    void RemoveLoop(int index);
    void Clear();
    int GetLoopCount() const;
    const WGLoop2d* GetLoop(int index) const;
    const WGBox2d& GetBox() const;
    WGPolygon* Clone() const;
    double GetJointEpsilon() const;
    double CalculateArea() const;
public:
    bool RepairSingularity(double distance_epsilon);
public:
    virtual void Move(const WGVector2d& vt);
public:
    static WGPolygon* BuildPolygon(WGCurve2d** curves, int curve_count);
private:
    std::vector<WGLoop2d*> m_loops;
    mutable WGBox2d* m_box;
    mutable double* m_joint_epsilon;
};

#endif
