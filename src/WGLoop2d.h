#ifndef _WG_LOOP_2D_
#define _WG_LOOP_2D_

#include "WGWire2d.h"

class WGLoop2d : public WGObject2d {
public:
    WG_OBJECT_TYPE_DEF(WGLoop2d, WGObject2d);
public:
    WGLoop2d();
    virtual ~WGLoop2d();
    const WGWire2d* GetWire() const;
    void SetWire(WGWire2d* wire, double sign_area);
    bool IsClockwise() const;
    double GetArea() const;
    double GetSignArea() const;
    WGBox2d CalculateBox() const;
    WGLoop2d* Clone() const;
public:
    bool RepairSingularity(double distance_epsilon);
public:
    virtual void Move(const WGVector2d& vt);
public:
    static double CalculateSignArea(const WGWire2d* wire);
private:
    WGWire2d* m_wire;
    double m_sign_area;
};

#endif
