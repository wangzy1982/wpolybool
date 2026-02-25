#ifndef _WG_CURVE_2D_
#define _WG_CURVE_2D_

#include "WGUtil.h"
#include "WGObject2d.h"
#include "WGBox2d.h"
#include "wsolver.h"
#include <vector>

class WGCurve2d : public WGObject2d {
public:
    WG_OBJECT_TYPE_DEF(WGCurve2d, WGObject2d);
public:
    virtual int GetPieceCount() const = 0;
    virtual WSInterval GetPieceDomain(int index) const = 0;
    virtual WGVector2d GetStartPoint() const = 0;
    virtual WGVector2d GetEndPoint() const = 0;
    virtual WGVector2d CalculateD0(int piece_index, double t) const = 0;
    virtual WGVector2d CalculateD1(int piece_index, double t) const = 0;
    virtual double CalculateVariableEpsilon(int piece_index, double distance_epsilon) const = 0;
    virtual WGBox2d CalculateBox() const = 0;
    virtual double CalculateAreaBelow(int piece_index, double start_t, double end_t) const = 0;
    virtual WGCurve2d* Clone() const = 0;
    virtual WGCurve2d* CreateSubCurve(int start_piece_index, double start_t, int end_piece_index, double end_t) const = 0;
    virtual void Reverse() = 0;
    //(distance_epsilon || equal andgle_epsilon) && limit_distance_epsilon && limit_angle_epsilon
    virtual void Linearize(std::vector<WGVector2d>& vertices, bool is_contiguous_with_prev,
        double limit_distance_epsilon, double limit_angle_epsilon, double distance_epsilon, double angle_epsilon) const = 0;
public:
    double CalculateAreaBelow();
};

#endif
