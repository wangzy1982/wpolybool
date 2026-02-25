#ifndef _WG_LINE_2D_
#define _WG_LINE_2D_

#include "WGCurve2d.h"

class WGLine2d : public WGCurve2d {
public:
    WG_OBJECT_TYPE_DEF(WGLine2d, WGCurve2d);
public:
    WGLine2d();
    WGLine2d(const WGVector2d& start_point, const WGVector2d& end_point);
    virtual WGVector2d GetStartPoint() const;
    virtual WGVector2d GetEndPoint() const;
    void SetStartPoint(const WGVector2d& point);
    void SetEndPoint(const WGVector2d& point);
public:
    virtual int GetPieceCount() const;
    virtual WSInterval GetPieceDomain(int index) const;
    virtual WGVector2d CalculateD0(int piece_index, double t) const;
    virtual WGVector2d CalculateD1(int piece_index, double t) const;
    virtual double CalculateVariableEpsilon(int piece_index, double distance_epsilon) const;
    virtual WGBox2d CalculateBox() const;
    virtual double CalculateAreaBelow(int piece_index, double start_t, double end_t) const;
    virtual WGCurve2d* Clone() const;
    virtual WGCurve2d* CreateSubCurve(int start_piece_index, double start_t, int end_piece_index, double end_t) const;
    virtual void Reverse();
    virtual void Linearize(std::vector<WGVector2d>& vertices, bool is_contiguous_with_prev,
        double limit_distance_epsilon, double limit_angle_epsilon, double distance_epsilon, double angle_epsilon) const;
public:
    virtual void Move(const WGVector2d& vt);
public:
    static double CalculateAreaBelow(const WGVector2d& start_point, const WGVector2d& end_point);
protected:
    WGVector2d m_start_point;
    WGVector2d m_end_point;
};

#endif
