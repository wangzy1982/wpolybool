#ifndef _WG_ARC_2D_
#define _WG_ARC_2D_

#include "WGCurve2d.h"

class WGArc2d : public WGCurve2d {
public:
    WG_OBJECT_TYPE_DEF(WGArc2d, WGCurve2d);
public:
    WGArc2d();
    WGArc2d(const WGVector2d& center, double radius, double start_angle, double delta_angle);
    virtual WGVector2d GetStartPoint() const;
    virtual WGVector2d GetEndPoint() const;
    const WGVector2d& GetCenter() const;
    double GetRadius() const;
    double GetStartAngle() const;
    double GetDeltaAngle() const;
    void SetCenter(const WGVector2d& center);
    void SetRadius(double radius);
    void SetStartAngle(double start_angle);
    void SetDeltaAngle(double delta_angle);
public:
    double GetT(double angle) const;
    double GetAngle(double t) const;
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
    static double CalculateT(double start_angle, double delta_angle, double angle);
    static WGVector2d CalculatePoint(const WGVector2d& center, double radius, double angle);
public:
    static bool BuildByPoints(const WGVector2d& point0, const WGVector2d& point1, const WGVector2d& point2,
        WGVector2d& center, double& radius, double& start_angle, double& delta_angle);
    static WGArc2d* BuildByPoints(const WGVector2d& point0, const WGVector2d& point1, const WGVector2d& point2);
protected:
    WGVector2d m_center;
    double m_radius;
    double m_start_angle;
    double m_delta_angle;
};

#endif
