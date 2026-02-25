#ifndef _WG_NURBS_CURVE_2D_
#define _WG_NURBS_CURVE_2D_

#include "WGCurve2d.h"

class WGNurbsCurve2d : public WGCurve2d {
public:
    WG_OBJECT_TYPE_DEF(WGNurbsCurve2d, WGCurve2d);
public:
    WGNurbsCurve2d();
    WGNurbsCurve2d(int degree, int control_point_count, bool weights);
    WGNurbsCurve2d(int degree, int control_point_count, const double* knots, const WGVector2d* control_points, const double* weights);
    virtual ~WGNurbsCurve2d();
    virtual WGVector2d GetStartPoint() const;
    virtual WGVector2d GetEndPoint() const;
    WGNurbsCurve2d* InsertKnot(int piece_index, double t, int insert_count) const;
    int GetDegree() const;
    int GetControlPointCount() const;
    int GetKnotCount() const;
    const double* GetKnots() const;
    const WGVector2d* GetControlPoints() const;
    const double* GetWeights() const;
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
    static WGNurbsCurve2d* BuildByControlPoints(int control_point_count, const WGVector2d* control_points);
public:
    class BezierIterator {
    public:
        BezierIterator(const WGNurbsCurve2d* curve, int piece_index, int piece_count);
        virtual ~BezierIterator();
        void First();
        bool Eof();
        void Next();
        int GetCurrentPieceIndex() const;
        WSInterval GetCurrentDomain() const;
        const WGVector2d* GetCurrentColtrolPoints() const;
        const double* GetCurrentWeights() const;
    private:
        void SkipZero();
        void InsertNext();
    private:
        const WGNurbsCurve2d* m_curve;
        int m_piece_index;
        int m_piece_count;
        bool m_is_bezier;
        int m_offset;
        int m_control_point_count;
        double* m_knot_buffer;
        WGVector2d* m_control_point_buffer;
        double* m_weight_buffer;
        mutable int m_current_piece_index;
        double* m_current_knots;
        WGVector2d* m_current_control_points;
        double* m_current_weights;
    };
protected:
    int m_degree;
    int m_control_point_count;
    double* m_knots;
    WGVector2d* m_control_points;
    double* m_weights;
};

#endif
