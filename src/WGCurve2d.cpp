#include "WGCurve2d.h"

WG_OBJECT_TYPE_IMP(WGCurve2d, WGObject2d);

double WGCurve2d::CalculateAreaBelow() {
    double d = 0;
    for (int i = 0, n = GetPieceCount(); i < n; ++i) {
        WSInterval domain = GetPieceDomain(i);
        d += CalculateAreaBelow(i, domain.Min, domain.Max);
    }
    return d;
}