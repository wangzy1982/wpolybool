#ifndef _WG_OBJECT_2D_
#define _WG_OBJECT_2D_

#include "WGObject.h"
#include "WGVector2d.h"

class WGObject2d : public WGObject {
public:
    WG_OBJECT_TYPE_DEF(WGObject2d, WGObject);
public:
    virtual void Move(const WGVector2d& vt) = 0;
};

#endif