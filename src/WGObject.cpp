#include "WGObject.h"

WGObject::Type::Type() :
    m_parent(nullptr) {
}

WGObject::Type::Type(Type* parent) :
    m_parent(parent) {
}

WGObject::Type::~Type() {
}

bool WGObject::Type::IsA(Type* type) {
    if (this == type) {
        return true;
    }
    if (m_parent) {
        return m_parent->IsA(type);
    }
    return false;
}

WGObject::Type* WGObject::Type::Instance() {
    return &m_instance;
}

WGObject::Type WGObject::Type::m_instance;

WGObject::Type* WGObject::GetType() const {
    return Type::Instance();
}

WGObject::~WGObject() {
}