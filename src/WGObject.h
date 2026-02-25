#ifndef _WG_OBJECT_
#define _WG_OBJECT_

class WGObject {
public:
    class Type {
    public:
        Type();
        Type(Type* parent);
        virtual ~Type();
        bool IsA(Type* type);
    protected:
        Type* m_parent;
    public:
        static Type* Instance();
    private:
        static Type m_instance;
    };
    virtual Type* GetType() const;
public:
    virtual ~WGObject();
};

#define WG_OBJECT_TYPE_DEF(class_name, parent_class_name) \
    class Type : public WGObject::Type { \
    public: \
        Type(WGObject::Type* parent); \
        static Type* Instance(); \
    private: \
        static Type m_instance; \
    }; \
    virtual WGObject::Type* GetType() const;

#define WG_OBJECT_TYPE_IMP(class_name, parent_class_name) \
    class_name::Type::Type(WGObject::Type* parent) : WGObject::Type(parent) { \
    } \
    class_name::Type* class_name::Type::Instance() { \
        return &m_instance; \
    } \
    class_name::Type class_name::Type::m_instance = class_name::Type(parent_class_name::Type::Instance()); \
    WGObject::Type* class_name::GetType() const { \
        return Type::Instance(); \
    }

#endif