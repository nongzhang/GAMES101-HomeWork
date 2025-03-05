#pragma once

#include "Vector.hpp"

/// <summary>
/// 点光源
/// </summary>
class Light
{
public:
    Light(const Vector3f& p, const Vector3f& i)
        : position(p)
        , intensity(i)
    {}
    virtual ~Light() = default;
    Vector3f position;   
    Vector3f intensity;   //光源强度也是三维向量
};
