//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_BOUNDS3_H
#define RAYTRACING_BOUNDS3_H
#include "Ray.hpp"
#include "Vector.hpp"
#include <limits>
#include <array>


/// <summary>
/// 定义了包围盒Bounds3类，即一个典型的AABB。由于轴向对齐，Bounds3的构造函数以长方体对角线上的一对点为入参。
/// 无参构造函数则将两点分别设置为(1.7976931348623158e+308, 1.7976931348623158e+308, 1.7976931348623158e+308)和
/// (-1.7976931348623158e+308,-1.7976931348623158e+308,-1.7976931348623158e+308)。
/// 这个超大的数字是标准库中双精度浮点型double的最大值。单入参构造函数则表示一个点，长宽高都是0。
/// maxExtent 函数用于计算包围盒哪个轴向的值最大，x轴的值最长则返回0，y轴的值最长则返回1，z轴的值最长则返回2。
/// Overlaps 和Inside 两个函数分别判断两个Bounds3对象是否存在重叠和包含关系。
/// 两个Union函数会返回合并后的Bounds3对象，入参可以是两个Bounds3对象，也可以是一个Bounds3对象和一个点(也可以看成是单点构造的Bounds3对象)

/// </summary>
class Bounds3
{
  public:
    Vector3f pMin, pMax; // two points to specify the bounding box
    //这样初始化的目的是确保后续更新边界框时，任何新的点都会调整边界框的大小。最初时，pMax 和 pMin 会分别表示极端的负无穷大和正无穷大。
    Bounds3()
    {
        double minNum = std::numeric_limits<double>::lowest();
        double maxNum = std::numeric_limits<double>::max();
        pMax = Vector3f(minNum, minNum, minNum);
        pMin = Vector3f(maxNum, maxNum, maxNum);
    }
    /// <summary>
    /// 使用一个点 p 初始化边界框，将 pMin 和 pMax 都设置为这个点。这表示边界框在初始化时只有一个点，pMin 和 pMax 都是相同的
    /// </summary>
    /// <param name="p"></param>
    Bounds3(const Vector3f p) : pMin(p), pMax(p) {}
    Bounds3(const Vector3f p1, const Vector3f p2)
    {
        pMin = Vector3f(fmin(p1.x, p2.x), fmin(p1.y, p2.y), fmin(p1.z, p2.z));
        pMax = Vector3f(fmax(p1.x, p2.x), fmax(p1.y, p2.y), fmax(p1.z, p2.z));
    }

    Vector3f Diagonal() const { return pMax - pMin; }
    int maxExtent() const
    {
        Vector3f d = Diagonal();
        if (d.x > d.y && d.x > d.z)
            return 0;
        else if (d.y > d.z)
            return 1;
        else
            return 2;
    }

    double SurfaceArea() const
    {
        Vector3f d = Diagonal();
        return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
    }

    Vector3f Centroid() { return 0.5 * pMin + 0.5 * pMax; }
    Bounds3 Intersect(const Bounds3& b)
    {
        return Bounds3(Vector3f(fmax(pMin.x, b.pMin.x), fmax(pMin.y, b.pMin.y),
                                fmax(pMin.z, b.pMin.z)),
                       Vector3f(fmin(pMax.x, b.pMax.x), fmin(pMax.y, b.pMax.y),
                                fmin(pMax.z, b.pMax.z)));
    }

    Vector3f Offset(const Vector3f& p) const
    {
        Vector3f o = p - pMin;
        if (pMax.x > pMin.x)
            o.x /= pMax.x - pMin.x;
        if (pMax.y > pMin.y)
            o.y /= pMax.y - pMin.y;
        if (pMax.z > pMin.z)
            o.z /= pMax.z - pMin.z;
        return o;
    }

    bool Overlaps(const Bounds3& b1, const Bounds3& b2)
    {
        bool x = (b1.pMax.x >= b2.pMin.x) && (b1.pMin.x <= b2.pMax.x);
        bool y = (b1.pMax.y >= b2.pMin.y) && (b1.pMin.y <= b2.pMax.y);
        bool z = (b1.pMax.z >= b2.pMin.z) && (b1.pMin.z <= b2.pMax.z);
        return (x && y && z);
    }

    bool Inside(const Vector3f& p, const Bounds3& b)
    {
        return (p.x >= b.pMin.x && p.x <= b.pMax.x && p.y >= b.pMin.y &&
                p.y <= b.pMax.y && p.z >= b.pMin.z && p.z <= b.pMax.z);
    }
    inline const Vector3f& operator[](int i) const
    {
        return (i == 0) ? pMin : pMax;
    }

    inline bool IntersectP(const Ray& ray, const Vector3f& invDir,
                           const std::array<int, 3>& dirisNeg) const;
};



inline bool Bounds3::IntersectP(const Ray& ray, const Vector3f& invDir,
                                const std::array<int, 3>& dirIsNeg) const
{
    // invDir: ray direction(x,y,z), invDir=(1.0/x,1.0/y,1.0/z), use this because Multiply is faster that Division
    // dirIsNeg: ray direction(x,y,z), dirIsNeg=[int(x>0),int(y>0),int(z>0)], use this to simplify your logic
    // TODO test if ray bound intersects
    Vector3f intersectMin = (pMin - ray.origin) * invDir;
    Vector3f intersectMax = (pMax - ray.origin) * invDir;
    float iMin[3] = { intersectMin.x, intersectMin.y, intersectMin.z };
    float iMax[3] = { intersectMax.x, intersectMax.y, intersectMax.z };

    float tEnterArray[3], tExitArray[3];
    for (int i = 0; i < 3; i++)
    {
        tEnterArray[i] = dirIsNeg[i] ? iMin[i] : iMax[i];
        tExitArray[i] = dirIsNeg[i] ? iMax[i] : iMin[i];
    }

    float tEnter = std::max(std::max(tEnterArray[0], tEnterArray[1]), tEnterArray[2]);
    float tExit = std::min(std::min(tExitArray[0], tExitArray[1]), tExitArray[2]);
    if (tEnter < tExit && tExit >= 0)
    {
        return true;
    }
    return false;
}

inline Bounds3 Union(const Bounds3& b1, const Bounds3& b2)
{
    Bounds3 ret;
    ret.pMin = Vector3f::Min(b1.pMin, b2.pMin);
    ret.pMax = Vector3f::Max(b1.pMax, b2.pMax);
    return ret;
}

inline Bounds3 Union(const Bounds3& b, const Vector3f& p)
{
    Bounds3 ret;
    ret.pMin = Vector3f::Min(b.pMin, p);
    ret.pMax = Vector3f::Max(b.pMax, p);
    return ret;
}

#endif // RAYTRACING_BOUNDS3_H
