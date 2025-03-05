#pragma once

#include <cmath>
#include <iostream>
#include <random>

#define M_PI 3.14159265358979323846

constexpr float kInfinity = std::numeric_limits<float>::max();

inline float clamp(const float& lo, const float& hi, const float& v)
{
    return std::max(lo, std::min(hi, v));
}

//求解一元二次方程
inline bool solveQuadratic(const float& a, const float& b, const float& c, float& x0, float& x1)
{
    float discr = b * b - 4 * a * c;
    if (discr < 0)
        return false;
    else if (discr == 0)
        x0 = x1 = -0.5 * b / a;
    else
    {
        float q = (b > 0) ? -0.5 * (b + sqrt(discr)) : -0.5 * (b - sqrt(discr));  //这里判断b>0而不是判断 discr>0 防止舍入误差
        x0 = q / a;
        x1 = c / q;                                                                                    //x0和x1是方程的两个解，那么x0*x1 = c/a, 这是由 x0 = (-b+sqrt(discr))/2a ;  x1 = (-b-sqrt(discr))/2a 计算得来的
    }
    if (x0 > x1)
        std::swap(x0, x1);
    return true;
}

enum MaterialType
{
    DIFFUSE_AND_GLOSSY,
    REFLECTION_AND_REFRACTION,
    REFLECTION
};

/// <summary>
/// 此函数的作用是返回一个随机浮点数，该浮点数在 0 和 1 之间（包括 0，但不包括 1）。
/// 它通过 std::random_device 来获取一个种子，再通过 std::mt19937 引擎生成伪随机数，
/// 最终通过 std::uniform_real_distribution<float> 保证生成的随机数在指定区间内均匀分布。
/// </summary>
/// <returns></returns>
inline float get_random_float()
{
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<float> dist(0.f, 1.f); 

    return dist(rng);
}

/// <summary>
/// 在命令行输出一些字符表示渲染进度
/// </summary>
/// <param name="progress"></param>
inline void UpdateProgress(float progress)
{
    int barWidth = 70;

    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i)
    {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}
