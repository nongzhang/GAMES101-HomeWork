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

//���һԪ���η���
inline bool solveQuadratic(const float& a, const float& b, const float& c, float& x0, float& x1)
{
    float discr = b * b - 4 * a * c;
    if (discr < 0)
        return false;
    else if (discr == 0)
        x0 = x1 = -0.5 * b / a;
    else
    {
        float q = (b > 0) ? -0.5 * (b + sqrt(discr)) : -0.5 * (b - sqrt(discr));  //�����ж�b>0�������ж� discr>0 ��ֹ�������
        x0 = q / a;
        x1 = c / q;                                                                                    //x0��x1�Ƿ��̵������⣬��ôx0*x1 = c/a, ������ x0 = (-b+sqrt(discr))/2a ;  x1 = (-b-sqrt(discr))/2a ���������
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
/// �˺����������Ƿ���һ��������������ø������� 0 �� 1 ֮�䣨���� 0���������� 1����
/// ��ͨ�� std::random_device ����ȡһ�����ӣ���ͨ�� std::mt19937 ��������α�������
/// ����ͨ�� std::uniform_real_distribution<float> ��֤���ɵ��������ָ�������ھ��ȷֲ���
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
/// �����������һЩ�ַ���ʾ��Ⱦ����
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
