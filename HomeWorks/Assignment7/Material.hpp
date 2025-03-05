//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_MATERIAL_H
#define RAYTRACING_MATERIAL_H
//#include <eigen3/Eigen/Eigen>
#include "Vector.hpp"
#define PI 3.14159265358979323846
enum MaterialType { DIFFUSE, MICROFACET};

class Material{
private:
    /// <summary>
    /// 计算向量模长
    /// </summary>
    /// <param name="v"></param>
    /// <returns></returns>
    float length(const Vector3f &v)
    {
        return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    }

    /// <summary>
    /// 计算半程向量与表面法线之间的夹角（弧度）
    /// </summary>
    /// <param name="omega_h">半程向量</param>
    /// <param name="N">法线</param>
    /// <returns></returns>
    float angleBetween(const Vector3f& omega_h, const Vector3f& N)
    {
        float dot = dotProduct(omega_h, N);
        float len_omega_h = length(omega_h);
        float len_N = length(N);

        float cos_theta = dot / (len_omega_h * len_N);
        return acos(cos_theta);
    }

    // Compute reflection direction
    Vector3f reflect(const Vector3f &I, const Vector3f &N) const
    {
        return I - 2 * dotProduct(I, N) * N;
    }

    // Compute refraction direction using Snell's law
    //
    // We need to handle with care the two possible situations:
    //
    //    - When the ray is inside the object
    //
    //    - When the ray is outside.
    //
    // If the ray is outside, you need to make cosi positive cosi = -N.I
    //
    // If the ray is inside, you need to invert the refractive indices and negate the normal N
    Vector3f refract(const Vector3f &I, const Vector3f &N, const float &ior) const
    {
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        Vector3f n = N;
        if (cosi < 0) { cosi = -cosi; } else { std::swap(etai, etat); n= -N; }
        float eta = etai / etat;
        float k = 1 - eta * eta * (1 - cosi * cosi);
        return k < 0 ? 0 : eta * I + (eta * cosi - sqrtf(k)) * n;
    }

    // Compute Fresnel equation
    //
    // \param I is the incident view direction
    //
    // \param N is the normal at the intersection point
    //
    // \param ior is the material refractive index
    //
    // \param[out] kr is the amount of light reflected
    void fresnel(const Vector3f &I, const Vector3f &N, const float &ior, float &kr) const
    {
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        if (cosi > 0) {  std::swap(etai, etat); }
        // Compute sini using Snell's law
        float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
        // Total internal reflection
        if (sint >= 1) {
            kr = 1;
        }
        else {
            float cost = sqrtf(std::max(0.f, 1 - sint * sint));
            cosi = fabsf(cosi);
            float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
            float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
            kr = (Rs * Rs + Rp * Rp) / 2;
        }
        // As a consequence of the conservation of energy, transmittance is given by:
        // kt = 1 - kr;
    }

    Vector3f toWorld(const Vector3f &a, const Vector3f &N){
        Vector3f B, C;
        if (std::fabs(N.x) > std::fabs(N.y)){
            float invLen = 1.0f / std::sqrt(N.x * N.x + N.z * N.z);
            C = Vector3f(N.z * invLen, 0.0f, -N.x *invLen);
        }
        else {
            float invLen = 1.0f / std::sqrt(N.y * N.y + N.z * N.z);
            C = Vector3f(0.0f, N.z * invLen, -N.y *invLen);
        }
        B = crossProduct(C, N);
        return a.x * B + a.y * C + a.z * N;
    }

    //法线分布，视频中的D(h)项， h是半程向量，这里使用blinn-phong分布，公式为D(h) = (α+2)/(2π) * (n·h)^α
    void distributionH(const Vector3f& normal, const Vector3f& half, const float& roughness, float& D)
    {
        D = (roughness + 2.0f) / (M_PI * 2.0f) * pow(std::max(dotProduct(normal, half), EPSILON), roughness);
    }

    /// <summary>
    /// 计算法线分布函数 D (Beckmann 分布)
    /// </summary>
    /// <param name="alpha">表面粗糙度，控制法线分布的宽度</param>
    /// <param name="theta_h">半角向量与表面法线的夹角，弧度</param>
    /// <returns></returns>
    void beckmannDistribution(float alpha, float theta_h, float& D) {
        // 计算 tan^2(theta_h)
        float tanThetaH2 = tan(theta_h) * tan(theta_h);

        // 计算贝克曼分布
        D = exp(-tanThetaH2 / (alpha * alpha)) / (PI * alpha * alpha * cos(theta_h) * cos(theta_h) * cos(theta_h) * cos(theta_h));
    }

    //shadowing-masking term G(i, o, h) 几何项，自遮挡，自阴影
    void shadowingmaskingG(const Vector3f& normal, const Vector3f& inDir, const Vector3f& outDir, const float& roughness, float& G)
    {
        /*Vector3f half = (inDir + outDir).normalized();
        float G1 = 2 * dotProduct(normal, half) * dotProduct(normal, outDir) / dotProduct(outDir, half);
        float G2 = 2 * dotProduct(normal, half) * dotProduct(normal, inDir) / dotProduct(outDir, half);
        G =  std::min(1.0f, std::min(G1, G2));*/
        Vector3f half = (inDir + outDir).normalized();
        float G1 = 2 * dotProduct(normal, half) * dotProduct(normal, inDir) / dotProduct(inDir, half);
        float G2 = 2 * dotProduct(normal, half) * dotProduct(normal, outDir) / dotProduct(inDir, half);
        G = std::min(1.0f, std::min(G1, G2));
    }

public:
    MaterialType m_type;
    //Vector3f m_color;
    Vector3f m_emission;
    float ior;
    Vector3f Kd, Ks;
    float specularExponent;
    //Texture tex;

    inline Material(MaterialType t=DIFFUSE, Vector3f e=Vector3f(0,0,0));
    inline MaterialType getType();
    //inline Vector3f getColor();
    inline Vector3f getColorAt(double u, double v);
    inline Vector3f getEmission();
    inline bool hasEmission();

    // sample a ray by Material properties
    inline Vector3f sample(const Vector3f &wi, const Vector3f &N);
    // given a ray, calculate the PdF of this ray
    inline float pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N);
    // given a ray, calculate the contribution of this ray
    inline Vector3f eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N);

};

Material::Material(MaterialType t, Vector3f e){
    m_type = t;
    //m_color = c;
    m_emission = e;
}

MaterialType Material::getType(){return m_type;}
///Vector3f Material::getColor(){return m_color;}
Vector3f Material::getEmission() {return m_emission;}
bool Material::hasEmission() {
    if (m_emission.norm() > EPSILON) return true;
    else return false;
}

Vector3f Material::getColorAt(double u, double v) {
    return Vector3f();
}


Vector3f Material::sample(const Vector3f &wi, const Vector3f &N){
    switch(m_type){
        case DIFFUSE:
        {
            // uniform sample on the hemisphere
            float x_1 = get_random_float(), x_2 = get_random_float();
            float z = std::fabs(1.0f - 2.0f * x_1);
            float r = std::sqrt(1.0f - z * z), phi = 2 * M_PI * x_2;
            Vector3f localRay(r*std::cos(phi), r*std::sin(phi), z);
            return toWorld(localRay, N);
            
            break;
        }
        case MICROFACET:
        {
            //return Vector3f(0.0f);
            float x_1 = get_random_float(), x_2 = get_random_float();
            float z = std::fabs(1.0f - 2.0f * x_1);
            float r = std::sqrt(1.0f - z * z), phi = 2 * M_PI * x_2;
            Vector3f localRay(r * std::cos(phi), r * std::sin(phi), z);
            return toWorld(localRay, N);
            break;               
        }
    }
}

float Material::pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N){
    switch(m_type){
        case DIFFUSE:
        {
            // uniform sample probability 1 / (2 * PI)
            if (dotProduct(wo, N) > 0.0f)
                return 0.5f / M_PI;
            else
                return 0.0f;
            break;
        }
        case MICROFACET:
        {
            if (dotProduct(wo, N) > 0.0f)
                return 0.5f / M_PI;
            else
                return 0.0f;
            break;
        }
    }
}

//用于计算在给定的入射光线、出射光线和法线情况下，漫反射材质的光照贡献 即BRDF项，比较简略
Vector3f Material::eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N){
    switch(m_type){
        case DIFFUSE:
        {
            // calculate the contribution of diffuse   model
            float cosalpha = dotProduct(N, wo);
            if (cosalpha > 0.0f) {
                Vector3f diffuse = Kd / M_PI;
                return diffuse;
            }
            else
                return Vector3f(0.0f);
            break;
        }
        case MICROFACET:
        {
            float cosalpha = dotProduct(N, wo);
            if (cosalpha > 0.0f)
            {
                float F, D, G;
                float eta = 1.5;      //介质η值 折射率
                float roughness = 0.35f;
                
                fresnel(wi, N, eta, F);
                Vector3f halfVec = (-wi + wo).normalized();       //半程向量
                float rad_h_n = angleBetween(halfVec, N);
                //distributionH(N, halfVec, roughness, D);
                beckmannDistribution(roughness, rad_h_n, D);
                shadowingmaskingG(N, -wi, wo, roughness, G);

                float specular = F * G *D / std::max((4 * std::max(dotProduct(N, -wi), 0.0f) * std::max(dotProduct(N, wo), 0.0f)), 0.001f);
                Vector3f diffuse = 1.0f / M_PI;
                return Ks * specular + (1 - F) * Kd * diffuse;
            }
            else
            {
                return Vector3f(0.0f);
            }
        }
    }
}

#endif //RAYTRACING_MATERIAL_H
