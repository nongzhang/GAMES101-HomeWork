//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    // TO DO Implement Path Tracing Algorithm here

    //1. 求射线与场景BVH的交点，没有交点则直接返回(0,0,0)
    Vector3f vector(0);
    Intersection intersection = Scene::intersect(ray);
    if (intersection.happened)
    {
        //2. 对光源采样，获取光源位置与光源pdf(概率密度函数函数)，同时获取交点处材质信息
        Vector3f L_dir(0.0), L_indir(0.0);
        Intersection lightPos;
        float lightPdf;
        sampleLight(lightPos, lightPdf);
        Vector3f N = intersection.normal;
        Material* material = intersection.m;

        //从目标交点向光源发出射线，并求出光源射线与场景内物体的交点
        Vector3f x = intersection.coords;    //着色点
        Vector3f light_p_x = (lightPos.coords - x).normalized();        //从光源指向着色点的单位向量，因为方向是向外，所以用lighPos-x

        Ray ray_p_x(x, light_p_x);
        Intersection intersection_p_x = Scene::intersect(ray_p_x);
        if (intersection_p_x.happened && intersection_p_x.m->hasEmission())
        {
            //GAMES101_Lecture_16 P41中的n'
            Vector3f NN = intersection_p_x.normal;  

            //直接来自光源的贡献 L_dir = L_i * f_r * cos θ * cos θ’ / |x’ - p|^2 / pdf_light
            //material->eval(ray.direction, light_p_x, N)指的是BRDF,代表着色点x处出射方向的Radiance/点x处接收到的Irradiance
            L_dir = lightPos.emit * material->eval(ray.direction, light_p_x, N) * dotProduct(light_p_x, N)
                * dotProduct(-light_p_x, NN) / intersection_p_x.distance / intersection_p_x.distance / lightPdf;

            /*L_dir = intersection_p_x.m->m_emission * material->eval(ray.direction, light_p_x, N) * dotProduct(light_p_x, N)
                * dotProduct(-light_p_x, NN) / intersection_p_x.distance / intersection_p_x.distance / lightPdf;*/

        }

        //使用俄罗斯轮盘赌来计算反射光
        if (get_random_float() <= RussianRoulette)
        {
            //wi = sample(wo, N)  
            Vector3f inputDir = material->sample(Vector3f(0), N).normalized();
            Ray inputDirRay(x, inputDir);
            Intersection intersectionInputDirRay = Scene::intersect(inputDirRay);

            // If ray r hit a non-emitting object at q
            if (intersectionInputDirRay.happened && !intersectionInputDirRay.m->hasEmission())
            {
                // L_indir = shade(q, -wi) * f_r * cos θ / pdf_hemi / P_RR     这里的wi是光源
                //L_indir = shade(q, wi) * eval(wo, wi, N) * dot(wi, N) / pdf(wo, wi, N ) / RussianRoulette
                //material->pdf(ray.direction, inputDir, N)  半球的概率密度函数 pdf_hemi
                
                L_indir = castRay(inputDirRay, depth + 1) * material->eval(ray.direction, inputDir, N) * dotProduct(inputDir, N) / material->pdf(ray.direction, inputDir, N) / RussianRoulette;
                //L_indir = castRay(inputDirRay, depth + 1) * material->eval(ray.direction, inputDir, N) * dotProduct(inputDir, N) / material->pdf(Vector3f(0), inputDir, N) / RussianRoulette;

            }
        }

        return material->getEmission() + L_dir + L_indir;
    }
    return vector;
}

