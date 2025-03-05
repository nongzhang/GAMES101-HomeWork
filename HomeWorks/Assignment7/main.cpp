#include "Renderer.hpp"
#include "Scene.hpp"
#include "Triangle.hpp"
#include "Sphere.hpp"
#include "Vector.hpp"
#include "global.hpp"
#include <chrono>
#include "Windows.h"
using namespace std;

// In the main function of the program, we create the scene (create objects and
// lights) as well as set the options for the render (image width and height,
// maximum recursion depth, field-of-view, etc.). We then call the render
// function().
int main(int argc, char** argv)
{

    // Change the definition here to change resolution
    Scene scene(784, 784);

    Material* red = new Material(DIFFUSE, Vector3f(0.0f));
    red->Kd = Vector3f(0.63f, 0.065f, 0.05f);
    Material* green = new Material(DIFFUSE, Vector3f(0.0f));
    green->Kd = Vector3f(0.14f, 0.45f, 0.091f);
    Material* yellow = new Material(DIFFUSE, Vector3f(0.0f));
    yellow->Kd = Vector3f(0.60f, 0.60f, 0.185f);
    Material* white = new Material(DIFFUSE, Vector3f(0.0f));
    white->Kd = Vector3f(0.725f, 0.71f, 0.68f);

    Material* microfacet = new Material(MICROFACET, Vector3f(0.0f));    
    /*microfacet->Ks = Vector3f(0.5f);
    microfacet->Kd = Vector3f(0.5f);*/
    microfacet->Ks = Vector3f(0.8, 0.8, 0.8);
    microfacet->Kd = Vector3f(0.6, 0.6, 0.6);
    //microfacet->ior = 12.85;
    microfacet->ior = 2.42;

    //光源材质 light 是通过多个颜色通道的加权计算来创建的，模拟了一个强度不同的光源。
    Material* light = new Material(DIFFUSE, (8.0f * Vector3f(0.747f+0.058f, 0.747f+0.258f, 0.747f) + 15.6f * Vector3f(0.740f+0.287f,0.740f+0.160f,0.740f) + 18.4f *Vector3f(0.737f+0.642f,0.737f+0.159f,0.737f)));
    light->Kd = Vector3f(0.65f);


    string obj_path = "\\..\\..\\models\\cornellbox\\";
    char pathBuf[MAX_PATH];
    char* p;
    if (GetModuleFileNameA(NULL, pathBuf, MAX_PATH))
    {
        p = strrchr(pathBuf, '\\');
        if (p)
        {
            *p = '\0';
            string exe_path = pathBuf;
            obj_path = exe_path + obj_path;
        }
    }

    //直接使用构造函数进行初始化是一种简洁且常见的做法,等同于
    //MeshTriangle floor;  // 默认构造对象
    //floor = MeshTriangle("../models/cornellbox/floor.obj", white);  // 调用带参数的构造函数
    MeshTriangle floor(obj_path + "floor.obj", white);
    MeshTriangle shortbox(obj_path + "shortbox.obj", white);
    MeshTriangle tallbox(obj_path + "tallbox.obj", white);
    MeshTriangle left(obj_path + "left.obj", red);
    MeshTriangle right(obj_path + "right.obj", green);
    MeshTriangle light_(obj_path + "light.obj", light);
   /* MeshTriangle floor("../models/cornellbox/floor.obj", white);
    MeshTriangle shortbox("../models/cornellbox/shortbox.obj", white);
    MeshTriangle tallbox("../models/cornellbox/tallbox.obj", white);
    MeshTriangle left("../models/cornellbox/left.obj", red);
    MeshTriangle right("../models/cornellbox/right.obj", green);
    MeshTriangle light_("../models/cornellbox/light.obj", light);*/

    Sphere sphere(Vector3f(150, 100, 300), 100, microfacet);
    //Sphere sphere(Vector3f(150, 100, 200), 100, microfacet);

    scene.Add(&floor);
    scene.Add(&sphere);
    /*scene.Add(&shortbox);
    scene.Add(&tallbox);*/
    scene.Add(&left);
    scene.Add(&right);
    scene.Add(&light_);

    scene.buildBVH();

    Renderer r;

    auto start = std::chrono::system_clock::now();
    r.Render(scene);
    auto stop = std::chrono::system_clock::now();

    std::cout << "Render complete: \n";
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::hours>(stop - start).count() << " hours\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::minutes>(stop - start).count() << " minutes\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << " seconds\n";

    return 0;
}