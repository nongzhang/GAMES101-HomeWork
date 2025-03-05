#include "Renderer.hpp"
#include "Scene.hpp"
#include "Triangle.hpp"
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
    Scene scene(1280, 960);

    string obj_path = "\\..\\..\\models\\bunny\\bunny.obj";
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


    //MeshTriangle bunny("../models/bunny/bunny.obj");
    MeshTriangle bunny(obj_path);

    scene.Add(&bunny);
    scene.Add(std::make_unique<Light>(Vector3f(-20, 70, 20), 1));
    scene.Add(std::make_unique<Light>(Vector3f(20, 70, 20), 1));
    scene.buildBVH();

    Renderer r;

    auto start = std::chrono::system_clock::now();               //SAH加速做效率比较
    r.Render(scene);
    auto stop = std::chrono::system_clock::now();

    std::cout << "Render complete: \n";
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::hours>(stop - start).count() << " hours\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::minutes>(stop - start).count() << " minutes\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << " seconds\n";

    return 0;
}