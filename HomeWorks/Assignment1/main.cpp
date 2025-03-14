#include "Triangle.hpp"
#include "rasterizer.hpp"
#include <eigen3/Eigen/Eigen>
#include <iostream>
#include <opencv2/opencv.hpp>

constexpr double MY_PI = 3.1415926;

Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1, 0, 0, -eye_pos[0], 0, 1, 0, -eye_pos[1], 0, 0, 1,
        -eye_pos[2], 0, 0, 0, 1;

    view = translate * view;

    return view;
}

//Model Matrix
Eigen::Matrix4f get_model_matrix(float rotation_angle)
{
    Eigen::Matrix4f model = Eigen::Matrix4f::Identity();

    // TODO: Implement this function
    // Create the model matrix for rotating the triangle around the Z axis.
    // Then return it.
    Matrix4f rotationZ;
    float rad = rotation_angle / 180.0 * MY_PI;
    rotationZ << cos(rad), -sin(rad), 0, 0,
                        sin(rad), cos(rad), 0, 0,
                        0, 0, 1, 0,
                        0, 0, 0, 1;
    model = rotationZ * model;
    return model;
}
//Projection Matrix
Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio,
                                      float zNear, float zFar)
{
    Eigen::Matrix4f projection = Eigen::Matrix4f::Identity();

    //透视投影矩阵=正交投影矩阵 * 变换矩阵(Perspective->Orthographic)
    //参数eye_fov就是fovY,以近平面中心点为例 tan(fovY/2) = t/abs(n) ，aspect_ratio = r/t，而且由于标准立方体是对称的，l=-r, b=-t, 参数zNear就是n, 参数zFar就是f
    float t = atan(eye_fov / 360.0 * acos(-1)) * (-zNear);
    float b = -t;
    float l = -aspect_ratio * t;
    float r = -l;

    Matrix4f scaleOrtho;
    Matrix4f transformOrtho;
    Matrix4f ortho;
    Matrix4f perspToOrtho;

    scaleOrtho << 2 / (r - l), 0, 0, 0,
                                0, 2 / (t - b), 0, 0,
                                0, 0, 2 / (zNear - zFar), 0,
                                0, 0, 0, 1;

    transformOrtho << 1, 0, 0, -(r + l) / 2,
                                    0, 1, 0, -(t + b) / 2,
                                    0, 0, 1, -(zNear + zFar) / 2,
                                    0, 0, 0, 1;

    perspToOrtho << zNear, 0, 0, 0,
                                0, zNear, 0, 0,
                                0, 0, zNear + zFar, -zNear * zFar,
                                0, 0, 1, 0;

    ortho = scaleOrtho * transformOrtho;
    projection = ortho * perspToOrtho * projection;

    return projection;
}

int main(int argc, const char** argv)
{
    float angle = 0;
    bool command_line = false;
    std::string filename = "output.png";

    if (argc >= 3) {
        command_line = true;
        angle = std::stof(argv[2]); // -r by default
        if (argc == 4) {
            filename = std::string(argv[3]);
        }
    }

    rst::rasterizer r(700, 700);

    Eigen::Vector3f eye_pos = {0, 0, 5};

    std::vector<Eigen::Vector3f> pos{{2, 0, -2}, {0, 2, -2}, {-2, 0, -2}};

    std::vector<Eigen::Vector3i> ind{{0, 1, 2}};

    auto pos_id = r.load_positions(pos);
    auto ind_id = r.load_indices(ind);

    int key = 0;
    int frame_count = 0;

    if (command_line) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);

        cv::imwrite(filename, image);

        return 0;
    }

    while (key != 27) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);

        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::imshow("image", image);
        key = cv::waitKey(10);

        std::cout << "frame count: " << frame_count++ << '\n';

        if (key == 'a') {
            angle += 10;
        }
        else if (key == 'd') {
            angle -= 10;
        }
    }

    return 0;
}
