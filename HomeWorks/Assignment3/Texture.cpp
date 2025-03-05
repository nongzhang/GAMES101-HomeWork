#include "Texture.hpp"
//
// Created by LEI XU on 4/27/19.
//
using namespace std;
Eigen::Vector3f Texture::getColorBilinear(float u, float v)
{
    //u 和 v 是纹理坐标，通常范围是 [0, 1]。
    //_img 和 v_img 将纹理坐标转换为图像的像素坐标。注意，v 需要反转，因为纹理坐标的 v=0 通常对应于图像的底部，而在图像中，v=0 通常对应于顶部。
    float u_img = u * width;
    float v_img = (1 - v) * height;
    int u_center = round(u_img);       //u_center 和 v_center 是 (u_img, v_img) 的四舍五入值，用来表示最接近的像素中心
    int v_center = round(v_img);

    //u0, u1, v0, v1 分别是围绕 (u_center, v_center) 的四个相邻像素坐标。clamp 函数保证像素索引在有效范围内，防止超出边界。
    int u0 = clamp(u_center - 1, 0, width);
    int u1 = clamp(u_center + 1, 0, width);
    int v0 = clamp(v_center - 1, 0, width);
    int v1 = clamp(v_center + 1, 0, width);

    //通过 cv::Vec3b，获取了四个相邻像素的颜色，存储在 colors 数组中。
    array<cv::Vec3b, 4> colors;
    colors[0] = image_data.at<cv::Vec3b>(v0, u0);
    colors[1] = image_data.at<cv::Vec3b>(v1, u0);
    colors[2] = image_data.at<cv::Vec3b>(v0, u1);
    colors[3] = image_data.at<cv::Vec3b>(v1, u1);

    //这里的 s 和 t 是对纹理坐标的偏移量，用于双线性插值的计算。s 和 t 计算了当前纹理坐标与最接近像素的偏差
    float s = u_img - clamp(u_center - 0.5, 0.5, width - 0.5);
    float t = v_img - clamp(v_center - 0.5, 0.5, height - 0.5);
    // TODO: Implement Bilinear Interpolation Sampler

    //auto color = image_data.at<cv::Vec3b>(v_img, u_img);
    //使用 lerp 函数进行双线性插值。lerp 是线性插值的缩写。这里首先在 t 方向对 colors[0] 和 colors[1] 进行插值，
    // 再对 colors[2] 和 colors[3] 进行插值。然后，使用 s 对这两次插值结果进行进一步的插值，得到最终的颜色
    cv::Vec3b color = lerp(lerp(colors[0], colors[1], t), lerp(colors[2], colors[3], t), s);
    return Eigen::Vector3f(color[0], color[1], color[2]);
}

cv::Vec3b Texture::lerp(cv::Vec3b color0, cv::Vec3b color1, float t)
{
    return color0 * (1 - t) + color1 * t;
    /*cv::Vec3b result_color;
    for (int i = 0; i < 3; i++)
    {
        result_color[i] = color0[i] + (color1[i] - color0[i]) * t;
    }
    return result_color;*/
}
