#include "Texture.hpp"
//
// Created by LEI XU on 4/27/19.
//
using namespace std;
Eigen::Vector3f Texture::getColorBilinear(float u, float v)
{
    //u �� v ���������꣬ͨ����Χ�� [0, 1]��
    //_img �� v_img ����������ת��Ϊͼ����������ꡣע�⣬v ��Ҫ��ת����Ϊ��������� v=0 ͨ����Ӧ��ͼ��ĵײ�������ͼ���У�v=0 ͨ����Ӧ�ڶ�����
    float u_img = u * width;
    float v_img = (1 - v) * height;
    int u_center = round(u_img);       //u_center �� v_center �� (u_img, v_img) ����������ֵ��������ʾ��ӽ�����������
    int v_center = round(v_img);

    //u0, u1, v0, v1 �ֱ���Χ�� (u_center, v_center) ���ĸ������������ꡣclamp ������֤������������Ч��Χ�ڣ���ֹ�����߽硣
    int u0 = clamp(u_center - 1, 0, width);
    int u1 = clamp(u_center + 1, 0, width);
    int v0 = clamp(v_center - 1, 0, width);
    int v1 = clamp(v_center + 1, 0, width);

    //ͨ�� cv::Vec3b����ȡ���ĸ��������ص���ɫ���洢�� colors �����С�
    array<cv::Vec3b, 4> colors;
    colors[0] = image_data.at<cv::Vec3b>(v0, u0);
    colors[1] = image_data.at<cv::Vec3b>(v1, u0);
    colors[2] = image_data.at<cv::Vec3b>(v0, u1);
    colors[3] = image_data.at<cv::Vec3b>(v1, u1);

    //����� s �� t �Ƕ����������ƫ����������˫���Բ�ֵ�ļ��㡣s �� t �����˵�ǰ������������ӽ����ص�ƫ��
    float s = u_img - clamp(u_center - 0.5, 0.5, width - 0.5);
    float t = v_img - clamp(v_center - 0.5, 0.5, height - 0.5);
    // TODO: Implement Bilinear Interpolation Sampler

    //auto color = image_data.at<cv::Vec3b>(v_img, u_img);
    //ʹ�� lerp ��������˫���Բ�ֵ��lerp �����Բ�ֵ����д������������ t ����� colors[0] �� colors[1] ���в�ֵ��
    // �ٶ� colors[2] �� colors[3] ���в�ֵ��Ȼ��ʹ�� s �������β�ֵ������н�һ���Ĳ�ֵ���õ����յ���ɫ
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
