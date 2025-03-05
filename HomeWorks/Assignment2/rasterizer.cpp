// clang-format off
//
// Created by goksu on 4/6/19.
//

#include <algorithm>
#include <vector>
#include "rasterizer.hpp"
#include <opencv2/opencv.hpp>
#include <math.h>
using namespace std;


rst::pos_buf_id rst::rasterizer::load_positions(const std::vector<Eigen::Vector3f> &positions)
{
    auto id = get_next_id();
    pos_buf.emplace(id, positions);

    return {id};
}

rst::ind_buf_id rst::rasterizer::load_indices(const std::vector<Eigen::Vector3i> &indices)
{
    auto id = get_next_id();
    ind_buf.emplace(id, indices);

    return {id};
}

rst::col_buf_id rst::rasterizer::load_colors(const std::vector<Eigen::Vector3f> &cols)
{
    auto id = get_next_id();
    col_buf.emplace(id, cols);

    return {id};
}

auto to_vec4(const Eigen::Vector3f& v3, float w = 1.0f)
{
    return Vector4f(v3.x(), v3.y(), v3.z(), w);
}


static bool insideTriangle(int x, int y, const Vector3f* _v)
{   
    // TODO : Implement this function to check if the point (x, y) is inside the triangle represented by _v[0], _v[1], _v[2]
    //�������ε���������Ϊp0,p1,p2,�ֱ��ǲ����е�_v[0]��_v[1]��_v[2]�������ص����ĵ�ΪQ(x,y)�����ÿγ��еķ��Ŵ�����ֹ���ƫ��
    //x��y�����͸�Ϊfloat��_v�������Ǹ�ָ�룬�ɵ�������ʹ��
    //������˿���������ʽչ������
    Vector2f p0p1, p1p2, p2p0;
    p0p1 << _v[1].x() - _v[0].x(), _v[1].y() - _v[0].y();
    p1p2 << _v[2].x() - _v[1].x(), _v[2].y() - _v[1].y();
    p2p0 << _v[0].x() - _v[2].x(), _v[0].y() - _v[2].y();

    Vector2f p0Q, p1Q, p2Q;
    p0Q << x - _v[0].x(), y - _v[0].y();
    p1Q << x - _v[1].x(), y - _v[1].y();
    p2Q << x - _v[2].x(), y - _v[2].y();

    float z0, z1, z2;
    z0 = p0p1.x() * p0Q.y() - p0Q.x() * p0p1.y();
    z1 = p1p2.x() * p1Q.y() - p1Q.x() * p1p2.y();
    z2 = p2p0.x() * p2Q.y() - p2Q.x() * p2p0.y();

    return (z0 > 0 && z1 > 0 && z2 > 0) || (z0 < 0 && z1 < 0 && z2 < 0);
}

static std::tuple<float, float, float> computeBarycentric2D(float x, float y, const Vector3f* v)
{
    float c1 = (x*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*y + v[1].x()*v[2].y() - v[2].x()*v[1].y()) / (v[0].x()*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*v[0].y() + v[1].x()*v[2].y() - v[2].x()*v[1].y());
    float c2 = (x*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*y + v[2].x()*v[0].y() - v[0].x()*v[2].y()) / (v[1].x()*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*v[1].y() + v[2].x()*v[0].y() - v[0].x()*v[2].y());
    float c3 = (x*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*y + v[0].x()*v[1].y() - v[1].x()*v[0].y()) / (v[2].x()*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*v[2].y() + v[0].x()*v[1].y() - v[1].x()*v[0].y());
    return {c1,c2,c3};
}

void rst::rasterizer::draw(pos_buf_id pos_buffer, ind_buf_id ind_buffer, col_buf_id col_buffer, Primitive type)
{
    auto& buf = pos_buf[pos_buffer.pos_id];       //�洢�˶���λ�õĻ�������
    auto& ind = ind_buf[ind_buffer.ind_id];          //�洢�������ζ��������Ļ�������
    auto& col = col_buf[col_buffer.col_id];           //�洢����ɫ�Ļ�����

    float f1 = (50 - 0.1) / 2.0;
    float f2 = (50 + 0.1) / 2.0;

    Eigen::Matrix4f mvp = projection * view * model;
    for (auto& i : ind)
    {
        Triangle t;
        Eigen::Vector4f v[] = {
                mvp * to_vec4(buf[i[0]], 1.0f),              //ͨ������ i ���������ε��������㡣��ÿ���������MVP����任
                mvp * to_vec4(buf[i[1]], 1.0f),
                mvp * to_vec4(buf[i[2]], 1.0f)
        };
        //Homogeneous division  ����ָ�,���������ָ�����׼��3D������
        for (auto& vec : v) {
            vec /= vec.w();
        }
        //Viewport transformation    �ӿڱ任  ��ϰ֪ʶ���ӿڱ任����Ϊ��
        //Mviewport << width/2, 0, 0, width/2,
         //                       0, height / 2, 0, height / 2,
         //                       0, 0, 1, 0
         //                       0, 0, 0, 1;
        for (auto & vert : v)
        {
            vert.x() = 0.5*width*(vert.x()+1.0);
            vert.y() = 0.5*height*(vert.y()+1.0);
            vert.z() = vert.z() * f1 + f2;                    //ִ�� �ӿڱ任 ʱ������ȣ�z ���꣩���е����ź�ƽ�Ʋ�����
        }

        //���������εĶ��㣬���任��Ķ������괫�ݸ������ζ��� t
        for (int i = 0; i < 3; ++i)
        {
            t.setVertex(i, v[i].head<3>());
            /*t.setVertex(i, v[i].head<3>());
            t.setVertex(i, v[i].head<3>());*/
        }

        auto col_x = col[i[0]];
        auto col_y = col[i[1]];
        auto col_z = col[i[2]];

        //����ɫ������ col �л�ȡÿ���������ɫ��������Щ��ɫֵӦ�õ������εĶ�����
        t.setColor(0, col_x[0], col_x[1], col_x[2]);
        t.setColor(1, col_y[0], col_y[1], col_y[2]);
        t.setColor(2, col_z[0], col_z[1], col_z[2]);

        rasterize_triangle(t);
    }
}

void rst::rasterizer::enable_MASS(int MSrate = 2)
{
    enableMSAA = true;
    MSAArate = MSrate;
    ms_depth_buf.resize(width * height * MSAArate * MSAArate);
    ms_frame_buf.resize(width * height * MSAArate * MSAArate);
    fill(ms_frame_buf.begin(), ms_frame_buf.end(), Eigen::Vector3f{ 0, 0, 0 });
    fill(ms_depth_buf.begin(), ms_depth_buf.end(), numeric_limits<float>::infinity());
}

void rst::rasterizer::set_ms_pixel(const Eigen::Vector3f& point, const Eigen::Vector3f& color)
{
    int index = get_ms_index(point.x(), point.y());
    ms_frame_buf[index] = color;
}

void rst::rasterizer::set_ms_depth(const Eigen::Vector3f& point)
{
    int index = get_ms_index(point.x(), point.y());
    ms_depth_buf[index] = point.z();
}


float minn(float a, float b, float c) {
    return min(a, min(b, c));
}
float maxx(float a, float b, float c) {
    return max(a, max(b, c));
}


//Screen space rasterization
void rst::rasterizer::rasterize_triangle(const Triangle& t) {
    auto v = t.toVector4();
    
    // TODO : Find out the bounding box of current triangle.
    // iterate through the pixel and find if the current pixel is inside the triangle
    //�������������������x,y�������Сֵ���Թ���AABB��Χ��
    float xmin = min(min(v[0].x(), v[1].x()), v[2].x());
    float xmax = max(max(v[0].x(), v[1].x()), v[2].x());
    float ymin = min(min(v[0].y(), v[1].y()), v[2].y());
    float ymax = max(max(v[0].y(), v[1].y()), v[2].y());

    xmin = floor(xmin);
    xmax = ceil(xmax);
    ymin = floor(ymin);
    ymax = ceil(ymax);

    // If so, use the following code to get the interpolated z value.
    //auto[alpha, beta, gamma] = computeBarycentric2D(x, y, t.v);
    //float w_reciprocal = 1.0/(alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
    //float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
    //z_interpolated *= w_reciprocal;

    // TODO : set the current pixel (use the set_pixel function) to the color of the triangle (use getColor function) if it should be painted.

    //�����ϵ�α��������
	/*for (int x = 0; x < xmax; x++)
	{
		for (int y = 0; y < ymax; y++)
		{
			image[x][y] = inside(tri, x + 0.5, y + 0.5);
		}
	}*/
    for (int i = xmin; i < xmax; i++)
    {
        for (int j = ymin; j < ymax; j++)
        {

            if (enableMSAA)
            {
                //������ʱ��ԭ����һ�����ر��ֳ�nxn����С�����أ���n=2������ԭ�����ص����½�����Ϊ(2,1),ϸ��֮���ĸ�С���ص����½�����Ϊ
                //(2,1) (2,1.5)  (2.5,1) (2.5, 1.5) ������С�������½�����Ϊ (x+k/n, y+l/n)����0<=k<n, 0<=l<n. Ȼ��С�������ĵ������Ϊ(x+k/n+1/(2n),  y+k/n+1/(2n))

                float n = (float)MSAArate;
                // ����С������ɫ
                vector<Vector3f> pixelCols;
                // ����������ֻ��һ����ɫ����ǰȡ��������Ƶ������
                Vector3f color = t.getColor();
                // �����ص�������ɫ
                Vector3f final_color;
                float final_depth;

                for (int l = 0; l < MSAArate; l++)
                {
                    for (int k = 0; k < MSAArate; k++)
                    {
                        if (insideTriangle(i + k / n + 0.5f / n, j + l / n + 0.5f / n, t.v))
                        {
                            // ͨ�������򻺴�С���ص���ɫ
                            pixelCols.push_back(color);
                            /*final_color /= (n * n);*/
                            auto [alpha, beta, gamma] = computeBarycentric2D(i, j, t.v);
                            float w_reciprocal = 1.0 / (alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
                            float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
                            z_interpolated *= w_reciprocal;
                            if (z_interpolated < ms_depth_buf[get_ms_index(i * MSAArate + k, j * MSAArate + l)])
                            {
                                set_ms_pixel(Vector3f(i * MSAArate + k, j * MSAArate + l, z_interpolated), color);
                                set_ms_depth(Vector3f(i * MSAArate + k, j * MSAArate + l, z_interpolated));
                            }
                        }
                    }
                }
                //�ж�С�����Ƿ����������ڡ����ͨ�����ԣ��͸���n��n��С���ص���ɫ�����ֵ��������ص���Ϣ����ɫȡƽ��ֵ�����ȡ��Сֵ
                if (!pixelCols.empty())
                {
                    for (int l = 0; l < n; l++)
                    {
                        for (int k = 0; k < n; k++)
                        {
                            final_color += ms_frame_buf[get_ms_index(i * MSAArate + k, j * MSAArate + l)];
                            final_depth = min(final_depth, ms_depth_buf[get_ms_index(i * MSAArate + k, j * MSAArate + l)]);
                        }
                    }
                    final_color /= (n * n);
                    
                    //if (z_interpolated < depth_buf[get_index(i, j)])
                    //{
                    //    set_pixel(Vector3f(i, j, z_interpolated), t.getColor());   //����������ɫֵ
                    //    depth_buf[get_index(i, j)] = z_interpolated;               //�������ֵ
                    //}
                    //�����棬������һ��С���ز���
                    set_pixel(Vector3f(i, j, 0), final_color);
                    set_depth(Vector3f(i, j, final_depth));
                    final_color.setZero();
                    final_depth = numeric_limits<float>::infinity();
                    pixelCols.clear();
                }
            }
            else
            {
                if (insideTriangle(i + 0.5, j + 0.5, t.v))
                {
                    auto [alpha, beta, gamma] = computeBarycentric2D(i, j, t.v);
                    float w_reciprocal = 1.0 / (alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
                    float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
                    z_interpolated *= w_reciprocal;

                    if (z_interpolated < depth_buf[get_index(i, j)])
                    {
                        set_pixel(Vector3f(i, j, z_interpolated), t.getColor());   //����������ɫֵ
                        depth_buf[get_index(i, j)] = z_interpolated;               //�������ֵ
                    }
                }
            }
           
        }
    }
}

void rst::rasterizer::set_model(const Eigen::Matrix4f& m)
{
    model = m;
}

void rst::rasterizer::set_view(const Eigen::Matrix4f& v)
{
    view = v;
}

void rst::rasterizer::set_projection(const Eigen::Matrix4f& p)
{
    projection = p;
}

void rst::rasterizer::clear(rst::Buffers buff)
{
    if ((buff & rst::Buffers::Color) == rst::Buffers::Color)
    {
        std::fill(frame_buf.begin(), frame_buf.end(), Eigen::Vector3f{0, 0, 0});
    }
    if ((buff & rst::Buffers::Depth) == rst::Buffers::Depth)
    {
        std::fill(depth_buf.begin(), depth_buf.end(), std::numeric_limits<float>::infinity());
    }
}

rst::rasterizer::rasterizer(int w, int h) : width(w), height(h)
{
    frame_buf.resize(w * h);
    depth_buf.resize(w * h);
}

int rst::rasterizer::get_index(int x, int y)
{
    return (height-1-y)*width + x;
}

int rst::rasterizer::get_ms_index(int x, int y)
{
    return (height * MSAArate - 1 - y) * width * MSAArate + x;
}

void rst::rasterizer::set_pixel(const Eigen::Vector3f& point, const Eigen::Vector3f& color)
{
    //old index: auto ind = point.y() + point.x() * width;
    //��ά����ֻʹ����x,y,���¸������Ӧ��֡������ɫֵ
    auto ind = (height-1-point.y())*width + point.x();
    frame_buf[ind] = color;

}

void rst::rasterizer::set_depth(const Eigen::Vector3f& point)
{
    int index = get_index(point.x(), point.y());
    depth_buf[index] = point.z();
}

// clang-format on