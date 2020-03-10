#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "svpng.inc"

const double PI = 3.1415926;
const double ONE_PI = 1.0 / PI;

// 产生 0-1 的随机数
double erand48(unsigned short xsubi[3]) 
{
	return (double) rand() / (double) RAND_MAX;
} 


// 向量：表示坐标、颜色 
struct Vec 
{

  double x, y, z;                  // 颜色 (r,g,b)
  
  Vec(double x_=0, double y_=0, double z_=0) { x=x_; y=y_; z=z_; }
  
  Vec operator+(const Vec &b) const { return Vec(x+b.x,y+b.y,z+b.z); }
  
  Vec operator-(const Vec &b) const { return Vec(x-b.x,y-b.y,z-b.z); }
  
  Vec operator*(double b) const { return Vec(x*b,y*b,z*b); }
  
  Vec mult(const Vec &b) const { return Vec(x*b.x,y*b.y,z*b.z); }
  
  Vec& norm(){ return *this = *this * (1/sqrt(x*x+y*y+z*z)); }
  
  double dot(const Vec &b) const { return x*b.x+y*b.y+z*b.z; } // 点乘 
  
  Vec operator%(Vec&b) { return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x); } // 叉乘 
  
};


// 射线：P(t) = O + t D 
struct Ray 
{ 
	Vec o, d; 
	Ray(Vec o_, Vec d_) : o(o_), d(d_) {} 
};


// 材质类型
enum Refl_type { DIFF, SPEC, REFR };


// 球（唯一支持的物体）
struct Sphere 
{
    double rad;       // 半径
    Vec pos, emi, col;      // 位置, 自发光, 颜色
    Refl_type reflType;      // 反射 (DIFFuse, SPECular, REFRactive)

    Sphere(double rad_, Vec pos_, Vec emi_, Vec col_, Refl_type reflType_)
        : rad(rad_), pos(pos_), emi(emi_), col(col_), reflType(reflType_) {}

    // 球与射线相交测试：返回交点 t，0表示
    double intersect(const Ray &r) const 
    {
        // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        Vec op = pos - r.o;
        double t, 
        eps=1e-4, b=op.dot(r.d), det=b*b-op.dot(op)+rad*rad;
        if (det < 0) 
            return 0; 
        else 
            det=sqrt(det);
        return (t = b-det) > eps ? t : ((t = b + det) > eps ? t : 0);
    }
  
};

// 场景（由很多球组成） 
Sphere spheres[] = {
    //球: 半径, 位置, 自发光, 颜色, 材质
    Sphere(1e5, Vec(1e5 + 1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF),//Left
    Sphere(1e5, Vec(-1e5 + 99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),//Rght
    Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF),//Back
    Sphere(1e5, Vec(50,40.8,-1e5 + 170), Vec(),Vec(),           DIFF),//Frnt
    Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF),//Botm
    Sphere(1e5, Vec(50,-1e5 + 81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top
    Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.999, SPEC),//Mirr
    Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, REFR),//Glas
    Sphere(600, Vec(50,681.6 - .27,81.6),Vec(12,12,12),  Vec(), DIFF) //Lite
};

// 球的个数
int numSpheres = sizeof(spheres) / sizeof(Sphere);

// 约束到 0-1 
inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }

// 颜色值转成0-255，其中1/2.2是gamma校正 
inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + 0.5); }

// 遍历场景中的球，找到射线撞击到的那个
inline bool intersect (const Ray &r, double &t, int &id) 
{
  double n = sizeof(spheres) / sizeof(Sphere), d, inf = t = 1e20;
  for(int i = int(n); i--; )
  { 
    if( (d = spheres[i].intersect(r)) && d<t ) 
    {
	  t = d;
      id = i;
	}
  }
  return t < inf;
}


// 计算辐射率
Vec radiance(const Ray &r, int depth, unsigned short *Xi, int E = 1) 
{
    double t; // 交点 t
    int index = 0; // 相交物体的索引
    if (!intersect(r, t, index))
        return Vec(); // 没有与物体相交则返回黑色
  
    const Sphere &obj = spheres[index]; //相交的球
    Vec inter = r.o + r.d * t; // 交点
    Vec n = (inter - obj.pos).norm(); // 法线 = 交点 - 球心
    Vec nl = n.dot(r.d) < 0 ? n : n * -1; // 光线和法线的点乘小于 0 说明光线在物体外部, 否则法线要取反指向内部
    Vec f = obj.col; // 物体颜色（反射率）

    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; // 俄罗斯轮盘的p取反射率最大分量
    p = 0.25 + 0.5 * p;
    if (++depth > 5)
    {
        if (erand48(Xi) < p) // 俄罗斯轮盘（Russian Roulette）
            f = f * (1 / p);
        else 
            return obj.emi * E; // 结束递归，返回物体的自发光
    }

    if (obj.reflType == DIFF) // 理想漫反射
    {
        Vec w = nl; // 法线
        Vec u = ((fabs(w.x) > 0.1 ? Vec(0, 1) : Vec(1)) % w).norm();
        Vec v = w % u;

        double r1 = 2 * PI * erand48(Xi); // 采样随机角度 0-2π
        double r2 = erand48(Xi), r2s=sqrt(r2); // 采样随机距离
        Vec d = ( u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2) ).norm(); // 采样半球

        // 遍历所有的光源
        Vec emi;
        for (int i = 0; i < numSpheres; i++) 
        {
            const Sphere &sphe = spheres[i];
            if (sphe.emi.x <= 0 && sphe.emi.y <= 0 && sphe.emi.z <= 0) // 略过没有自发光的非光源
                continue;

            // 构建坐标系
            Vec sw = sphe.pos - inter;
            Vec su = ((fabs(sw.x) > 0.1 ? Vec(0, 1) : Vec(1)) % sw).norm();
            Vec sv = sw % su;

            // 在球表示的光源上采样
            double cos_a_max = sqrt(1 - sphe.rad * sphe.rad / (inter - sphe.pos).dot(inter - sphe.pos));
            double eps1 = erand48(Xi), eps2 = erand48(Xi);
            double cos_a = 1 - eps1 + eps1 * cos_a_max;
            double sin_a = sqrt(1 - cos_a * cos_a);
            double phi = 2 * PI * eps2;
            Vec inte2light = su * cos(phi) * sin_a + sv * sin(phi) * sin_a + sw * cos_a;
            inte2light.norm();

            // 发射阴影光线（shadow ray）
            if (intersect(Ray(inter, inte2light), t, index) && index == i)
            {
                double omega = 2 * PI * (1 - cos_a_max);
                emi = emi + f.mult(sphe.emi * inte2light.dot(nl) * omega) * ONE_PI;  // Lambertian 漫反射
            }
        }

        return obj.emi + emi + f.mult( radiance(Ray(inter,d), depth, Xi, 0) );
    } 
    else if (obj.reflType == SPEC) // 理想镜面反射
    {
        // 镜面反射光线 R = D - 2 (N · D) N
        return obj.emi + f.mult( radiance( Ray(inter, r.d - n * 2 * n.dot(r.d)), depth, Xi) );
    }
    else if (obj.reflType == REFR) // 理想介质折射
    {
        Ray reflRay(inter, r.d - n * 2 * n.dot(r.d));
        bool into = n.dot(nl) > 0; // 判断光线是不是从外部进入内部
        double nc = 1, nt = 1.5, nnt = into ? nc/nt : nt/nc; // 空气的折射率约等于1,水的折射率1.33,玻璃的折射率1.4-1.6
        double ddn = r.d.dot(nl), cos2t;

        if ( (cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0 ) // 全内反射
        {
            return obj.emi + f.mult(radiance(reflRay,depth,Xi));
        }

        Vec tdir = (r.d * nnt - n * ( (into ? 1 : -1) * (ddn * nnt + sqrt(cos2t))) ).norm(); // 计算折射光线
    
        double a = nt - nc, b = nt + nc;
        double F0 = a * a / (b * b); // F0 = ((n1-n2)/(n1+n2))^2 表示入射光垂直于表面时，菲涅尔反射率的值
    
        double c = 1 - (into ? -ddn : tdir.dot(n));
        double Re = F0 + (1 - F0) * c * c * c * c * c; // 光线反射部分所占比例（Fresnel-Schlick 近似法）
        double Tr = 1 - Re; // 光线折射部分所占比例
        double P = 0.25 + 0.5 * Re; // TODO：这么取的意义，随便取的？
        double RP = Re / P, TP = Tr / (1 - P);
        return obj.emi + 
            f.mult(depth > 2 ?
            (erand48(Xi) < P ? radiance(reflRay, depth, Xi) * RP : radiance(Ray(inter,tdir), depth, Xi) * TP) :
            radiance(reflRay, depth, Xi) * Re + radiance(Ray(inter, tdir), depth, Xi) * Tr);
    }
}


// 主函数 
int main(int argc, char *argv[])
{

    int w = 1024, h = 768; // 图像大小 
    int samps = argc == 2 ? atoi(argv[1]) / 4 : 1; // 采样数

    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // 相机位置和方向 
    Vec cx = Vec(w * 0.5135 / h);
    Vec cy = (cx % cam.d).norm() * 0.5135;

    Vec r;
    unsigned char *img = new unsigned char[w * h * 3]; // 保存图像

    int k = 0;
#pragma omp parallel for schedule(dynamic, 1) private(r) // CPU并行：OpenMP
    for (int y = 0; y < h; y++) // 遍历图像行
    {
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100.*y / (h - 1));
        unsigned short Xi[3] = { 0, 0, y * y * y };
        for (unsigned short x = 0; x < w; x++)   // 遍历图像列
        {
            Vec color;
            int inde = (h - y - 1) * w + x; // 计算当前像素位置
            for (int sy = 0; sy < 2; sy++)     // 遍历 2 x 2 子像素 行
            {
                for (int sx = 0; sx < 2; sx++, r = Vec()) // 遍历 2 x 2 子像素 列
                {
                    for (int s = 0; s < samps; s++)
                    {
                        double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                            cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
                        r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0, Xi) * (1.0 / samps);
                    }
                    color = color + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * 0.25;
                }
            }
            img[inde * 3 + 0] = unsigned char(toInt(color.x)); // r
            img[inde * 3 + 1] = unsigned char(toInt(color.y)); // g
            img[inde * 3 + 2] = unsigned char(toInt(color.z)); // b
        }
    }

    // 结果输出到 PPM 文件中
    svpng(fopen("ismallpt.png", "wb"), w, h, img, 0); // 将结果输出到png图片中
    delete[] img;

}
