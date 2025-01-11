
#include "mrcmx/mrcstack.h"
#include "opts.h"
#include <fstream>
#include <iostream>
#include <vector>
// #include <eigen3/Eigen/Dense>
#include <opencv2/opencv.hpp>
// using namespace Eigen;
#include <opencv2/opencv.hpp>
#include <opencv2/xphoto.hpp>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "../lwt/lwave_test.h"
#include "../filter/filter_prj.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"

#define MILLION 1000000

#define PI_180 0.01745329252f

#ifndef PI
#define PI 3.14159265358979323846
#endif

#define D2R(__ANGLE__) ((__ANGLE__) * PI_180)

struct Coeff
{ // 20个双精度数，前10个是a后10个是b

	double a[10];
	double b[10];
};
/*三元组定义*/
typedef struct
{			 // 都是从0开始
	int r;	 /*行号*/
	int c;	 /*列号*/
	float d; /*元素值*/
} TupNode;

//  /*三元组存储结构顺序表定义*/
//  size_t TS_MAXSIZE=10;
// typedef struct
// {
// 	int rows; // 行数值
// 	int cols; // 列数值
// 	int nums; // 非零元素个数
// 	TupNode data[3000];//这个选多少啊？？能输入的把？
// }TSMatrix;

bool ReadAngles(std::vector<float> &angles, const char *name)
{
	std::ifstream in(name);
	if (!in.good())
	{
		return false;
	}

	while (in.good())
	{
		float val;
		in >> val;
		if (in.fail())
		{
			break;
		}

		angles.push_back(val);
	}
	in.close();
	return true;
}

extern "C" void TranslateAngleToCoefficients(float *angles, Coeff *coeffs, int anglesize)
{
	// coeffs.resize(angles.size());coeffs和angles一样大
	for (int i = 0; i < anglesize; i++)
	{
		memset(coeffs[i].a, 0, sizeof(double) * 10);
		memset(coeffs[i].b, 0, sizeof(double) * 10);
		float beta = D2R(angles[i]);

		coeffs[i].a[0] = 0;			 //
		coeffs[i].a[1] = cos(beta);	 // x
		coeffs[i].a[2] = 0;			 // y
		coeffs[i].a[3] = -sin(beta); // z
		coeffs[i].b[0] = 0;			 //
		coeffs[i].b[1] = 0;			 // x
		coeffs[i].b[2] = 1;			 // y
		coeffs[i].b[3] = 0;			 // z
	}
}

/*solve inverse transfroms defined by Geometry; substitute the inversion into
 * coefficients*/
extern "C" void DecorateCoefficients(Coeff *coeffs, float pitch_angle, float offset, float zshift, int anglesize)
{
	double alpha = -D2R(pitch_angle), beta = D2R(offset), t = -zshift;
	double ca = cos(alpha), sa = sin(alpha), cb = cos(beta), sb = sin(beta);
	double ca2 = ca * ca, sa2 = sa * sa, cb2 = cb * cb, sb2 = sb * sb;

	for (int i = 0; i < anglesize; i++)
	{
		double a[10], b[10];
		memcpy(a, coeffs[i].a, sizeof(double) * 10);
		memcpy(b, coeffs[i].b, sizeof(double) * 10);
		coeffs[i].a[0] = a[0];
		coeffs[i].a[1] = (a[2] * sa * sb + a[3] * ca * sb + a[1] * cb); //*x
		coeffs[i].a[2] = (a[2] * ca - a[3] * sa);						//*y
		coeffs[i].a[3] = (a[3] * ca * cb + a[2] * cb * sa - a[1] * sb); //*z
		coeffs[i].a[4] = (a[4] * ca * cb - a[5] * cb * sa - a[6] * sa2 * sb -
						  2 * a[9] * ca * sa * sb + a[6] * ca2 * sb +
						  2 * a[8] * ca * sa * sb); //*x*y
		coeffs[i].a[5] =
			(a[4] * cb2 * sa - a[5] * ca * sb2 + a[5] * ca * cb2 -
			 2 * a[7] * cb * sb - a[4] * sa * sb2 + 2 * a[9] * ca2 * cb * sb +
			 2 * a[8] * cb * sa2 * sb + 2 * a[6] * ca * cb * sa * sb); //*x*z
		coeffs[i].a[6] =
			(2 * a[8] * ca * cb * sa - a[4] * ca * sb + a[5] * sa * sb +
			 a[6] * ca2 * cb - a[6] * cb * sa2 - 2 * a[9] * ca * cb * sa); //*y*z
		coeffs[i].a[7] =
			(a[7] * cb2 + a[5] * ca * cb * sb + a[9] * ca2 * sb2 +
			 a[8] * sa2 * sb2 + a[4] * cb * sa * sb + a[6] * ca * sa * sb2); //*x^2
		coeffs[i].a[8] = (a[8] * ca2 + a[9] * sa2 - a[6] * ca * sa);		 //*y^2
		coeffs[i].a[9] = (a[7] * sb2 + a[9] * ca2 * cb2 + a[8] * cb2 * sa2 -
						  a[4] * cb * sa * sb + a[6] * ca * cb2 * sa -
						  a[5] * ca * cb * sb); //*z^2

		coeffs[i].b[0] = b[0];
		coeffs[i].b[1] = (b[2] * sa * sb + b[3] * ca * sb + b[1] * cb); //*x
		coeffs[i].b[2] = (b[2] * ca - b[3] * sa);						//*y
		coeffs[i].b[3] = (b[3] * ca * cb + b[2] * cb * sa - b[1] * sb); //*z
		coeffs[i].b[4] = (b[4] * ca * cb - b[5] * cb * sa - b[6] * sa2 * sb -
						  2 * b[9] * ca * sa * sb + b[6] * ca2 * sb +
						  2 * b[8] * ca * sa * sb); //*x*y
		coeffs[i].b[5] =
			(b[4] * cb2 * sa - b[5] * ca * sb2 + b[5] * ca * cb2 -
			 2 * b[7] * cb * sb - b[4] * sa * sb2 + 2 * b[9] * ca2 * cb * sb +
			 2 * b[8] * cb * sa2 * sb + 2 * b[6] * ca * cb * sa * sb); //*x*z
		coeffs[i].b[6] =
			(2 * b[8] * ca * cb * sa - b[4] * ca * sb + b[5] * sa * sb +
			 b[6] * ca2 * cb - b[6] * cb * sa2 - 2 * b[9] * ca * cb * sa); //*y*z
		coeffs[i].b[7] =
			(b[7] * cb2 + b[5] * ca * cb * sb + b[9] * ca2 * sb2 +
			 b[8] * sa2 * sb2 + b[4] * cb * sa * sb + b[6] * ca * sa * sb2); //*x^2
		coeffs[i].b[8] = (b[8] * ca2 + b[9] * sa2 - b[6] * ca * sa);		 //*y^2
		coeffs[i].b[9] = (b[7] * sb2 + b[9] * ca2 * cb2 + b[8] * cb2 * sa2 -
						  b[4] * cb * sa * sb + b[6] * ca * cb2 * sa -
						  b[5] * ca * cb * sb); //*z^2
	}

	// considering z_shift
	for (int i = 0; i < anglesize; i++)
	{
		double a[10], b[10];
		memcpy(a, coeffs[i].a, sizeof(double) * 10);
		memcpy(b, coeffs[i].b, sizeof(double) * 10);

		coeffs[i].a[0] = a[0] + a[3] * t + a[9] * t * t;
		coeffs[i].a[1] = a[1] + a[5] * t;	  //*x
		coeffs[i].a[2] = a[2] + a[6] * t;	  //*y
		coeffs[i].a[3] = a[3] + 2 * a[9] * t; //*z

		coeffs[i].b[0] = b[0] + b[3] * t + b[9] * t * t;
		coeffs[i].b[1] = b[1] + b[5] * t;	  //*x
		coeffs[i].b[2] = b[2] + b[6] * t;	  //*y
		coeffs[i].b[3] = b[3] + 2 * b[9] * t; //*z
	}
}

void (*funL)(const Coeff &, double, double, double, double *);

void WarpPosition(const Coeff &coeff, double X, double Y, double Z, double *n)
{
	n[0] = coeff.a[0] + coeff.a[1] * X + coeff.a[2] * Y + coeff.a[3] * Z +
		   coeff.a[4] * X * Y + coeff.a[5] * X * Z + coeff.a[6] * Y * Z +
		   coeff.a[7] * X * X + coeff.a[8] * Y * Y + coeff.a[9] * Z * Z;
	n[1] = coeff.b[0] + coeff.b[1] * X + coeff.b[2] * Y + coeff.b[3] * Z +
		   coeff.b[4] * X * Y + coeff.b[5] * X * Z + coeff.b[6] * Y * Z +
		   coeff.b[7] * X * X + coeff.b[8] * Y * Y + coeff.b[9] * Z * Z;
}

void LinearPosition(const Coeff &coeff, double X, double Y, double Z,
					double *n)
{
	n[0] = coeff.a[0] + coeff.a[1] * X + coeff.a[2] * Y + coeff.a[3] * Z;
	n[1] = coeff.b[0] + coeff.b[1] * X + coeff.b[2] * Y + coeff.b[3] * Z;
}

void ValCoef(const Point3DF &origin, const Point3D &coord, const Coeff &coeff,
			 Weight *wt)
{
	double x, y;

	double X, Y, Z, n[2];
	X = coord.x - origin.x;
	Y = coord.y - origin.y;
	Z = coord.z - origin.z;

	// 	funL(coeff, X, Y, Z, n);
	// 投影坐标
	n[0] = coeff.a[0] + coeff.a[1] * X + coeff.a[2] * Y + coeff.a[3] * Z;
	n[1] = coeff.b[0] + coeff.b[1] * X + coeff.b[2] * Y + coeff.b[3] * Z;

	x = n[0] + origin.x;
	y = n[1] + origin.y;

	wt->x_min = floor(x);
	wt->y_min = floor(y);

	wt->x_min_del = x - wt->x_min;
	wt->y_min_del = y - wt->y_min;
}

// 这个不是repro 是backpro
void Reproject_admm_htb(const Point3DF &origin, int width, int length, int height, const Coeff &coeff, float *htb, const float *slcdata)
{
	Point3D coord;
	int n;
	size_t volsize = width * length * height;
	for (int z = 0; z < height; z++)
	{

		float *htb_z = htb + z * (size_t)width * length;
		coord.z = z;

		for (int y = 0; y < length; y++)
		{

			float *htb_y = htb_z + y * width;
			coord.y = y;

			for (int x = 0; x < width; x++)
			{

				float *htb_x = htb_y + x;
				coord.x = x;
				Weight wt;
				ValCoef(origin, coord, coeff, &wt);
				if (wt.x_min >= 0 && wt.x_min < width && wt.y_min >= 0 &&
					wt.y_min < length)
				{																	//(x_min, y_min)
					n = wt.x_min + wt.y_min * width;								// index in reproj
					*htb_x += (1 - wt.x_min_del) * (1 - wt.y_min_del) * slcdata[n]; // htb
				}
				if ((wt.x_min + 1) >= 0 && (wt.x_min + 1) < width &&
					wt.y_min >= 0 && wt.y_min < length)
				{										 //(x_min+1, y_min)
					n = wt.x_min + 1 + wt.y_min * width; // index in reproj
					*htb_x += wt.x_min_del * (1 - wt.y_min_del) * slcdata[n];
				}
				if (wt.x_min >= 0 && wt.x_min < width && (wt.y_min + 1) >= 0 &&
					(wt.y_min + 1) < length)
				{										   //(x_min, y_min+1)
					n = wt.x_min + (wt.y_min + 1) * width; // index in reproj
					*htb_x += (1 - wt.x_min_del) * wt.y_min_del * slcdata[n];
				}
				if ((wt.x_min + 1) >= 0 && (wt.x_min + 1) < width &&
					(wt.y_min + 1) >= 0 &&
					(wt.y_min + 1) < length)
				{												 //(x_min+1, y_min+1)
					n = (wt.x_min + 1) + (wt.y_min + 1) * width; // index in reproj
					*htb_x += wt.x_min_del * wt.y_min_del * slcdata[n];
				}
			}
		}
	}
}
// repro+backpro
void Reproject_admm_atax(const Point3DF &origin, int width, int length, int height, const float *x0, Coeff coeffv[],
						 float *atax, int proj)

{
	Point3D coord;
	int n;
	size_t volsize = width * length * height;
	float *ax = (float *)malloc(sizeof(float) * width * length);
	float *w = (float *)malloc(sizeof(float) * width * length);
	int pxsize = width * length;
	float a = 0;

	for (int idx = 0; idx < proj; idx++)
	{
		memset(ax, 0, sizeof(float) * width * length);
		memset(w, 0, sizeof(float) * width * length);
		for (int z = 0; z < height; z++)
		{

			const float *vdrefz = x0 + z * (size_t)width * length;
			coord.z = z;

			for (int y = 0; y < length; y++)
			{
				const float *vdrefy = vdrefz + y * width;
				coord.y = y;

				for (int x = 0; x < width; x++)
				{
					const float *vdrefx = vdrefy + x;
					coord.x = x;
					Weight wt;
					ValCoef(origin, coord, coeffv[idx], &wt);
					if (wt.x_min >= 0 && wt.x_min < width && wt.y_min >= 0 &&
						wt.y_min < length)
					{									 //(x_min, y_min)
						n = wt.x_min + wt.y_min * width; // index in reproj

						ax[n] += (1 - wt.x_min_del) * (1 - wt.y_min_del) * (*vdrefx);
						w[n] += (1 - wt.x_min_del) * (1 - wt.y_min_del);
						// if(id==0){
						// file<<"n: "<<n<<"a: "<<a<<std::endl;
						// }
					}
					if ((wt.x_min + 1) >= 0 && (wt.x_min + 1) < width &&
						wt.y_min >= 0 && wt.y_min < length)
					{										 //(x_min+1, y_min)
						n = wt.x_min + 1 + wt.y_min * width; // index in reproj
						ax[n] += wt.x_min_del * (1 - wt.y_min_del) * (*vdrefx);
						w[n] += wt.x_min_del * (1 - wt.y_min_del);
						// if(id==0){
						// file<<"n: "<<n<<"a: "<<a<<std::endl;
						// }
					}
					if (wt.x_min >= 0 && wt.x_min < width && (wt.y_min + 1) >= 0 &&
						(wt.y_min + 1) < length)
					{										   //(x_min, y_min+1)
						n = wt.x_min + (wt.y_min + 1) * width; // index in reproj

						ax[n] += (1 - wt.x_min_del) * wt.y_min_del * (*vdrefx);
						w[n] += (1 - wt.x_min_del) * wt.y_min_del;
						// if(id==0){
						// file<<"n: "<<n<<"a: "<<a<<std::endl;
						// }
					}
					if ((wt.x_min + 1) >= 0 && (wt.x_min + 1) < width &&
						(wt.y_min + 1) >= 0 &&
						(wt.y_min + 1) < length)
					{												 //(x_min+1, y_min+1)
						n = (wt.x_min + 1) + (wt.y_min + 1) * width; // index in reproj
						ax[n] += wt.x_min_del * wt.y_min_del * (*vdrefx);
						w[n] += wt.x_min_del * wt.y_min_del;
						// if(id==0){
						// file<<"n: "<<n<<"a: "<<a<<std::endl;
						// }
					}
				}
			}
		}

		// file.close();
		for (int z = 0; z < height; z++)
		{
			// float *vdrefz = ax + z * (size_t) width *  length;
			float *atax_z = atax + z * (size_t)width * length;

			coord.z = z;

			for (int y = 0; y < length; y++)
			{
				// float *vdrefy = vdrefz + y *  width;
				float *atax_y = atax_z + y * width;

				coord.y = y;

				for (int x = 0; x < width; x++)
				{
					// float *vdrefx = vdrefy + x;
					float *atax_x = atax_y + x;

					coord.x = x;
					Weight wt;
					ValCoef(origin, coord, coeffv[idx], &wt);

					if (wt.x_min >= 0 && wt.x_min < width && wt.y_min >= 0 &&
						wt.y_min < length)
					{									 //(x_min, y_min)
						n = wt.x_min + wt.y_min * width; // index in reproj
						*atax_x += (1 - wt.x_min_del) * (1 - wt.y_min_del) * ax[n] / w[n];
					}
					if ((wt.x_min + 1) >= 0 && (wt.x_min + 1) < width &&
						wt.y_min >= 0 && wt.y_min < length)
					{										 //(x_min+1, y_min)
						n = wt.x_min + 1 + wt.y_min * width; // index in reproj
						*atax_x += wt.x_min_del * (1 - wt.y_min_del) * ax[n] / w[n];
					}
					if (wt.x_min >= 0 && wt.x_min < width && (wt.y_min + 1) >= 0 &&
						(wt.y_min + 1) < length)
					{										   //(x_min, y_min+1)
						n = wt.x_min + (wt.y_min + 1) * width; // index in reproj
						*atax_x += (1 - wt.x_min_del) * wt.y_min_del * ax[n] / w[n];
					}
					if ((wt.x_min + 1) >= 0 && (wt.x_min + 1) < width &&
						(wt.y_min + 1) >= 0 &&
						(wt.y_min + 1) < length)
					{												 //(x_min+1, y_min+1)
						n = (wt.x_min + 1) + (wt.y_min + 1) * width; // index in reproj
						*atax_x += wt.x_min_del * wt.y_min_del * ax[n] / w[n];
					}
				}
			}
		}
	}
	free(ax);
	free(w);
}

// ATb+mu*I^T(u1-d1,u2-d2,u3-d3)
// 结果在atb_lt
void ATbmuIT(float *atb, float *atb_lt, const float *uk, const float *dk, int width, int length, int height, float mu)
{
	size_t volsize = width * length * height;
	//  std::ofstream file;
	//  file.open("atb.txt",std::ios::out);
	for (int z = 0; z < height; z++)
	{
		float *drez_atb = atb + z * (size_t)width * length; // 体的值
		const float *drez_u1 = uk + z * (size_t)width * length;
		const float *drez_d1 = dk + z * (size_t)width * length;
		float *drez_atblt = atb_lt + z * (size_t)width * length; // 体的值

		for (int y = 0; y < length; y++)
		{
			float *drey_atb = drez_atb + y * width;
			const float *drey_u1 = drez_u1 + y * width;
			const float *drey_d1 = drez_d1 + y * width;
			float *drey_atblt = drez_atblt + y * width;

			for (int x = 0; x < width; x++)
			{
				float *drex_atb = drey_atb + x;
				const float *drex_u1 = drey_u1 + x;
				const float *drex_d1 = drey_d1 + x;
				float *drex_atblt = drey_atblt + x; // xn

				*drex_atblt = (*drex_u1 - *drex_d1) * mu + *drex_atb;
			}
		}
	}
}

// ATA(X)+mu*LTL(X,L,LT);
void ATAmuLTL(float *vol, float *ata, float *ltl, float mu,
			  int width, int length, int height)
{
	size_t volsize = width * length * height;
	// float *ata = (float*)malloc(sizeof(float)*volsize);//ata(x)
	// memset(ata,0,sizeof(float)*volsize);
	float *Lx = (float *)malloc(sizeof(float) * volsize); // 和体一样大
	float *Ly = (float *)malloc(sizeof(float) * volsize);
	float *Lz = (float *)malloc(sizeof(float) * volsize);
	size_t n;

	// ATA(vol,t,t_b,ata,num,width,length,height);

	for (int z = 0; z < height; z++)
	{
		float *drez = vol + z * (size_t)width * length; // 体的值
		float *drez_lx = Lx + z * (size_t)width * length;
		float *drez_ly = Ly + z * (size_t)width * length;
		float *drez_lz = Lz + z * (size_t)width * length;
		float *drez_l = ltl + z * (size_t)width * length;
		float *drez_ATA = ata + z * (size_t)width * length;

		for (int y = 0; y < length; y++)
		{
			float *drey = drez + y * width; // 值
			float *drey_lx = drez_lx + y * width;
			float *drey_ly = drez_ly + y * width;
			float *drey_lz = drez_lz + y * width;
			float *drey_l = drez_l + y * width;
			float *drey_ATA = drez_ATA + y * width;

			for (int x = 0; x < width; x++)
			{
				float *drex = drey + x;
				float *L_x = drey_lx + x;
				float *L_y = drey_ly + x;
				float *L_z = drey_lz + x;
				float *dreL = drey_l + x;
				float *dreATA = drey_ATA + x;

				// L
				if (x + 1 == width)
				{
					*L_x = 0;
				}
				else
				{
					*L_x = *drex - *(drex + 1);
				}
				if (y + 1 == length)
				{
					*L_y = 0;
				}
				else
				{
					*L_y = *drex - *(drex + width); // yn-yn+1
				}
				if (z + 1 == height)
				{
					*L_z = 0;
				}
				else
				{
					*L_z = *drex - *(drex + width * length); // zn-zn+1
				}
				// if (y == 1) {//*L_y的第二行
				//	printf("vol:%f   ",*L_y );
				// }
				if (x > 0)
				{
					*dreL += *L_x - *(L_x - 1); // xn-xn+1
				}
				else
				{
					*dreL += *L_x;
				}
				// L_y
				if (y > 0)
				{
					*dreL += *L_y - *(L_y - width);
				}
				else
				{
					*dreL += *L_y;
				}
				// L_z
				if (z > 0)
				{
					*dreL += *L_z - *(L_z - width * length);
				}
				else
				{
					*dreL += *L_z;
				}
				// test
				//  if (y == 7) {//q的第x+1列
				//  	printf("vol:%f   ", *dreL);
				//  }
				*dreL = *dreATA + mu * (*dreL);
			}
		}
	}
	// memcpy( data,L,volsize);
	// free(ata);
	free(Lx);
	free(Ly);
	free(Lz);
}

// ATA(X)+mu(x);
void ATAmuITI(const float *vol, float *ata, float *ltl, float mu,
			  int width, int length, int height)
{
	size_t volsize = width * length * height;
	size_t n;

	// ATA(vol,t,t_b,ata,num,width,length,height);

	for (int z = 0; z < height; z++)
	{
		float *drez_l = ltl + z * (size_t)width * length;
		const float *drez_vol = vol + z * (size_t)width * length;
		const float *drez_ATA = ata + z * (size_t)width * length;

		for (int y = 0; y < length; y++)
		{
			float *drey_l = drez_l + y * width;
			const float *drey_vol = drez_vol + y * width;
			const float *drey_ATA = drez_ATA + y * width;

			for (int x = 0; x < width; x++)
			{
				float *dreL = drey_l + x;
				const float *drex_vol = drey_vol + x;
				const float *dreATA = drey_ATA + x;
				*dreL = *dreATA + mu * (*drex_vol);
			}
		}
	}
}

void applycg(int width, int length, int height, float *voldata, float *atax0, float *ATb, int numberIteration, float mu, const Point3DF &origin, Coeff coeffv[], int proj, float *htb)
{
	// 初始化
	size_t volsize = length * width * height;
	float *r0, *p0;

	if ((r0 = (float *)malloc(sizeof(float) * volsize)) == NULL)
	{
		printf("false2");
	}
	if ((p0 = (float *)malloc(sizeof(float) * volsize)) == NULL)
	{
		printf("false3");
	}
	float *ax = (float *)malloc(sizeof(float) * volsize);

	// float *h0=(float *)malloc(sizeof(float)*volsize);//==A(x)
	// //float *r1=(float *)malloc(sizeof(float)*volsize);
	// float *r0=(float *)malloc(sizeof(float)*volsize);//     =  ATb- data   initial gradient
	// float *p0=(float *)malloc(sizeof(float)*volsize);//    =  r0 ;   initial search direction 初始搜索方向

	// float *x1=(float *)malloc(sizeof(float)*volsize);

	double d1; // r(k+1)T*r(k+1)
	double d0; // rT*r
	double beta = 0;

	for (int k = 0; k < volsize; k++)
	{
		r0[k] = atax0[k] - ATb[k];
		p0[k] = -r0[k];
	}
	// memcpy(p0,r0,sizeof(float)*volsize);
	float *h0;
	if ((h0 = (float *)malloc(sizeof(float) * volsize)) == NULL)
	{
		printf("false1");
	}

	for (int i = 0; i < numberIteration; i++)
	{ // 最大迭代次数
		// p0输入 h0输出=Apk
		printf("cg i:%d\n", i);
		memset(h0, 0, sizeof(float) * volsize); // h0是作为Apk 是一个中间变量
		{
			float *atax = (float *)malloc(sizeof(float) * volsize);
			memset(atax, 0, sizeof(float) * volsize);
			Reproject_admm_atax(origin, width, length, height, p0, coeffv, atax, proj); // p0输入 atax输出
			ATAmuITI(p0, atax, h0, mu, width, length, height);							// ATA=ATAmuLTLlambda
			free(atax);
		}
		// 计算alpha
		double alpha = 0; // rT*r/pk*A*p
		double d2 = 0;
		d1 = 0;
		d0 = 0;
		for (int k = 0; k < volsize; k++)
		{
			d0 += r0[k] * p0[k];
			d2 += h0[k] * p0[k];
		}
		// MPI_Allreduce(MPI_IN_PLACE, &d0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		// MPI_Allreduce(MPI_IN_PLACE, &d2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		if (d2 > 10e-6 || -d2 > 10e-6)
		{
			alpha = -d0 / d2;
		}
		else
		{
			alpha = 0;
		}

		// 计算残差x1、r1、d1
		for (int k = 0; k < volsize; k++)
		{
			voldata[k] = voldata[k] + alpha * p0[k];
			// if(  data[k]<0){
			// 	printf("k:%d\n",k);
			// }
		}
		{
			//   data输入 ax输出=Ax
			float *atax = (float *)malloc(sizeof(float) * volsize);
			memset(atax, 0, sizeof(float) * volsize);
			Reproject_admm_atax(origin, width, length, height, voldata, coeffv, atax, proj);
			ATAmuITI(voldata, atax, ax, mu, width, length, height);
			free(atax);
		}
		for (int k = 0; k < volsize; k++)
		{
			// r1[k]=r0[k]-al0pha * h0[k];//rk+1
			r0[k] = ax[k] - ATb[k]; // rk+1
			d1 += r0[k] * h0[k];	//
									// r02 += r0[k]*r0[k];
		}
		// MPI_Allreduce(MPI_IN_PLACE, &d1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		//  printf("r02:%f\n",r02);
		if (d2 > 10e-6 || -d2 > 10e-6)
		{
			beta = d1 / d2;
		}
		else
		{
			beta = 0;
		}
		// 计算p1 更新变量
		for (int k = 0; k < volsize; k++)
		{
			p0[k] = -r0[k] + beta * p0[k]; //			p1[k]=r1[k]+beta*p0[k]; p0[k]  =  p1[k];
		}

		// costfun(origin,vol,htb, coeffv,proj, projs);
	}
	free(h0);
	// free(x1);
	// free(r1);
	free(r0);
	free(p0);
	free(ax);
}

extern "C" void PnPDataFidelity(float *proj, float *voldata, Coeff *coeffv, int projsZ, int width, int length, int height,
								float lamb, //  lamb控制收敛速度
								const float *u_k, const float *d_k)
// 可以用稀疏矩阵存放A：数据：4*体的大小， 但是切分的话其实不用太考虑大小
{
	Point3DF origin;
	origin.x = width * .5;
	origin.y = length * .5;
	origin.z = height * .5;
	// printf("origin.x :%f\n", origin.x);
	int numberIteration = 5; // cg循环次数
	size_t pxsize = width * length;
	size_t volsize = length * width * height;

	float *htb = (float *)malloc(sizeof(float) * volsize); // 体
	memset(htb, 0, sizeof(float) * volsize);			   // 把体置为0

	// coeffv[idx] 将会是一个 Coeff 类型
	for (int idx = 0; idx < projsZ; idx++)
	{
		// printf("coeffv[idx].a[0]:%f\n", coeffv[idx].a[0]);
		//  printf("idx:%d\n",idx);
		//  printf("proj[idx * pxsize]:%f\n",proj[idx * pxsize]);
		Reproject_admm_htb(origin, width, length, height, coeffv[idx], htb, &(proj[idx * pxsize])); // 相当于backpro 计算结果是H^t*b
	}

	printf("**  PnP with Data Fidelity item iteration,iter **\n");
	float *x0 = (float *)malloc(sizeof(float) * volsize);
	float *htb_lt = (float *)malloc(sizeof(float) * volsize);
	memset(htb_lt, 0, sizeof(float) * volsize);
	memset(x0, 0, sizeof(float) * volsize);
	// ATb+muKT(u-d)
	// ATbmuIT(htb, htb_lt, u_k, d_k,  width,  length,  height, lamb);
	ATbmuIT(htb, htb_lt, u_k, d_k, width, length, height, lamb); // 在循环内 每次更新ud
	// printf("ATAmuLTL");
	// vol不变 x0ATA(vol)+mu*LTL(vol,L,LTL)
	{
		float *atax = (float *)malloc(sizeof(float) * volsize);
		memset(atax, 0, sizeof(float) * volsize);
		Reproject_admm_atax(origin, width, length, height, voldata, coeffv, atax, projsZ);
		ATAmuITI(voldata, atax, x0, lamb, width, length, height); // x0里面放的是ATA(X)+LTLX 是一个初始值
		free(atax);
	}
	// vol就是ATAx0
	printf("**  Conjugate Gradient Method  **\n");
	applycg(width, length, height, voldata, x0, htb_lt, numberIteration, lamb, origin, coeffv, projsZ, htb); // vol已经更新->x_k 也会作为下一次cg迭代的初始值
	free(x0);
	free(htb_lt);
	free(htb);
}

// int ATOM(options &opt, int myid, int procs)
// {
// MrcStackM projs, mrcvol;
// if (!projs.ReadFile(opt.input))
// {
// 	printf("File %s cannot access.\n", opt.input);

// 	return -1;
// }

// if (myid == 0)
// {
// 	projs.ReadHeader();
// }

// MPI_Bcast(&(projs.header), sizeof(MRCheader), MPI_CHAR, 0, MPI_COMM_WORLD);

// mrc  InitializeHeader();
// mrc  SetSize(projs.X(), projs.Y(), opt.thickness);

// std::vector<float> angles;
// ReadAngles(angles, opt.angle);

// std::vector<float> xangles;

// if (opt.xangle[0] != '\0')
// {
// 	ReadAngles(xangles, opt.xangle);
// }
// else
// {
// 	xangles.resize(angles.size(), 0.0);
// }

// std::vector<Coeff> params;
// TranslateAngleToCoefficients(angles, xangles, params);

// Geometry geo;
//  geo.offset = opt.offset;
//  geo.pitch_angle = opt.pitch_angle;
//  geo.zshift = opt.zshift;

// DecorateCoefficients(params, geo);
//}
void UpdateU(MrcStackM &projs, float *voldata, int width, int length, int height, float *u_k, float *d_k, float eta, float lamb)
{
	size_t volsize = width * length * height;
	// std::vector<float> v_vol;
	// v_vol.reserve(volsize); // 和原图一样大 为了放原图
	// v_vol.assign(voldata, voldata + volsize);
	//lwt3<float> lift3(v_vol, width, length, height, lft);

	// std::vector<float> aLLL;
	// aLLL.reserve(size); // 放alll
	// lift3.getLLH(aLLL);
	// std::vector<float> vectors[8];
	// for (int i = 0; i < 8; i++)
	// {
	// 	int size = l_length[0 + i * 3] * l_length[1 + i * 3] * l_length[2 + i * 3];
	// 	vectors[i].reserve(size);
	// }
	// lift3.getCoef(vectors[0], vectors[1], vectors[2], vectors[3], vectors[4], vectors[5], vectors[6], vectors[7]);



	//for (int i = 0; i < 8; i++)
	//{
		//int size = l_length[0 + i * 3] * l_length[1 + i * 3] * l_length[2 + i * 3];
		float *temp = (float *)malloc(sizeof(float) * volsize);
		for (int j = 0; j < volsize; j++)
		{
			temp[j] = voldata[j] + d_k[j]*1;
		}

		// 计算缩放因子
		float *max_it = std::max_element(temp, temp + volsize);
		float *min_it = std::min_element(temp, temp + volsize);

		MPI_Allreduce(MPI_IN_PLACE, &(*max_it), 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &(*min_it), 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
		float scale = 255.0 / (*max_it - *min_it);

		float offset = -(*min_it) * scale;
		// 计算逆向的 scale 和 offset
		float inv_scale = 1.0 / scale;
		float inv_offset = -offset / scale;

		for (int z = 0; z < height; z++) // 每一层
		{
			cv::Mat slice(width, length, CV_32F, &temp[z * width * length]);
			// 只接受0-255的输入
			slice.convertTo(slice, CV_8U, scale, offset);
			cv::Mat sliceFiltered(width, length, CV_8U);

			cv::xphoto::bm3dDenoising(slice, sliceFiltered, eta / lamb); // eta/lamb是滤波强度
			// 转换回去
			sliceFiltered.convertTo(sliceFiltered, CV_32F, inv_scale, inv_offset);

			memcpy(&u_k[ z *width * length], sliceFiltered.data, width * length * sizeof(float));
		}

		free(temp);
	//}
}
// void UpdateD(float *voldata, int width, int length, int height, const float *u_k, float *d_k, size_t volsize)
// {
// 	size_t volsize = width * length * height;
// 	// std::vector<float> v_vol;
// 	// v_vol.reserve(width * length * height);
// 	// v_vol.assign(voldata, voldata + width * length * height);
// 	//lwt3<float> lift3(v_vol, width, length, height, lft);
// 	// std::vector<float> vectors[8];
// 	// for (int i = 0; i < 8; i++)
// 	// {
// 	// 	int size = l_length[0 + i * 3] * l_length[1 + i * 3] * l_length[2 + i * 3];
// 	// 	vectors[i].reserve(size);
// 	// }
// 	// lift3.getCoef(vectors[0], vectors[1], vectors[2], vectors[3], vectors[4], vectors[5], vectors[6], vectors[7]);
// 	//int index = 0;
// 	// for (int i = 0; i < 8; i++)
// 	// {
// 		//int size = l_length[0 + i * 3] * l_length[1 + i * 3] * l_length[2 + i * 3];

// 		for (int j = 0; j < volsize; j++)
// 		{
// 			d_k[j] = d_k[j] + vectors[j] - u_k[j];
// 			//index++;
// 		}
// 	//}
// }
void PnPDataFidelity(MrcStackM &projs, float *voldata, Coeff *coeffv, int projsZ, int width, int length, int height,
					 float lamb, //  lamb控制收敛速度
					 float *u_k, float *d_k)
// 可以用稀疏矩阵存放A：数据：4*体的大小， 但是切分的话其实不用太考虑大小
{
	int eta = 1; //?
	Point3DF origin;
	origin.x = width * .5;
	origin.y = length * .5;
	origin.z = height * .5;
	// printf("origin.x :%f\n", origin.x);
	int numberIteration = 10; // cg循环次数
	int pxsize = width * length;
	size_t volsize = length * width * height;

	float *htb = (float *)malloc(sizeof(float) * volsize); // 体
	memset(htb, 0, sizeof(float) * volsize);			   // 把体置为0

	// coeffv[idx] 将会是一个 Coeff 类型
	for (int idx = 0; idx < projsZ; idx++)
	{
		Slice proj(width, length);
		projs.ReadSliceZ(idx, proj.data); // 全部读入到proj,data中
		// printf("coeffv[idx].a[0]:%f\n", coeffv[idx].a[0]);
		//  printf("idx:%d\n",idx);
		//  printf("proj[idx * pxsize]:%f\n",proj[idx * pxsize]);
		Reproject_admm_htb(origin, width, length, height, coeffv[idx], htb, proj.data); // 相当于backpro 计算结果是H^t*b
	}
	// std::vector<float> v_vol;
	// v_vol.reserve(volsize);
	// v_vol.assign(voldata, voldata + volsize);

	// string name = "haar";
	// int J = 1;
	// liftscheme blift(name);
	// std::vector<int> v_length;
	// lwt3<float> lift3(v_vol, width, length, height, blift, J);
	// lift3.getDim(v_length);

	printf("**  PnP with Data Fidelity item iteration,iter **\n");

		float *x0 = (float *)malloc(sizeof(float) * volsize);
		float *htb_lt = (float *)malloc(sizeof(float) * volsize);
		memset(htb_lt, 0, sizeof(float) * volsize);
		memset(x0, 0, sizeof(float) * volsize);
		// ATb+muKT(u-d)
		// ATbmuIT(htb, htb_lt, u_k, d_k,  width,  length,  height, lamb);
		ATbmuIT(htb, htb_lt, u_k, d_k, width, length, height, lamb); // 在循环内 每次更新ud
		// printf("ATAmuLTL");
		// vol不变 x0ATA(vol)+mu*LTL(vol,L,LTL)
		{
			float *atax = (float *)malloc(sizeof(float) * volsize);
			memset(atax, 0, sizeof(float) * volsize);
			Reproject_admm_atax(origin, width, length, height, voldata, coeffv, atax, projsZ);
			ATAmuITI(voldata, atax, x0, lamb, width, length, height); // x0里面放的是ATA(X)+LTLX 是一个初始值
			free(atax);
		}
		// vol就是ATAx0
		printf("**  Conjugate Gradient Method  **\n");
		applycg(width, length, height, voldata, x0, htb_lt, numberIteration, lamb, origin, coeffv, projsZ, htb); // vol已经更新->x_k 也会作为下一次cg迭代的初始值
																												 // 更新u 使用去噪算子NLM
	//	UpdateU(projs, voldata, width, length, height, u_k, d_k, eta, lamb, blift, v_length);
		// 更新d
	//	UpdateD(voldata, width, length, height, u_k, d_k, volsize, blift, v_length);
	
	free(x0);
	free(htb_lt);
	free(htb);
}
// void PnPDataFidelity_BM3D(MrcStackM &projs, float *voldata, Coeff *coeffv, int projsZ, int width, int length, int height,
// 					 float lamb, //  lamb控制收敛速度
// 					 float *u_k, float *d_k)
// // 可以用稀疏矩阵存放A：数据：4*体的大小， 但是切分的话其实不用太考虑大小
// {
// 	int eta = 1; //?
// 	Point3DF origin;
// 	origin.x = width * .5;
// 	origin.y = length * .5;
// 	origin.z = height * .5;
// 	// printf("origin.x :%f\n", origin.x);
// 	int numberIteration = 2; // cg循环次数
// 	int pxsize = width * length;
// 	size_t volsize = length * width * height;

// 	float *htb = (float *)malloc(sizeof(float) * volsize); // 体
// 	memset(htb, 0, sizeof(float) * volsize);			   // 把体置为0

// 	// coeffv[idx] 将会是一个 Coeff 类型
// 	for (int idx = 0; idx < projsZ; idx++)
// 	{
// 		Slice proj(width, length);
// 		projs.ReadSliceZ(idx, proj.data); // 全部读入到proj,data中
// 		// printf("coeffv[idx].a[0]:%f\n", coeffv[idx].a[0]);
// 		//  printf("idx:%d\n",idx);
// 		//  printf("proj[idx * pxsize]:%f\n",proj[idx * pxsize]);
// 		Reproject_admm_htb(origin, width, length, height, coeffv[idx], htb, proj.data); // 相当于backpro 计算结果是H^t*b
// 	}
// 	// std::vector<float> v_vol;
// 	// v_vol.reserve(volsize);
// 	// v_vol.assign(voldata, voldata + volsize);

// 	// string name = "haar";
// 	// int J = 1;
// 	// liftscheme blift(name);
// 	// std::vector<int> v_length;
// 	// lwt3<float> lift3(v_vol, width, length, height, blift, J);
// 	// lift3.getDim(v_length);
// 	for (int i = 0; i < maxOutIter; i++)
// 	{

// 		printf("**  PnP with Conjugate Gradient method ,iter:%d  **\n", i);
// 		float *x0 = (float *)malloc(sizeof(float) * volsize);
// 		float *htb_lt = (float *)malloc(sizeof(float) * volsize);
// 		memset(htb_lt, 0, sizeof(float) * volsize);
// 		memset(x0, 0, sizeof(float) * volsize);
// 		// ATb+muKT(u-d)
// 		// ATbmuIT(htb, htb_lt, u_k, d_k,  width,  length,  height, lamb);
// 		ATbmuIT(htb, htb_lt, u_k, d_k, width, length, height, lamb); // 在循环内 每次更新ud
// 		// printf("ATAmuLTL");
// 		// vol不变 x0ATA(vol)+mu*LTL(vol,L,LTL)
// 		{
// 			float *atax = (float *)malloc(sizeof(float) * volsize);
// 			memset(atax, 0, sizeof(float) * volsize);
// 			Reproject_admm_atax(origin, width, length, height, voldata, coeffv, atax, projsZ);
// 			ATAmuITI(voldata, atax, x0, lamb, width, length, height); // x0里面放的是ATA(X)+LTLX 是一个初始值
// 			free(atax);
// 		}
// 		// vol就是ATAx0
// 		printf("**  Conjugate Gradient Method  **\n");
// 		applycg(width, length, height, voldata, x0, htb_lt, numberIteration, lamb, origin, coeffv, projsZ, htb); // vol已经更新->x_k 也会作为下一次cg迭代的初始值
// 																												 // 更新u 使用去噪算子NLM
// 		UpdateU(projs, voldata, width, length, height, u_k, d_k, eta);
// 		// 更新d
// 		UpdateD(voldata, width, length, height, u_k, d_k, volsize);
// 	}
// 	free(x0);
// 	free(htb_lt);
// 	free(htb);
// }
void PnP(const Point3DF &origin, MrcStackM &projs, Volume &vol, Coeff coeffv[],
		 int maxOutIter, // 最大外循环迭代次数
		 float eta)		 // 松弛变量
// 可以用稀疏矩阵存放A：数据：4*体的大小， 但是切分的话其实不用太考虑大小
{
	float lamb = 0.4;		 // 和收敛速度有关
	int numberIteration = 2; // cg循环次数
	int pxsize = projs.X() * projs.Y();
	size_t volsize = vol.length * vol.width * vol.height;

	float *u_k = (float *)malloc(sizeof(float) * volsize); // 如果不够
	float *d_k = (float *)malloc(sizeof(float) * volsize);
	memset(u_k, 0, sizeof(float) * volsize);
	memset(d_k, 0, sizeof(float) * volsize);

	PnPDataFidelity(projs, vol.data, coeffv, projs.Z(), vol.width, vol.length, vol.height, lamb, u_k, d_k);
	// 更新u、d
	// 更新u 使用去噪算子NLM
	// UpdateU(projs, vol, u_k, d_k, eta, lamb, blift, v_length);
	// 更新d
	// UpdateD(vol, u_k, d_k, volsize, blift, v_length);

	free(u_k);
	free(d_k);
}
struct SysInfo
{
	int id;
	int procs;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int namelen;
};

extern "C" void ATOM(float *voldata, char *optinput, char *optangle, char *output, int thickness, float lamb, float *u_k, float *d_k)
{

	SysInfo info;

	int err;
	// info.id = rank;
	// info.procs = size;
	// printf("rank:%d\n", rank);

	err = MPI_Init(NULL, NULL); // 一般不会有影响
	if (err != MPI_SUCCESS)
	{
		fprintf(stderr, "MPI_Init failed with error code %d\n", err);
		exit(EXIT_FAILURE);
	}
	err = MPI_Comm_rank(MPI_COMM_WORLD, &(info.id));
	if (err != MPI_SUCCESS)
	{
		fprintf(stderr, "MPI_Comm_rank failed with error code %d\n", err);
		exit(EXIT_FAILURE);
	}
	// info.id=rank;
	MPI_Comm_size(MPI_COMM_WORLD, &(info.procs));
	MPI_Get_processor_name(info.processor_name, &(info.namelen));

	MrcStackM projs, mrcvol;

	if (!projs.ReadFile(optinput))
	{
		printf("File %s cannot access.\n", optinput);
	}

	if (info.id == 0)
	{
		projs.ReadHeader();
	}

	MPI_Bcast(&(projs.header), sizeof(MRCheader), MPI_CHAR, 0, MPI_COMM_WORLD);

	mrcvol.InitializeHeader();
	mrcvol.SetSize(projs.X(), projs.Y(), thickness);

	std::vector<float> angles;
	ReadAngles(angles, optangle);

	std::vector<Coeff> params;
	params.resize(angles.size());
	TranslateAngleToCoefficients(angles.data(), params.data(), angles.size());

	DecorateCoefficients(params.data(), 0, 0, 0, angles.size());

	int zrem = mrcvol.Z() % info.procs;		  // 得到的余数
	int volz;								  // the start slice of reproject per process
	int baseHeight = mrcvol.Z() / info.procs; // 每个线程本应该得到的数
	int remainder = 0;
	int height = baseHeight;
	if (zrem != 0) // 如果有剩下的
	{
		if (baseHeight % 2 != 0)
		{
			baseHeight = baseHeight - 1;
			remainder = mrcvol.Z() - baseHeight * info.procs; // 现在剩下的余数
		}
		else
		{
			baseHeight = baseHeight;
			remainder = zrem;
		}

		if (remainder % 2 == 0) // 如果剩余层数是偶数，可以平均分配给一些处理器
		{
			if (info.id < remainder / 2) // 给开始的多分两张
			{
				height = baseHeight + 2;
			}
			else
			{
				height = baseHeight;
			}
		}
		else // 如果剩余层数是奇数，需要特殊处理
		{

			if (info.id < remainder / 2)
			{
				height = baseHeight + 2;
			}
			else if (info.id == remainder / 2)
			{
				height = baseHeight + 1;
			}
			else
			{
				height = baseHeight;
			}
		}
	}

	// 计算volz
	if (info.id == 0) // 第一个线程的起始位置总是 0
	{
		volz = 0;
	}
	else
	{
		// 其他线程的起始位置是前面所有线程的层数之和
		volz = baseHeight * info.id;

		// 如果余数是偶数
		if (remainder % 2 == 0)
		{
			if (info.id > remainder / 2)
			{
				volz += remainder;
			}
			else
			{
				volz += 2 * info.id;
			}
		}
		else // 如果余数是奇数
		{

			if (info.id > remainder / 2)
			{
				volz += remainder;
			}
			else
			{
				volz += 2 * info.id;
			}
		}
	}

	// Volume vol(0, 0, volz, mrcvol.Y(), mrcvol.X(), height);
	std::cout << info.id << std::endl;
	// 		  << "&(" << vol.width << "," << vol.length << "," << vol.height
	// 		  << ")" << std::endl;

	Point3DF origin;

	// AssignValue 原来的会丢失小数
	origin.x = mrcvol.X() * .5;
	origin.y = mrcvol.Y() * .5;
	origin.z = mrcvol.Z() * .5;

	if (info.id == 0)
	{
		printf("origin.x is %f, origin.y is %f, origin.z is %f\n", origin.x,
			   origin.y, origin.z);
	}

	PnPDataFidelity(projs, voldata, params.data(), projs.Z(), mrcvol.X(), mrcvol.Y(), height,
					lamb, //  lamb控制收敛速度
					u_k, d_k);
	// vol.data = voldata; // 都会改变
	//  PnP(origin, projs, vol, &params[0], 1, gamma);
	//   vol.height=vol.height-1;

	mrcvol.WriteToFile(output);

	if (info.id == 0)
	{
		mrcvol.WriteHeader();
	}

	{
		// if (myid == 0)
		// { // 第一个切片
		// 	// 起点不变volz=volz;
		// 	vol.height -= 1; // 后面多切一片
		// }
		// else if (myid == procs - 1) // 最后一个切片
		// {
		// 	vol.z = vol.z + 1; // 前面多切一片
		// 	vol.height -= 1;
		// }
		// else
		// {
		// 	// vol.z=vol.z+1;
		// 	vol.height -= 1;
		// }
		// std::cout << info.id << ": (" << vol.x << "," << vol.y << "," << vol.z << ")"
		// 		  << "&(" << vol.width << "," << vol.length << "," << vol.height
		// 		  << ")" << std::endl;

		// 分片写入
		int j = 0;
		for (j = 0; volz + j + 5 < volz + height; j += 5)
		{
			mrcvol.WriteBlock(volz + j, volz + j + 5, 'z', (voldata + (size_t)mrcvol.X() * mrcvol.Y() * j));
		}
		mrcvol.WriteBlock(volz + j, volz + height, 'z', (voldata + (size_t)mrcvol.X() * mrcvol.Y() * j));
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (info.id == 0)
	{
		mrcvol.UpdateHeader();
	}

	projs.Close();
	mrcvol.Close();

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize(); // parallel finish
}
