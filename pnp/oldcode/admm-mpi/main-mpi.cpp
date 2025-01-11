
#include "mrcmx/mrcstack.h"
#include "opts.h"
#include <fstream>
#include <iostream>
#include <vector>
// #include <eigen3/Eigen/Dense>
// #include <opencv4/opencv2/imgproc/imgproc_c.h>
// #include <opencv4/opencv2/imgproc/types_c.h>
// #include <opencv4/opencv2/highgui/highgui_c.h>
// #include <opencv4/opencv2/core/core_c.h>
// #include <opencv4/opencv2/photo/photo.hpp>
#include <opencv2/opencv.hpp>
// using namespace Eigen;

#define MILLION 1000000

#define PI_180 0.01745329252f

#ifndef PI
#define PI 3.14159265358979323846
#endif

#define D2R(__ANGLE__) ((__ANGLE__) * PI_180)

struct Coeff
{ // 20个双精度数，前10个是a后10个是b
	union
	{
		double p[20];
		struct
		{
			double a[10];
			double b[10];
		};
	};
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

void TranslateAngleToCoefficients(const std::vector<float> &angles,
								  const std::vector<float> &xangles,
								  std::vector<Coeff> &coeffs)
{
	coeffs.resize(angles.size());
	for (int i = 0; i < angles.size(); i++)
	{
		memset(coeffs[i].p, 0, sizeof(double) * 20);
		float beta = D2R(angles[i]);
		float alpha = D2R(xangles[i]);

		coeffs[i].a[0] = 0;						  //
		coeffs[i].a[1] = cos(beta);				  // x
		coeffs[i].a[2] = sin(alpha) * sin(beta);  // y
		coeffs[i].a[3] = -cos(alpha) * sin(beta); // z
		coeffs[i].b[0] = 0;						  //
		coeffs[i].b[1] = 0;						  // x
		coeffs[i].b[2] = cos(alpha);			  // y
		coeffs[i].b[3] = sin(alpha);			  // z
	}
}

/*solve inverse transfroms defined by Geometry; substitute the inversion into
 * coefficients*/
void DecorateCoefficients(std::vector<Coeff> &coeffs, const Geometry &geo)
{
	double alpha = -D2R(geo.pitch_angle), beta = D2R(geo.offset), t = -geo.zshift;
	double ca = cos(alpha), sa = sin(alpha), cb = cos(beta), sb = sin(beta);
	double ca2 = ca * ca, sa2 = sa * sa, cb2 = cb * cb, sb2 = sb * sb;

	for (int i = 0; i < coeffs.size(); i++)
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
	for (int i = 0; i < coeffs.size(); i++)
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

void Reproject(const Point3DF &origin, const Volume &vol, const Coeff &coeff,
			   Slice &reproj_val, Slice &reproj_wt)
{

	Point3D coord;
	int n;

	for (int z = 0; z < vol.height; z++)
	{
		float *vdrefz = vol.data + z * (size_t)vol.width * vol.length;
		coord.z = z + vol.z;

		for (int y = 0; y < vol.length; y++)
		{
			float *vdrefy = vdrefz + y * vol.width;
			coord.y = y + vol.y;

			for (int x = 0; x < vol.width; x++)
			{
				float *vdrefx = vdrefy + x;
				coord.x = x + vol.x;
				Weight wt;
				ValCoef(origin, coord, coeff, &wt);

				if (wt.x_min >= 0 && wt.x_min < vol.width && wt.y_min >= 0 &&
					wt.y_min < vol.length)
				{										 //(x_min, y_min)
					n = wt.x_min + wt.y_min * vol.width; // index in reproj 索引
					reproj_val.data[n] +=
						(1 - wt.x_min_del) * (1 - wt.y_min_del) * (*vdrefx);
					reproj_wt.data[n] += (1 - wt.x_min_del) * (1 - wt.y_min_del);
				}
				if ((wt.x_min + 1) >= 0 && (wt.x_min + 1) < vol.width &&
					wt.y_min >= 0 && wt.y_min < vol.length)
				{											 //(x_min+1, y_min)
					n = wt.x_min + 1 + wt.y_min * vol.width; // index in reproj
					reproj_val.data[n] += wt.x_min_del * (1 - wt.y_min_del) * (*vdrefx);
					reproj_wt.data[n] += wt.x_min_del * (1 - wt.y_min_del);
				}
				if (wt.x_min >= 0 && wt.x_min < vol.width && (wt.y_min + 1) >= 0 &&
					(wt.y_min + 1) < vol.length)
				{											   //(x_min, y_min+1)
					n = wt.x_min + (wt.y_min + 1) * vol.width; // index in reproj
					reproj_val.data[n] += (1 - wt.x_min_del) * wt.y_min_del * (*vdrefx);
					reproj_wt.data[n] += (1 - wt.x_min_del) * wt.y_min_del;
				}
				if ((wt.x_min + 1) >= 0 && (wt.x_min + 1) < vol.width &&
					(wt.y_min + 1) >= 0 &&
					(wt.y_min + 1) < vol.length)
				{													 //(x_min+1, y_min+1)
					n = (wt.x_min + 1) + (wt.y_min + 1) * vol.width; // index in reproj
					reproj_val.data[n] += wt.x_min_del * wt.y_min_del * (*vdrefx);
					reproj_wt.data[n] += wt.x_min_del * wt.y_min_del;
				}
			}
		}
	}
}
// 这个不是repro 是backpro
void Reproject_admm_htb(const Point3DF &origin, const Volume &vol, const Coeff &coeff, float *htb, const Slice &slc)
{
	Point3D coord;
	int n;
	size_t volsize = vol.width * vol.length * vol.height;
	for (int z = 0; z < vol.height; z++)
	{
		float *vdrefz = vol.data + z * (size_t)vol.width * vol.length;
		float *htb_z = htb + z * (size_t)vol.width * vol.length;
		coord.z = z + vol.z;

		for (int y = 0; y < vol.length; y++)
		{
			float *vdrefy = vdrefz + y * vol.width;
			float *htb_y = htb_z + y * vol.width;
			coord.y = y + vol.y;

			for (int x = 0; x < vol.width; x++)
			{
				float *vdrefx = vdrefy + x;
				float *htb_x = htb_y + x;
				coord.x = x + vol.x;
				Weight wt;
				ValCoef(origin, coord, coeff, &wt);
				if (wt.x_min >= 0 && wt.x_min < slc.width && wt.y_min >= 0 &&
					wt.y_min < slc.height)
				{																	 //(x_min, y_min)
					n = wt.x_min + wt.y_min * vol.width;							 // index in reproj
					*htb_x += (1 - wt.x_min_del) * (1 - wt.y_min_del) * slc.data[n]; // htb
				}
				if ((wt.x_min + 1) >= 0 && (wt.x_min + 1) < slc.width &&
					wt.y_min >= 0 && wt.y_min < slc.height)
				{											 //(x_min+1, y_min)
					n = wt.x_min + 1 + wt.y_min * vol.width; // index in reproj
					*htb_x += wt.x_min_del * (1 - wt.y_min_del) * slc.data[n];
				}
				if (wt.x_min >= 0 && wt.x_min < slc.width && (wt.y_min + 1) >= 0 &&
					(wt.y_min + 1) < slc.height)
				{											   //(x_min, y_min+1)
					n = wt.x_min + (wt.y_min + 1) * vol.width; // index in reproj
					*htb_x += (1 - wt.x_min_del) * wt.y_min_del * slc.data[n];
				}
				if ((wt.x_min + 1) >= 0 && (wt.x_min + 1) < slc.width &&
					(wt.y_min + 1) >= 0 &&
					(wt.y_min + 1) < slc.height)
				{													 //(x_min+1, y_min+1)
					n = (wt.x_min + 1) + (wt.y_min + 1) * vol.width; // index in reproj
					*htb_x += wt.x_min_del * wt.y_min_del * slc.data[n];
				}
			}
		}
	}
}
// repro+backpro
void Reproject_admm_atax(const Point3DF &origin, const Volume &vol, float *x0, Coeff coeffv[],
						 float *atax, int proj)

{
	Point3D coord;
	int n;
	size_t volsize = vol.width * vol.length * vol.height;
	float *ax = (float *)malloc(sizeof(float) * vol.width * vol.length);
	float *w = (float *)malloc(sizeof(float) * vol.width * vol.length);
	int pxsize = vol.width * vol.length;
	float a = 0;
	// std::ofstream file;
	// file.open("A.txt",std::ios::out);
	// int id;
	// MPI_Comm_rank(MPI_COMM_WORLD, &id);
	for (int idx = 0; idx < proj; idx++)
	{
		memset(ax, 0, sizeof(float) * vol.width * vol.length);
		memset(w, 0, sizeof(float) * vol.width * vol.length);
		for (int z = 0; z < vol.height; z++)
		{

			float *vdrefz = x0 + z * (size_t)vol.width * vol.length;
			coord.z = z + vol.z;

			for (int y = 0; y < vol.length; y++)
			{
				float *vdrefy = vdrefz + y * vol.width;
				coord.y = y + vol.y;

				for (int x = 0; x < vol.width; x++)
				{
					float *vdrefx = vdrefy + x;
					coord.x = x + vol.x;
					Weight wt;
					ValCoef(origin, coord, coeffv[idx], &wt);
					if (wt.x_min >= 0 && wt.x_min < vol.width && wt.y_min >= 0 &&
						wt.y_min < vol.length)
					{										 //(x_min, y_min)
						n = wt.x_min + wt.y_min * vol.width; // index in reproj

						ax[n] += (1 - wt.x_min_del) * (1 - wt.y_min_del) * (*vdrefx);
						w[n] += (1 - wt.x_min_del) * (1 - wt.y_min_del);
						// if(id==0){
						// file<<"n: "<<n<<"a: "<<a<<std::endl;
						// }
					}
					if ((wt.x_min + 1) >= 0 && (wt.x_min + 1) < vol.width &&
						wt.y_min >= 0 && wt.y_min < vol.length)
					{											 //(x_min+1, y_min)
						n = wt.x_min + 1 + wt.y_min * vol.width; // index in reproj
						ax[n] += wt.x_min_del * (1 - wt.y_min_del) * (*vdrefx);
						w[n] += wt.x_min_del * (1 - wt.y_min_del);
						// if(id==0){
						// file<<"n: "<<n<<"a: "<<a<<std::endl;
						// }
					}
					if (wt.x_min >= 0 && wt.x_min < vol.width && (wt.y_min + 1) >= 0 &&
						(wt.y_min + 1) < vol.length)
					{											   //(x_min, y_min+1)
						n = wt.x_min + (wt.y_min + 1) * vol.width; // index in reproj

						ax[n] += (1 - wt.x_min_del) * wt.y_min_del * (*vdrefx);
						w[n] += (1 - wt.x_min_del) * wt.y_min_del;
						// if(id==0){
						// file<<"n: "<<n<<"a: "<<a<<std::endl;
						// }
					}
					if ((wt.x_min + 1) >= 0 && (wt.x_min + 1) < vol.width &&
						(wt.y_min + 1) >= 0 &&
						(wt.y_min + 1) < vol.length)
					{													 //(x_min+1, y_min+1)
						n = (wt.x_min + 1) + (wt.y_min + 1) * vol.width; // index in reproj
						ax[n] += wt.x_min_del * wt.y_min_del * (*vdrefx);
						w[n] += wt.x_min_del * wt.y_min_del;
						// if(id==0){
						// file<<"n: "<<n<<"a: "<<a<<std::endl;
						// }
					}
				}
			}
		}
		MPI_Allreduce(MPI_IN_PLACE, ax, pxsize, MPI_FLOAT, MPI_SUM,
					  MPI_COMM_WORLD);
		// MPI_IN_PLACE:说明当前进程既发送又接受数据，而且要发送的数据和在要接收的数据的保存在同一内存。
		MPI_Allreduce(MPI_IN_PLACE, w, pxsize, MPI_FLOAT, MPI_SUM,
					  MPI_COMM_WORLD);

		// file.close();
		for (int z = 0; z < vol.height; z++)
		{
			// float *vdrefz = ax + z * (size_t)vol.width * vol.length;
			float *atax_z = atax + z * (size_t)vol.width * vol.length;

			coord.z = z + vol.z;

			for (int y = 0; y < vol.length; y++)
			{
				// float *vdrefy = vdrefz + y * vol.width;
				float *atax_y = atax_z + y * vol.width;

				coord.y = y + vol.y;

				for (int x = 0; x < vol.width; x++)
				{
					// float *vdrefx = vdrefy + x;
					float *atax_x = atax_y + x;

					coord.x = x + vol.x;
					Weight wt;
					ValCoef(origin, coord, coeffv[idx], &wt);

					if (wt.x_min >= 0 && wt.x_min < vol.width && wt.y_min >= 0 &&
						wt.y_min < vol.length)
					{										 //(x_min, y_min)
						n = wt.x_min + wt.y_min * vol.width; // index in reproj
						*atax_x += (1 - wt.x_min_del) * (1 - wt.y_min_del) * ax[n] / w[n];
					}
					if ((wt.x_min + 1) >= 0 && (wt.x_min + 1) < vol.width &&
						wt.y_min >= 0 && wt.y_min < vol.length)
					{											 //(x_min+1, y_min)
						n = wt.x_min + 1 + wt.y_min * vol.width; // index in reproj
						*atax_x += wt.x_min_del * (1 - wt.y_min_del) * ax[n] / w[n];
					}
					if (wt.x_min >= 0 && wt.x_min < vol.width && (wt.y_min + 1) >= 0 &&
						(wt.y_min + 1) < vol.length)
					{											   //(x_min, y_min+1)
						n = wt.x_min + (wt.y_min + 1) * vol.width; // index in reproj
						*atax_x += (1 - wt.x_min_del) * wt.y_min_del * ax[n] / w[n];
					}
					if ((wt.x_min + 1) >= 0 && (wt.x_min + 1) < vol.width &&
						(wt.y_min + 1) >= 0 &&
						(wt.y_min + 1) < vol.length)
					{													 //(x_min+1, y_min+1)
						n = (wt.x_min + 1) + (wt.y_min + 1) * vol.width; // index in reproj
						*atax_x += wt.x_min_del * wt.y_min_del * ax[n] / w[n];
					}
				}
			}
		}
	}
	free(ax);
	free(w);
}
inline void BilinearValue(const Slice &slc, const Weight &wt, float *val,
						  float *vwt)
{
	int n;
	if (wt.x_min >= 0 && wt.x_min < slc.width && wt.y_min >= 0 &&
		wt.y_min < slc.height)
	{ //(x_min, y_min)
		n = wt.x_min + wt.y_min * slc.width;
		*val += (1 - wt.x_min_del) * (1 - wt.y_min_del) * slc.data[n];
		*vwt += (1 - wt.x_min_del) * (1 - wt.y_min_del);
	}
	if ((wt.x_min + 1) >= 0 && (wt.x_min + 1) < slc.width && wt.y_min >= 0 &&
		wt.y_min < slc.height)
	{ //(x_min+1, y_min)
		n = wt.x_min + 1 + wt.y_min * slc.width;
		*val += wt.x_min_del * (1 - wt.y_min_del) * slc.data[n];
		*vwt += wt.x_min_del * (1 - wt.y_min_del);
	}
	if (wt.x_min >= 0 && wt.x_min < slc.width && (wt.y_min + 1) >= 0 &&
		(wt.y_min + 1) < slc.height)
	{ //(x_min, y_min+1)
		n = wt.x_min + (wt.y_min + 1) * slc.width;
		*val += (1 - wt.x_min_del) * wt.y_min_del * slc.data[n];
		*vwt += (1 - wt.x_min_del) * wt.y_min_del;
	}
	if ((wt.x_min + 1) >= 0 && (wt.x_min + 1) < slc.width &&
		(wt.y_min + 1) >= 0 && (wt.y_min + 1) < slc.height)
	{ //(x_min+1, y_min+1)
		n = wt.x_min + 1 + (wt.y_min + 1) * slc.width;
		*val += wt.x_min_del * wt.y_min_del * slc.data[n];
		*vwt += wt.x_min_del * wt.y_min_del;
	}
}
// t改变
inline void BilinearValue_admm(const Slice &slc, const Weight &wt, float *val,
							   TupNode *&t, int *num)
{
	int n;
	if (wt.x_min >= 0 && wt.x_min < slc.width && wt.y_min >= 0 &&
		wt.y_min < slc.height)
	{ //(x_min, y_min)
		n = wt.x_min + wt.y_min * slc.width;
		t[*num].c = n; /*列号*/
		*val += (1 - wt.x_min_del) * (1 - wt.y_min_del) * slc.data[n];
		t[*num].d += (1 - wt.x_min_del) * (1 - wt.y_min_del); //*vwt
		(*num)++;
	}
	if ((wt.x_min + 1) >= 0 && (wt.x_min + 1) < slc.width && wt.y_min >= 0 &&
		wt.y_min < slc.height)
	{ //(x_min+1, y_min)
		n = wt.x_min + 1 + wt.y_min * slc.width;
		t[*num].c = n;
		*val += wt.x_min_del * (1 - wt.y_min_del) * slc.data[n];
		t[*num].d += wt.x_min_del * (1 - wt.y_min_del);
		(*num)++;
	}
	if (wt.x_min >= 0 && wt.x_min < slc.width && (wt.y_min + 1) >= 0 &&
		(wt.y_min + 1) < slc.height)
	{ //(x_min, y_min+1)
		n = wt.x_min + (wt.y_min + 1) * slc.width;
		t[*num].c = n;
		*val += (1 - wt.x_min_del) * wt.y_min_del * slc.data[n];
		t[*num].d += (1 - wt.x_min_del) * wt.y_min_del;
		(*num)++;
	}
	if ((wt.x_min + 1) >= 0 && (wt.x_min + 1) < slc.width &&
		(wt.y_min + 1) >= 0 && (wt.y_min + 1) < slc.height)
	{ //(x_min+1, y_min+1)
		n = wt.x_min + 1 + (wt.y_min + 1) * slc.width;
		t[*num].c = n;
		*val += wt.x_min_del * wt.y_min_del * slc.data[n];
		t[*num].d += wt.x_min_del * wt.y_min_del;
		(*num)++;
	}
}

void BackProject(const Point3DF &origin, MrcStackM &projs, Volume &vol,
				 Coeff coeffv[])
{
	Slice proj(projs.X(), projs.Y());
	Point3D coord;

	memset(vol.data, 0, sizeof(float) * vol.length * vol.width * vol.height);

	for (int idx = 0; idx < projs.Z(); idx++)
	{
		printf("BPT begin to read %d projection for %d z-coordinate\n", idx, vol.z);
		projs.ReadSliceZ(idx, proj.data);

		for (int z = 0; z < vol.height; z++)
		{
			float *vdrefz = vol.data + z * (size_t)vol.width * vol.length;
			coord.z = z + vol.z;

			for (int y = 0; y < vol.length; y++)
			{
				float *vdrefy = vdrefz + y * vol.width;
				coord.y = y + vol.y;

				for (int x = 0; x < vol.width; x++)
				{
					coord.x = x + vol.x;
					Weight wt;
					float s = 0, c = 0;

					ValCoef(origin, coord, coeffv[idx], &wt);
					BilinearValue(proj, wt, &s, &c); // vol应该是HTB

					if (c)
					{
						*(vdrefy + x) += (float)(s / c);
					}
				}
			}
		}
	}
}
// void BackProject_admm(const Point3DF &origin, const Slice &proj, Volume &vol,
// 				 const Coeff &coeff)
// {
// 	//Slice proj(projs.X(), projs.Y());
// 	Point3D coord;

// 	//memset(vol.data, 0, sizeof(float) * vol.length * vol.width * vol.height);

// 	// for (int idx = 0; idx < projs.Z(); idx++)
// 	// {
// 		//printf("BPT begin to read %d projection for %d z-coordinate\n", idx, vol.z);
// 		//projs.ReadSliceZ(idx, proj.data);

// 		for (int z = 0; z < vol.height; z++)
// 		{
// 			float *vdrefz = vol.data + z * (size_t)vol.width * vol.length;
// 			coord.z = z + vol.z;

// 			for (int y = 0; y < vol.length; y++)
// 			{
// 				float *vdrefy = vdrefz + y * vol.width;
// 				coord.y = y + vol.y;

// 				for (int x = 0; x < vol.width; x++)
// 				{
// 					coord.x = x + vol.x;
// 					Weight wt;
// 					float s = 0, c = 0;

// 					ValCoef(origin, coord, coeff, &wt);
// 					BilinearValue(proj, wt, &s, &c);//vol应该是HTB

// 					if (c)
// 					{
// 						*(vdrefy + x) += (float)(s / c);
// 					}
// 				}
// 			}
// 		//}
// 	}
// }
void UpdateVolumeByProjDiff(const Point3DF &origin, const Slice &diff,
							Volume &vol, float gamma, const Coeff &coeff)
{

	Point3D coord;

	for (int z = 0; z < vol.height; z++)
	{
		float *vdrefz = vol.data + z * (size_t)vol.width * vol.length;
		coord.z = z + vol.z;

		for (int y = 0; y < vol.length; y++)
		{
			float *vdrefy = vdrefz + y * vol.width;
			coord.y = y + vol.y;

			for (int x = 0; x < vol.width; x++)
			{
				coord.x = x + vol.x;
				Weight wt;
				float s = 0, c = 0;

				ValCoef(origin, coord, coeff, &wt);
				BilinearValue(diff, wt, &s, &c);

				if (c)
				{
					*(vdrefy + x) += (float)(s / c) * gamma;
				}
			}
		}
	}
}

void UpdateWeightsByProjDiff(const Point3DF &origin, const Slice &diff,
							 Volume &values, Volume &weights,
							 const Coeff &coeff)
{

	Point3D coord;

	for (int z = 0; z < values.height; z++)
	{
		float *vdrefz = values.data + z * (size_t)values.width * values.length;
		float *wdrefz = weights.data + z * (size_t)weights.width * weights.length;
		coord.z = z + values.z;

		for (int y = 0; y < values.length; y++)
		{
			float *vdrefy = vdrefz + y * values.width;
			float *wdrefy = wdrefz + y * weights.width;
			coord.y = y + values.y;

			for (int x = 0; x < values.width; x++)
			{
				coord.x = x + values.x;
				Weight wt;
				float s = 0, c = 0;

				ValCoef(origin, coord, coeff, &wt);
				BilinearValue(diff, wt, &s, &c);

				*(vdrefy + x) += s;
				*(wdrefy + x) += c;
			}
		}
	}
}

void UpdateVolumeByWeights(Volume &vol, Volume &values, Volume &weights,
						   float gamma)
{
	size_t pxsize = vol.height * vol.width * vol.length;

	for (size_t i = pxsize; i--;)
	{
		if (weights.data[i])
		{
			vol.data[i] += values.data[i] / weights.data[i] * gamma;
		}
	}
}

void L(Volume &vol, float *L)
{
	size_t volsize = vol.width * vol.length * vol.height;
	float *Lx = (float *)malloc(sizeof(float) * volsize);
	float *Ly = (float *)malloc(sizeof(float) * volsize);
	float *Lz = (float *)malloc(sizeof(float) * volsize);

	for (int z = 0; z < vol.height; z++)
	{
		float *drez = vol.data + z * (size_t)vol.width * vol.length;
		float *drez_lx = Lx + z * (size_t)vol.width * vol.length;
		float *drez_ly = Ly + z * (size_t)vol.width * vol.length;
		float *drez_lz = Lz + z * (size_t)vol.width * vol.length;
		// float* drez_l = L + z * (size_t)vol.width * vol.length;

		for (int y = 0; y < vol.length; y++)
		{
			float *drey = drez + y * vol.width;
			float *drey_lx = drez_lx + y * vol.width;
			float *drey_ly = drez_ly + y * vol.width;
			float *drey_lz = drez_lz + y * vol.width;
			// float* drey_l = drez_l + y * vol.width;

			for (int x = 0; x < vol.width; x++)
			{
				float *drex = drey + x;
				float *L_x = drey_lx + x;
				float *L_y = drey_ly + x;
				float *L_z = drey_lz + x;
				// float* dreL = drey_l + x;

				// L

				if (x + 1 == vol.width)
				{
					*L_x = 0;
				}
				else
				{
					*L_x = *drex - *(drex + 1);
				}
				if (y + 1 == vol.length)
				{
					*L_y = 0;
				}
				else
				{
					*L_y = *drex - *(drex + vol.width);
				}
				if (z + 1 == vol.height)
				{
					*L_z = 0;
				}
				else
				{
					*L_z = *drex - *(drex + vol.width * vol.length);
				}
			}
		}
	}
	memcpy(L, Lx, sizeof(float) * volsize);
	memcpy(L + volsize, Ly, sizeof(float) * volsize);
	memcpy(L + 2 * volsize, Lz, sizeof(float) * volsize);
	free(Lx);
	free(Ly);
	free(Lz);
}

// ATb+mu*LT(u1-d1,u2-d2,u3-d3)
// 结果在htb_lt
void ATbmuLT(float *atb, float *atb_lt, float *uk, float *dk, int width, int length, int height, float mu)
{
	size_t volsize = width * length * height;
	//  std::ofstream file;
	//  file.open("atb.txt",std::ios::out);
	for (int z = 0; z < height; z++)
	{
		float *drez_atb = atb + z * (size_t)width * length; // 体的值
		float *drez_u1 = uk + z * (size_t)width * length;
		float *drez_d1 = dk + z * (size_t)width * length;
		float *drez_u2 = uk + z * (size_t)width * length + volsize;
		float *drez_d2 = dk + z * (size_t)width * length + volsize;
		float *drez_u3 = uk + z * (size_t)width * length + 2 * volsize;
		float *drez_d3 = dk + z * (size_t)width * length + 2 * volsize;
		float *drez_atblt = atb_lt + z * (size_t)width * length; // 体的值

		for (int y = 0; y < length; y++)
		{
			float *drey_atb = drez_atb + y * width;
			float *drey_u1 = drez_u1 + y * width;
			float *drey_d1 = drez_d1 + y * width;
			float *drey_u2 = drez_u2 + y * width;
			float *drey_d2 = drez_d2 + y * width;
			float *drey_u3 = drez_u3 + y * width;
			float *drey_d3 = drez_d3 + y * width;
			float *drey_atblt = drez_atblt + y * width;

			for (int x = 0; x < width; x++)
			{
				float *drex_atb = drey_atb + x;
				float *drex_u1 = drey_u1 + x;
				float *drex_d1 = drey_d1 + x;
				float *drex_u2 = drey_u2 + x;
				float *drex_d2 = drey_d2 + x;
				float *drex_u3 = drey_u3 + x;
				float *drex_d3 = drey_d3 + x;
				float *drex_atblt = drey_atblt + x; // xn

				if (x > 0)
				{																		  // LT(U)-LT(D)
					*drex_atblt += *drex_u1 - *(drex_u1 - 1) - *drex_d1 + *(drex_d1 - 1); // xn-xn+1
				}
				else
				{
					*drex_atblt += *drex_u1 - *drex_d1;
				}
				// L_y
				if (y > 0)
				{
					*drex_atblt += *drex_u2 - *(drex_u2 - width) - *drex_d2 + *(drex_d2 - width);
				}
				else
				{
					*drex_atblt += *drex_u2 - *drex_d2;
				}
				// L_z
				if (z > 0)
				{
					*drex_atblt += *drex_u3 - *drex_d3 - *(drex_u3 - width * length) + *(drex_d3 - width * length);
				}
				else
				{
					*drex_atblt += *drex_u3 - *drex_d3;
				}

				*drex_atblt = (*drex_atblt) * mu + *drex_atb;
				// test
				//  if (y == 7) {//q的第x+1列
				//  	printf("vol:%f   ", *dreL);
				//  }
			}
		}
	}
}

// ATb+mu*I^T(u1-d1,u2-d2,u3-d3)
// 结果在htb_lt
void ATbmuIT(float *atb, float *atb_lt, float *uk, float *dk, int width, int length, int height, float mu)
{
	size_t volsize = width * length * height;
	//  std::ofstream file;
	//  file.open("atb.txt",std::ios::out);
	for (int z = 0; z < height; z++)
	{
		float *drez_atb = atb + z * (size_t)width * length; // 体的值
		float *drez_u1 = uk + z * (size_t)width * length;
		float *drez_d1 = dk + z * (size_t)width * length;
		float *drez_atblt = atb_lt + z * (size_t)width * length; // 体的值

		for (int y = 0; y < length; y++)
		{
			float *drey_atb = drez_atb + y * width;
			float *drey_u1 = drez_u1 + y * width;
			float *drey_d1 = drez_d1 + y * width;
			float *drey_atblt = drez_atblt + y * width;

			for (int x = 0; x < width; x++)
			{
				float *drex_atb = drey_atb + x;
				float *drex_u1 = drey_u1 + x;
				float *drex_d1 = drey_d1 + x;
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
	// memcpy(vol.data,L,volsize);
	// free(ata);
	free(Lx);
	free(Ly);
	free(Lz);
}

// ATA(X)+mu;
void ATAmuITI(float *vol, float *ata, float *ltl, float mu,
			  int width, int length, int height)
{
	size_t volsize = width * length * height;
	size_t n;

	// ATA(vol,t,t_b,ata,num,width,length,height);

	for (int z = 0; z < height; z++)
	{
		float *drez_l = ltl + z * (size_t)width * length;
		float *drez_ATA = ata + z * (size_t)width * length;

		for (int y = 0; y < length; y++)
		{
			float *drey_l = drez_l + y * width;
			float *drey_ATA = drez_ATA + y * width;

			for (int x = 0; x < width; x++)
			{
				float *dreL = drey_l + x;
				float *dreATA = drey_ATA + x;
				*dreL = *dreATA + mu;
			}
		}
	}
}

// vol改变
/*ATAx0:ATAmuLTLlambda(x0)
  ATB:ATb+mu*LT(u1-d1,u2-d2,u3-d3)
*/
void costfun(const Point3DF &origin, Volume &vol, float *htb, Coeff coeffv[], int projnum, MrcStackM &projs)
{
	// size_t volsize = vol.length * vol.width * vol.height;
	// float *val = (float *)malloc(sizeof(float) * vol.length * vol.width);
	// float *vwt = (float *)malloc(sizeof(float) * vol.length * vol.width);
	// Slice proj(projs.X(), projs.Y());
	// Map<VectorXf> eig_val(&val[0], vol.length * vol.width);

	// // Map<VectorXf> eig_atax(&atax[0],volsize);
	// // Reproject_admm_atax(origin,vol,vol.data, coeffv,atax,proj);

	// // double cost;
	// // cost = (eig_atax-eig_htb).norm();
	// // printf("cost: %f \n",cost);
	// // free(atax);
	// Point3D coord;
	// int n;
	// double cost = 0;
	// for (int idx = 0; idx < projnum; idx++)
	// {
	// 	memset(val, 0, sizeof(float) * vol.length * vol.width);
	// 	memset(vwt, 0, sizeof(float) * vol.length * vol.width);
	// 	projs.ReadSliceZ(idx, proj.data);
	// 	Map<VectorXf> eig_b(&proj.data[0], vol.length * vol.width);
	// 	for (int z = 0; z < vol.height; z++)
	// 	{
	// 		float *vdrefz = vol.data + z * (size_t)vol.width * vol.length;
	// 		coord.z = z + vol.z;

	// 		for (int y = 0; y < vol.length; y++)
	// 		{
	// 			float *vdrefy = vdrefz + y * vol.width;
	// 			coord.y = y + vol.y;

	// 			for (int x = 0; x < vol.width; x++)
	// 			{
	// 				float *vdrefx = vdrefy + x;
	// 				coord.x = x + vol.x;
	// 				Weight wt;
	// 				ValCoef(origin, coord, coeffv[idx], &wt);

	// 				if (wt.x_min >= 0 && wt.x_min < vol.width && wt.y_min >= 0 &&
	// 					wt.y_min < vol.length)
	// 				{										 //(x_min, y_min)
	// 					n = wt.x_min + wt.y_min * vol.width; // index in reproj
	// 					val[n] += (1 - wt.x_min_del) * (1 - wt.y_min_del) * (*vdrefx);
	// 					vwt[n] += (1 - wt.x_min_del) * (1 - wt.y_min_del);
	// 				}
	// 				if ((wt.x_min + 1) >= 0 && (wt.x_min + 1) < vol.width &&
	// 					wt.y_min >= 0 && wt.y_min < vol.length)
	// 				{											 //(x_min+1, y_min)
	// 					n = wt.x_min + 1 + wt.y_min * vol.width; // index in reproj
	// 					val[n] += wt.x_min_del * (1 - wt.y_min_del) * (*vdrefx);
	// 					vwt[n] += wt.x_min_del * (1 - wt.y_min_del);
	// 				}
	// 				if (wt.x_min >= 0 && wt.x_min < vol.width && (wt.y_min + 1) >= 0 &&
	// 					(wt.y_min + 1) < vol.length)
	// 				{											   //(x_min, y_min+1)
	// 					n = wt.x_min + (wt.y_min + 1) * vol.width; // index in reproj
	// 					val[n] += (1 - wt.x_min_del) * wt.y_min_del * (*vdrefx);
	// 					vwt[n] += (1 - wt.x_min_del) * wt.y_min_del;
	// 				}
	// 				if ((wt.x_min + 1) >= 0 && (wt.x_min + 1) < vol.width &&
	// 					(wt.y_min + 1) >= 0 &&
	// 					(wt.y_min + 1) < vol.length)
	// 				{													 //(x_min+1, y_min+1)
	// 					n = (wt.x_min + 1) + (wt.y_min + 1) * vol.width; // index in reproj
	// 					val[n] += wt.x_min_del * wt.y_min_del * (*vdrefx);
	// 					vwt[n] += wt.x_min_del * wt.y_min_del;
	// 				}
	// 			}
	// 		}
	// 	}

	// 	for (int k = 0; k < vol.length * vol.width; k++)
	// 	{
	// 		val[k] = val[k] / vwt[k];
	// 		val[k] = (val[k] - proj.data[k]) * 1e-3;
	// 	}
	// 	// eig_val=eig_val-eig_b;
	// 	// eig_val=eig_val*1e-3;
	// 	// cost += eig_val.lpNorm<1>()/volsize;
	// 	cost += eig_val.norm();
	// }
	// printf("cost:%f\n", cost);
	// free(val);
	// free(vwt);
}
void applycg(Volume &vol, float *atax0, float *ATb, int numberIteration, float mu, const Point3DF &origin, Coeff coeffv[], int proj, MrcStackM &projs, float *htb)
{
	// 初始化
	size_t volsize = vol.length * vol.width * vol.height;
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
	// float *r0=(float *)malloc(sizeof(float)*volsize);//     =  ATb-vol.data   initial gradient
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
		memset(h0, 0, sizeof(float) * volsize); // h0是作为Apk 是一个中间变量
		{
			float *atax = (float *)malloc(sizeof(float) * volsize);
			memset(atax, 0, sizeof(float) * volsize);
			Reproject_admm_atax(origin, vol, p0, coeffv, atax, proj);	   // p0输入 atax输出
			ATAmuITI(p0, atax, h0, mu, vol.width, vol.length, vol.height); // ATA=ATAmuLTLlambda
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
		MPI_Allreduce(MPI_IN_PLACE, &d0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &d2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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
			vol.data[k] = vol.data[k] + alpha * p0[k];
			// if(vol.data[k]<0){
			// 	printf("k:%d\n",k);
			// }
		}
		{
			// vol.data输入 ax输出=Ax
			float *atax = (float *)malloc(sizeof(float) * volsize);
			memset(atax, 0, sizeof(float) * volsize);
			Reproject_admm_atax(origin, vol, vol.data, coeffv, atax, proj);
			ATAmuITI(vol.data, atax, ax, mu, vol.width, vol.length, vol.height);
			free(atax);
		}
		for (int k = 0; k < volsize; k++)
		{
			// r1[k]=r0[k]-al0pha * h0[k];//rk+1
			r0[k] = ax[k] - ATb[k]; // rk+1
			d1 += r0[k] * h0[k];	//
									// r02 += r0[k]*r0[k];
		}
		MPI_Allreduce(MPI_IN_PLACE, &d1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		// printf("r02:%f\n",r02);
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

void ADMM(const Point3DF &origin, MrcStackM &projs, Volume &vol, Coeff coeffv[],
		  int maxOutIter, // 最大外循环迭代次数
		  float mu)		  // 松弛变量
						  // // 可以用稀疏矩阵存放A：数据：4*体的大小， 但是切分的话其实不用太考虑大小
{
	// 	int numberIteration = 2; // cg循环次数
	// 	int pxsize = projs.X() * projs.Y();
	// 	size_t volsize = vol.length * vol.width * vol.height;

	// 	Slice proj(projs.X(), projs.Y());
	// 	Point3D coord;

	// 	float *htb = (float *)malloc(sizeof(float) * volsize); // 体
	// 	memset(htb, 0, sizeof(float) * volsize);			   // 把体置为0

	// 	for (int idx = 0; idx < projs.Z(); idx++)
	// 	{
	// 		// printf("prijection:%d\n",idx);
	// 		projs.ReadSliceZ(idx, proj.data);						 // 全部读入到proj,data中
	// 		Reproject_admm_htb(origin, vol, coeffv[idx], htb, proj); // 相当于backpro 计算结果是H^t*b
	// 	}

	// 	float *u_k = (float *)malloc(sizeof(float) * volsize * 3);
	// 	float *d_k = (float *)malloc(sizeof(float) * volsize * 3);
	// 	memset(u_k, 0, sizeof(float) * volsize * 3);
	// 	memset(d_k, 0, sizeof(float) * volsize * 3);

	// 	for (int i = 0; i < maxOutIter; i++)
	// 	{

	// 		printf("**  ADMM with Conjugate Gradient method ,iter:%d  **\n", i);
	// 		float *x0 = (float *)malloc(sizeof(float) * volsize);
	// 		float *htb_lt = (float *)malloc(sizeof(float) * volsize);

	// 		memset(htb_lt, 0, sizeof(float) * volsize);
	// 		memset(x0, 0, sizeof(float) * volsize);
	// 		// ATb+muLT(u1-d1,u2-d2,u3-d3)
	// 		ATbmuLT(htb, htb_lt, u_k, d_k, vol.width, vol.length, vol.height, mu); // 在循环内 每次更新ud
	// 		// printf("ATAmuLTL");
	// 		// vol不变 x0ATA(vol)+mu*LTL(vol,L,LTL)
	// 		{
	// 			float *atax = (float *)malloc(sizeof(float) * volsize);
	// 			memset(atax, 0, sizeof(float) * volsize);
	// 			Reproject_admm_atax(origin, vol, vol.data, coeffv, atax, projs.Z());
	// 			ATAmuLTL(vol.data, atax, x0, mu, vol.width, vol.length, vol.height); // x0里面放的是ATA(X)+LTLX 是一个初始值
	// 			free(atax);
	// 		}
	// 		// int id;
	// 		// MPI_Comm_rank(MPI_COMM_WORLD, &id);
	// 		// if(idx==0 && i==0&& id==0){
	// 		// 	std::ofstream file;
	// 		// 	file.opvoid ADMM(const Point3DF &origin, MrcStackM &projs, Volume &vol, Coeff coeffv[],
	// 		  int maxOutIter, // 最大外循环迭代次数
	// 		  float mu)		  // 松弛变量
	// // 可以用稀疏矩阵存放A：数据：4*体的大小， 但是切分的话其实不用太考虑大小
	// {
	// 	int numberIteration = 2; // cg循环次数
	// 	int pxsize = projs.X() * projs.Y();
	// 	size_t volsize = vol.length * vol.width * vol.height;

	// 	Slice proj(projs.X(), projs.Y());
	// 	Point3D coord;

	// 	float *htb = (float *)malloc(sizeof(float) * volsize); // 体
	// 	memset(htb, 0, sizeof(float) * volsize);			   // 把体置为0

	// 	for (int idx = 0; idx < projs.Z(); idx++)
	// 	{
	// 		// printf("prijection:%d\n",idx);
	// 		projs.ReadSliceZ(idx, proj.data);						 // 全部读入到proj,data中
	// 		Reproject_admm_htb(origin, vol, coeffv[idx], htb, proj); // 相当于backpro 计算结果是H^t*b
	// 	}

	// 	float *u_k = (float *)malloc(sizeof(float) * volsize * 3);
	// 	float *d_k = (float *)malloc(sizeof(float) * volsize * 3);
	// 	memset(u_k, 0, sizeof(float) * volsize * 3);
	// 	memset(d_k, 0, sizeof(float) * volsize * 3);

	// 	for (int i = 0; i < maxOutIter; i++)
	// 	{

	// 		printf("**  ADMM with Conjugate Gradient method ,iter:%d  **\n", i);
	// 		float *x0 = (float *)malloc(sizeof(float) * volsize);
	// 		float *htb_lt = (float *)malloc(sizeof(float) * volsize);

	// 		memset(htb_lt, 0, sizeof(float) * volsize);
	// 		memset(x0, 0, sizeof(float) * volsize);
	// 		// ATb+muLT(u1-d1,u2-d2,u3-d3)
	// 		ATbmuLT(htb, htb_lt, u_k, d_k, vol.width, vol.length, vol.height, mu); // 在循环内 每次更新ud
	// 		// printf("ATAmuLTL");
	// 		// vol不变 x0ATA(vol)+mu*LTL(vol,L,LTL)
	// 		{
	// 			float *atax = (float *)malloc(sizeof(float) * volsize);
	// 			memset(atax, 0, sizeof(float) * volsize);
	// 			Reproject_admm_atax(origin, vol, vol.data, coeffv, atax, projs.Z());
	// 			ATAmuLTL(vol.data, atax, x0, mu, vol.width, vol.length, vol.height); // x0里面放的是ATA(X)+LTLX 是一个初始值
	// 			free(atax);
	// 		}
	// 		// int id;
	// 		// MPI_Comm_rank(MPI_COMM_WORLD, &id);
	// 		// if(idx==0 && i==0&& id==0){
	// 		// 	std::ofstream file;
	// 		// 	file.open("vol.txt",std::ios::out);
	// 		// 	for(int k=0;k<volsize;k++){
	// 		// 		//printf("%f ",htb_lt[k]);
	// 		// 		file<<vol.data[k]<<std::endl;
	// 		// 	}
	// 		// 	file.close();
	// 		// }u_kx0, htb_lt, numberIteration, mu, origin, coeffv, projs.Z(), projs, htb); // vol已经更新->x_k 也会作为下一次cg迭代的初始值
	// 		free(x0);
	// 		free(htb_lt);
	// 		// 更新u、d
	// 		L(vol, u_k); //=L(x_k)
	// 					 //*soft: if abs(x)<a => x=0 else x=sgn(x)*(abs(x)-a)*
	// 		// u_k+1    =  softthreshold(L(x_k)+d_k,alhpha/mu)=softthreshold(u_k+d_k,10)
	// 		// d_k+1 = d_k+u_k-u_k+10
	// 		float soft = 10;
	// 		for (int k = 0; k < volsize; k++)
	// 		{
	// 			float u; // 中间变量u_k+1
	// 			if (u_k[k] + d_k[k] < -soft)
	// 			{
	// 				u = u_k[k] + d_k[k] + soft;
	// 				d_k[k] = d_k[k] + u_k[k] - u;
	// 				u_k[k] = u;
	// 			}
	// 			else if (u_k[k] + d_k[k] > soft)
	// 			{
	// 				u = u_k[k] + d_k[k] - soft;
	// 				d_k[k] = d_k[k] + u_k[k] - u;
	// 				u_k[k] = u;
	// 			}
	// 			else
	// 			{
	// 				u = 0;
	// 				d_k[k] = d_k[k] + u_k[k] - u;
	// 				u_k[k] = u;
	// 			}
	// 		}
	// 	}
	// 	free(u_k);
	// 	free(d_k);
	// 	free(htb);
	// }en("vol.txt",std::ios::out);
	// 		// 	for(int k=0;k<volsize;k++){
	// 		// 		//printf("%f ",htb_lt[k]);
	// 		// 		file<<vol.data[k]<<std::endl;
	// 		// 	}
	// 		// 	file.close();
	// 		// }u_kx0, htb_lt, numberIteration, mu, origin, coeffv, projs.Z(), projs, htb); // vol已经更新->x_k 也会作为下一次cg迭代的初始值
	// 		free(x0);
	// 		free(htb_lt);
	// 		// 更新u、d
	// 		L(vol, u_k); //=L(x_k)
	// 					 //*soft: if abs(x)<a => x=0 else x=sgn(x)*(abs(x)-a)*
	// 		// u_k+1    =  softthreshold(L(x_k)+d_k,alhpha/mu)=softthreshold(u_k+d_k,10)
	// 		// d_k+1 = d_k+u_k-u_k+10
	// 		float soft = 10;
	// 		for (int k = 0; k < volsize; k++)
	// 		{
	// 			float u; // 中间变量u_k+1
	// 			if (u_k[k] + d_k[k] < -soft)
	// 			{
	// 				u = u_k[k] + d_k[k] + soft;
	// 				d_k[k] = d_k[k] + u_k[k] - u;
	// 				u_k[k] = u;
	// 			}
	// 			else if (u_k[k] + d_k[k] > soft)
	// 			{
	// 				u = u_k[k] + d_k[k] - soft;
	// 				d_k[k] = d_k[k] + u_k[k] - u;
	// 				u_k[k] = u;
	// 			}
	// 			else
	// 			{
	// 				u = 0;
	// 				d_k[k] = d_k[k] + u_k[k] - u;
	// 				u_k[k] = u;
	// 			}
	// 		}
	// 	}
	// 	free(u_k);
	// 	free(d_k);
	// 	free(htb);
}
// void UpdateU( Volume &vol, float *u_k, const MrcStackM &projs, float eta, float lamb)
// {
// 	for (int z = 0; z < vol.height; z++)
// 	{
// 		cv::Mat slice(projs.Y(), projs.X(), CV_32F, &vol.data[z * projs.X() * projs.Y()]);
// 		slice.convertTo(slice, CV_8U, 255.0); //只接受0-255的输入

// 		cv::Mat sliceFiltered(projs.Y(), projs.X(), CV_32F, &u_k[z * projs.X() * projs.Y()]);
// 		sliceFiltered.convertTo(sliceFiltered, CV_8U, 255.0);

// 		cv::fastNlMeansDenoising(slice, sliceFiltered, eta / lamb); // eta/lamb是滤波强度
// 		//转换回去
// 		slice.convertTo(slice, CV_32F, 1.0/255.0);
// 		sliceFiltered.convertTo(sliceFiltered, CV_32F, 1.0/255.0);
// 	}
// }
// bm3d
void UpdateU(Volume &vol, float *u_k, const MrcStackM &projs, float eta, float lamb)
{
	// 计算缩放因子
	float *max_it = std::max_element(temp, temp + size);
	float *min_it = std::min_element(temp, temp + size);
	float scale = 255.0 / (*max_it - *min_it);
	float offset = -(*min_it) * scale;
	for (int z = 0; z < vol.height; z++)
	{
		cv::Mat slice(projs.Y(), projs.X(), CV_32F, &vol.data[z * projs.X() * projs.Y()]);
		slice.convertTo(slice, CV_8U, 255.0); // 只接受0-255的输入

		cv::Mat sliceFiltered(projs.Y(), projs.X(), CV_32F, &u_k[z * projs.X() * projs.Y()]);
		sliceFiltered.convertTo(sliceFiltered, CV_8U, 255.0);

		cv::fastNlMeansDenoising(slice, sliceFiltered, eta / lamb); // eta/lamb是滤波强度
		// 转换回去
		slice.convertTo(slice, CV_32F, 1.0 / 255.0);
		sliceFiltered.convertTo(sliceFiltered, CV_32F, 1.0 / 255.0);
	}
}
void UpdateD(const Volume &vol, const float *u_k, float *d_k, size_t volsize)
{
	for (int k = 0; k < volsize; k++)
	{
		d_k[k] = d_k[k] + vol.data[k] - u_k[k];
	}
}
void PnP(const Point3DF &origin, MrcStackM &projs, Volume &vol, Coeff coeffv[],
		 int maxOutIter, // 最大外循环迭代次数
		 float eta)		 // 松弛变量
// 可以用稀疏矩阵存放A：数据：4*体的大小， 但是切分的话其实不用太考虑大小
{
	float lamb = 0.05;
	int numberIteration = 2; // cg循环次数
	int pxsize = projs.X() * projs.Y();
	size_t volsize = vol.length * vol.width * vol.height;

	Slice proj(projs.X(), projs.Y());
	Point3D coord;

	float *htb = (float *)malloc(sizeof(float) * volsize); // 体
	memset(htb, 0, sizeof(float) * volsize);			   // 把体置为0

	for (int idx = 0; idx < projs.Z(); idx++)
	{
		// printf("prijection:%d\n",idx);
		projs.ReadSliceZ(idx, proj.data);						 // 全部读入到proj,data中
		Reproject_admm_htb(origin, vol, coeffv[idx], htb, proj); // 相当于backpro 计算结果是H^t*b
	}

	float *u_k = (float *)malloc(sizeof(float) * volsize);
	float *d_k = (float *)malloc(sizeof(float) * volsize);
	memset(u_k, 0, sizeof(float) * volsize);
	memset(d_k, 0, sizeof(float) * volsize);

	for (int i = 0; i < maxOutIter; i++)
	{

		printf("**  PnP with Conjugate Gradient method ,iter:%d  **\n", i);
		float *x0 = (float *)malloc(sizeof(float) * volsize);
		float *htb_lt = (float *)malloc(sizeof(float) * volsize);

		memset(htb_lt, 0, sizeof(float) * volsize);
		memset(x0, 0, sizeof(float) * volsize);
		// ATb+muIT(u-d)
		ATbmuIT(htb, htb_lt, u_k, d_k, vol.width, vol.length, vol.height, lamb); // 在循环内 每次更新ud
		// printf("ATAmuLTL");
		// vol不变 x0ATA(vol)+mu*LTL(vol,L,LTL)
		{
			float *atax = (float *)malloc(sizeof(float) * volsize);
			memset(atax, 0, sizeof(float) * volsize);
			Reproject_admm_atax(origin, vol, vol.data, coeffv, atax, projs.Z());
			ATAmuITI(vol.data, atax, x0, lamb, vol.width, vol.length, vol.height); // x0里面放的是ATA(X)+LTLX 是一个初始值
			free(atax);
		}
		// int id;
		// MPI_Comm_rank(MPI_COMM_WORLD, &id);
		// if(idx==0 && i==0&& id==0){
		// 	std::ofstream file;
		// 	file.open("vol.txt",std::ios::out);
		// 	for(int k=0;k<volsize;k++){
		// 		//printf("%f ",htb_lt[k]);
		// 		file<<vol.data[k]<<std::endl;
		// 	}
		// 	file.close();
		// }
		// vol就是ATAx0
		printf("**  Conjugate Gradient Method  **\n");
		applycg(vol, x0, htb_lt, numberIteration, lamb, origin, coeffv, projs.Z(), projs, htb); // vol已经更新->x_k 也会作为下一次cg迭代的初始值
		free(x0);
		free(htb_lt);
		// 更新u、d
		// 更新u 使用去噪算子NLM
		UpdateU(vol, u_k, projs, eta, lamb);
		// 更新d
		UpdateD(vol, u_k, d_k, volsize);
	}
	free(u_k);
	free(d_k);
	free(htb);
}

void SART(const Point3DF &origin, MrcStackM &projs, Volume &vol, Coeff coeffv[],
		  int iteration // 迭代次数
		  ,
		  float gamma)
{ // 松弛变量ADMMs.Y());

	// for (int i = 0; i < iteration; i++)
	// {
	// 	for (int idx = 0; idx < projs.Z(); idx++)
	// 	{
	// 		memset(reproj_val.data, 0, sizeof(float) * pxsize);
	// 		memset(reproj_wt.data, 0, sizeof(float) * pxsize);

	// 		Reproject(origin, vol, coeffv[idx], reproj_val, reproj_wt); // 一个reproj_wt

	// 		MPI_Allreduce(MPI_IN_PLACE, reproj_val.data, pxsize, MPI_FLOAT, MPI_SUM,
	// 					  MPI_COMM_WORLD);
	// 		// MPI_IN_PLACE:说明当前进程既发送又接受数据，而且要发送的数据和在要接收的数据的保存在同一内存。
	// 		MPI_Allreduce(MPI_IN_PLACE, reproj_wt.data, pxsize, MPI_FLOAT, MPI_SUM,
	// 					  MPI_COMM_WORLD);

	// 		printf("SART begin to read %d projection for %d z-coordinate (iteraion "
	// 			   "%d)\n",
	// 			   idx, vol.z, i);
	// 		projs.ReadSliceZ(idx, projection.data);

	// 		for (int n = 0; n < pxsize; n++)
	// 		{
	// 			if (reproj_wt.data[n] != 0)
	// 			{
	// 				reproj_val.data[n] /= reproj_wt.data[n];
	// 			}
	// 			reproj_val.data[n] = projection.data[n] - reproj_val.data[n];
	// 		}

	// 		UpdateVolumeByProjDiff(origin, reproj_val, vol, gamma, coeffv[idx]);
	// 	}
	// }
}

void SIRT(const Point3DF &origin, MrcStackM &projs, Volume &vol, Coeff coeffv[],
		  int iteration, float gamma)
{

	int pxsize = projs.X() * projs.Y();
	Slice reproj_val(projs.X(), projs.Y()); // reprojection value
	Slice reproj_wt(projs.X(), projs.Y());	// reprojection weight
	Slice projection(projs.X(), projs.Y());

	Volume valvol(vol.x, vol.y, vol.z, vol.length, vol.width, vol.height);
	Volume wtvol(vol.x, vol.y, vol.z, vol.length, vol.width, vol.height);

	for (int i = 0; i < iteration; i++)
	{ // 输入的迭代次数

		memset(valvol.data, 0,
			   sizeof(float) * valvol.length * valvol.width * valvol.height);
		memset(wtvol.data, 0,
			   sizeof(float) * wtvol.length * wtvol.width * wtvol.height);
		// printf("3");//这里不行

		for (int idx = 0; idx < projs.Z(); idx++) // projs.Z 是角度的数量
		{
			memset(reproj_val.data, 0, sizeof(float) * pxsize);
			memset(reproj_wt.data, 0, sizeof(float) * pxsize);

			Reproject(origin, vol, coeffv[idx], reproj_val,
					  reproj_wt); // vol is not changed during iteration

			MPI_Allreduce(MPI_IN_PLACE, reproj_val.data, pxsize, MPI_FLOAT, MPI_SUM,
						  MPI_COMM_WORLD);
			MPI_Allreduce(MPI_IN_PLACE, reproj_wt.data, pxsize, MPI_FLOAT, MPI_SUM,
						  MPI_COMM_WORLD);

			printf("SIRT begin to read %d projection for %d z-coordinate (iteraion "
				   "%d)\n",
				   idx, vol.z, i);

			projs.ReadSliceZ(idx, projection.data);

			for (int n = 0; n < pxsize; n++)
			{
				if (reproj_wt.data[n])
				{
					reproj_val.data[n] /= reproj_wt.data[n];
				}
				reproj_val.data[n] = projection.data[n] - reproj_val.data[n];
			}

			UpdateWeightsByProjDiff(origin, reproj_val, valvol, wtvol, coeffv[idx]);
		}

		UpdateVolumeByWeights(vol, valvol, wtvol, gamma);
	}
}

struct SysInfo
{
	int id;
	int procs;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int namelen;
};

int ATOM(options &opt, int myid, int procs)
{
	MrcStackM projs, mrcvol;
	if (!projs.ReadFile(opt.input))
	{
		printf("File %s cannot access.\n", opt.input);

		return -1;
	}

	if (myid == 0)
	{
		projs.ReadHeader();
	}
	MPI_Bcast(&(projs.header), sizeof(MRCheader), MPI_CHAR, 0, MPI_COMM_WORLD);

	mrcvol.InitializeHeader();
	mrcvol.SetSize(projs.X(), projs.Y(), opt.thickness);

	std::vector<float> angles;
	ReadAngles(angles, opt.angle);

	std::vector<float> xangles;

	if (opt.xangle[0] != '\0')
	{
		ReadAngles(xangles, opt.xangle);
	}
	else
	{
		xangles.resize(angles.size(), 0.0);
	}

	std::vector<Coeff> params;
	TranslateAngleToCoefficients(angles, xangles, params);

	Geometry geo;
	geo.offset = opt.offset;
	geo.pitch_angle = opt.pitch_angle;
	geo.zshift = opt.zshift;

	DecorateCoefficients(params, geo);

	int height;
	int zrem = mrcvol.Z() % procs;
	int volz; // the start slice of reproject per process

	if (myid < zrem)
	{
		height = mrcvol.Z() / procs + 1;
		volz = height * myid;
	}
	else
	{
		height = mrcvol.Z() / procs;
		volz = height * myid + zrem;
	}

	// if (myid == 0)
	// { // 第一个切片
	// 	// 起点不变volz=volz;
	// 	height += 1; // 后面多切一片
	// }
	// else if (myid == procs - 1) // 最后一个切片
	// {
	// 	volz = volz - 1; // 前面多切一片
	// 	height += 1;
	// }
	// else
	// {
	// 	// volz=volz-1;
	// 	height += 1;, 0, vol.width * vol.l
	// }

	Volume vol(0, 0, volz, mrcvol.Y(), mrcvol.X(), height);
	std::cout << myid << ": (" << vol.x << "," << vol.y << "," << vol.z << ")"
			  << "&(" << vol.width << "," << vol.length << "," << vol.height
			  << ")" << std::endl;

	Point3DF origin;

	// AssignValue 原来的会丢失小数
	origin.x = mrcvol.X() * .5;
	origin.y = mrcvol.Y() * .5;
	origin.z = mrcvol.Z() * .5;

	if (myid == 0)
	{
		printf("origin.x is %f, origin.y is %f, origin.z is %f\n", origin.x,
			   origin.y, origin.z);
	}

	/***************************reprojection************************/
	if (opt.method == "RP")
	{

		{ // 分片读入数据
			int j = 0;
			for (j = 0; vol.z + j + 5 < vol.z + vol.height; j += 5)
				projs.ReadBlock(vol.z + j, vol.z + j + 5, 'z', (vol.data + (size_t)vol.width * vol.length * j));
			projs.ReadBlock(vol.z + j, vol.z + vol.height, 'z', (vol.data + (size_t)vol.width * vol.length * j));

			// projs.ReadBlock(vol.z, vol.z + vol.height, 'z', (vol.data));
		}

		mrcvol.SetSize(projs.X(), projs.Y(), params.size());
		if (opt.output[0] == '\0')
		{
			const char *name = "reproj.mrc";
			memcpy(opt.output, name, sizeof(char) * 11);
		}

		mrcvol.WriteToFile(opt.output);

		if (myid == 0)
		{
			mrcvol.WriteHeader();
		}

		int pxsize = projs.X() * projs.Y();
		Slice reproj_val(projs.X(), projs.Y()); // reprojection value
		Slice reproj_wt(projs.X(), projs.Y());	// reprojection weight

		for (int idx = 0; idx < params.size(); idx++)
		{
			memset(reproj_val.data, 0, sizeof(float) * pxsize);
			memset(reproj_wt.data, 0, sizeof(float) * pxsize);

			Reproject(origin, vol, params[idx], reproj_val, reproj_wt); // reproject是输入vol输出投影的
			MPI_Allreduce(MPI_IN_PLACE, reproj_val.data, pxsize, MPI_FLOAT, MPI_SUM,
						  MPI_COMM_WORLD);
			MPI_Allreduce(MPI_IN_PLACE, reproj_wt.data, pxsize, MPI_FLOAT, MPI_SUM,
						  MPI_COMM_WORLD);

			for (int n = 0; n < pxsize; n++)
			{
				if (reproj_wt.data[n] != 0)
				{
					reproj_val.data[n] /= reproj_wt.data[n];
				}
			}

			if (myid == 0)
			{
				printf("RP begin to write %d projection\n", idx);
				mrcvol.WriteSlice(idx, 'z', reproj_val.data);
			}
		}

		if (myid == 0)
		{
			mrcvol.UpdateHeader();
		}

		projs.Close();
		mrcvol.Close();

		return 0;
	}

	/**********************reconstruction along Z-axis*******************/
	if (opt.method == "BPT")
	{
		BackProject(origin, projs, vol, &params[0]);
	}
	else if (opt.method == "SART")
	{
		const char *reference = opt.initial;
		if (reference[0] == '\0')
		{
			memset(vol.data, 0, vol.width * vol.length * vol.height * sizeof(float));
		}
		else
		{
			MrcStackM init;
			init.ReadFile(reference);
			init.ReadHeader();
			{
				int j = 0;
				for (j = 0; vol.z + j + 5 < vol.z + vol.height; j += 5)
					projs.ReadBlock(vol.z + j, vol.z + j + 5, 'z', (vol.data + (size_t)vol.width * vol.length * j));
				projs.ReadBlock(vol.z + j, vol.z + vol.height, 'z', (vol.data + (size_t)vol.width * vol.length * j));
			}
			init.Close();
		}

		SART(origin, projs, vol, &params[0], opt.iteration, opt.gamma);
	}
	else if (opt.method == "SIRT")
	{
		const char *reference = opt.initial;
		if (reference[0] == '\0')
		{
			memset(vol.data, 0, vol.width * vol.length * vol.height * sizeof(float));
			// printf("1");
		}
		else
		{
			MrcStackM init;
			init.ReadFile(reference);
			init.ReadHeader();
			{
				int j = 0;
				for (j = 0; vol.z + j + 5 < vol.z + vol.height; j += 5)
					projs.ReadBlock(vol.z + j, vol.z + j + 5, 'z', (vol.data + (size_t)vol.width * vol.length * j));
				projs.ReadBlock(vol.z + j, vol.z + vol.height, 'z', (vol.data + (size_t)vol.width * vol.length * j));
				// printf("2");
				//  projs.ReadBlock(vol.z, vol.z + vol.height, 'z', (vol.data));
			}
			init.Close();
		}
		SIRT(origin, projs, vol, &params[0], opt.iteration, opt.gamma);
	}
	else if (opt.method == "ADMM")
	{
		// printf("opt.gamma:%f\n",opt.gamma);
		const char *reference = opt.initial;
		// vol.height=vol.height+1;
		if (reference[0] == '\0')
		{
			memset(vol.data, 0, vol.width * vol.length * vol.height * sizeof(float));
			// printf("1");
		}
		else
		{
			MrcStackM init;
			init.ReadFile(reference);
			init.ReadHeader();
			{
				int j = 0;
				for (j = 0; vol.z + j + 5 < vol.z + vol.height; j += 5)
					projs.ReadBlock(vol.z + j, vol.z + j + 5, 'z', (vol.data + (size_t)vol.width * vol.length * j));
				projs.ReadBlock(vol.z + j, vol.z + vol.height, 'z', (vol.data + (size_t)vol.width * vol.length * j));
				// printf("2");
				//  projs.ReadBlock(vol.z, vol.z + vol.height, 'z', (vol.data));
			}
			init.Close();
		}
		// BackProject(origin, projs, vol, &params[0]);
		ADMM(origin, projs, vol, &params[0], opt.iteration, opt.gamma);
		// vol.height=vol.height-1;
	}
	else if (opt.method == "PnP")
	{
		// printf("opt.gamma:%f\n",opt.gamma);
		const char *reference = opt.initial;
		// vol.height=vol.height+1;
		if (reference[0] == '\0')
		{
			memset(vol.data, 0, vol.width * vol.length * vol.height * sizeof(float));
			// printf("1");
		}
		else
		{
			MrcStackM init;
			init.ReadFile(reference);
			init.ReadHeader();
			{
				int j = 0;
				for (j = 0; vol.z + j + 5 < vol.z + vol.height; j += 5)
					projs.ReadBlock(vol.z + j, vol.z + j + 5, 'z', (vol.data + (size_t)vol.width * vol.length * j));
				projs.ReadBlock(vol.z + j, vol.z + vol.height, 'z', (vol.data + (size_t)vol.width * vol.length * j));
				// printf("2");
				//  projs.ReadBlock(vol.z, vol.z + vol.height, 'z', (vol.data));
			}
			init.Close();
		}
		// BackProject(origin, projs, vol, &params[0]);
		PnP(origin, projs, vol, &params[0], opt.iteration, opt.gamma);
		// vol.height=vol.height-1;
	}

	mrcvol.WriteToFile(opt.output);

	if (myid == 0)
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
		std::cout << myid << ": (" << vol.x << "," << vol.y << "," << vol.z << ")"
				  << "&(" << vol.width << "," << vol.length << "," << vol.height
				  << ")" << std::endl;

		// 分片写入
		int j = 0;
		for (j = 0; vol.z + j + 5 < vol.z + vol.height; j += 5)
		{
			mrcvol.WriteBlock(vol.z + j, vol.z + j + 5, 'z', (vol.data + (size_t)vol.width * vol.length * j));
		}
		mrcvol.WriteBlock(vol.z + j, vol.z + vol.height, 'z', (vol.data + (size_t)vol.width * vol.length * j));
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (myid == 0)
	{
		mrcvol.UpdateHeader();
	}

	projs.Close();
	mrcvol.Close();

	return 0;
}

int main(int argc, char *argv[])
{
	SysInfo info;

	int err;
	err = MPI_Init(&argc, &argv);
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
	MPI_Comm_size(MPI_COMM_WORLD, &(info.procs));
	MPI_Get_processor_name(info.processor_name, &(info.namelen));

	options opts;
	InitOpts(&opts);

	if (GetOpts(argc, argv, &opts) <= 0)
	{
		EX_TRACE("***WRONG INPUT.\n");
		return -1;
	}

	if (info.id == 0)
	{
		PrintOpts(opts);
	}

	ATOM(opts, info.id, info.procs);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize(); // parallel finish

	return 0;
}
