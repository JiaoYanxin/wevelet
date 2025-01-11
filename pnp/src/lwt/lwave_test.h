//============================================================================
// Name : 1D/2D Wavelet Transform
// Author : Rafat Hussain
// Version :
// Copyright : GNU GPL License
// Description : LiftWave Wavelet Library Component
//============================================================================
/*
 * Copyright (c) 2012 Rafat Hussain
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#ifndef LWAVE_H
#define LWAVE_H
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

#include "lift.h"
#include "alg.h"

using namespace std;

template <class T>
struct IsInt
{
	static const bool value = false;
};

template <>
struct IsInt<int>
{
	static const bool value = true;
};

template <>
struct IsInt<short>
{
	static const bool value = true;
};

template <>
struct IsInt<long>
{
	static const bool value = true;
};

// 是固定横着算的
template <typename T>
class lwt
{
	vector<T> cA, cD;	   // ca 低频 cd 高频
	vector<int> cD_length; // 存储每个高频成分的长度
	int level;			   // 表示小波变换的层数或尺度

public:
	lwt(vector<T> &signal, liftscheme &lft) // 构造函数
	{
		level = 1;			  // 初始化
		vector<double> coeff; // 1 -0.5
		vector<int> lenv;	  //{1, 0, 1, 0};
		string lat;			  // dp
		double K;
		lft.getScheme(coeff, lenv, lat, K);

		// Number Of Liftin Stages N
		int N = lat.size(); // 2
		vector<T> sl, dl;	// sl奇数 dl偶数
		split(signal, sl, dl);
		int cume_coeff = 0;

		for (int i = 0; i < N; i++)
		{
			char lft_type = lat.at(i);	   // d/p
			vector<double> filt;		   // 用于存储当前阶段的滤波器系数。
			int len_filt = lenv[2 * i];	   // 当前阶段滤波器系数的长度。i=0 1
			int max_pow = lenv[2 * i + 1]; // 获取当前阶段的最大指数 i=0 0

			for (int j = 0; j < len_filt; j++) // 1
			{
				filt.push_back(coeff[cume_coeff + j]); // filt=1 -0.5
			}
			cume_coeff = cume_coeff + len_filt; // 更新累积系数索引，为下一个阶段的滤波器系数准备。 1

			if (lft_type == 'd') // 如果当前阶段的类型为d 对奇数进行处理 放入偶数 是update
			{

				for (int len_dl = 0; len_dl < (int)dl.size(); len_dl++) // 对于奇数序列 遍历里面每一个值
				{
					double temp = 0.0;
					for (int lf = 0; lf < len_filt; lf++) // 1
					{									  // 检查是否在有效范围内
						if ((len_dl + max_pow - lf) >= 0 && (len_dl + max_pow - lf) < (int)sl.size())
						{														// 和滤波器的对应位置相乘 累加
							temp = temp + filt[lf] * sl[len_dl + max_pow - lf]; // 1*sl[1]-0.5*sl[0]
						}
					}
					dl[len_dl] = dl[len_dl] - (T)temp;
					// dl[0]=dl[0]-sl[0]
					// dl[1]=dl[1]-(1*sl[1]-0.5*sl[0])
					//
				}
			}
			else if (lft_type == 'p') // 如果当前阶段的类型为p 是预测 对偶数部分进行处理 放入奇数里
			{

				for (int len_sl = 0; len_sl < (int)dl.size(); len_sl++)
				{
					double temp = 0.0;
					for (int lf = 0; lf < len_filt; lf++)
					{
						if ((len_sl + max_pow - lf) >= 0 && (len_sl + max_pow - lf) < (int)dl.size())
						{
							temp = temp + filt[lf] * dl[len_sl + max_pow - lf];
						}
					}
					sl[len_sl] = sl[len_sl] - (T)temp;
				}
			}
		}
		double K1 = 1.0 / K;
		if (!IsInt<T>::value)
		{
			vecmult(sl, K1);
			vecmult(dl, K);
		}

		cA = sl;
		cD = dl;
		cD_length.clear();
		cD_length.push_back((int)cD.size());
	}

	lwt(vector<T> &signal, string &name)
	{
		level = 1;
		liftscheme lft(name);		  // 决定阶段dp
		lwt<T> wavelift(signal, lft); // 调用第一个构造函数
		vector<T> sx, dx;
		wavelift.getCoeff(sx, dx); // wavelift计算的CA CD拷贝到 sx dx
		cA = sx;				   // 拷贝到这个对象的ca cd
		cD = dx;
		cD_length.clear();
		cD_length.push_back((int)cD.size());
	}

	lwt(vector<T> &signal, liftscheme &lft, int &J)
	{
		/*	int Max_Iter;
			Max_Iter = (int) ceil(log( double(signal.size()))/log (2.0)) - 1;

			if ( Max_Iter < J) {
				J = Max_Iter;

			}*/

		vector<T> temp = signal;
		vector<T> det, temp2;
		vector<int> len_det;

		for (int iter = 0; iter < J; iter++)
		{
			lwt jlevel(temp, lft);
			jlevel.getCoeff(temp, temp2);
			int len_d = temp2.size();
			det.insert(det.begin(), temp2.begin(), temp2.end());
			len_det.insert(len_det.begin(), len_d);
		}
		cA = temp;
		cD = det;
		cD_length = len_det;
		level = J;
	}

	lwt(vector<T> &signal, string &name, int &J)
	{
		liftscheme lft(name);
		lwt<T> wavelift(signal, lft, J);
		vector<T> sx, dx;
		wavelift.getCoeff(sx, dx);
		cA = sx;
		cD = dx;
		vector<int> cdlen;
		wavelift.getDetailVec(cdlen);
		cD_length = cdlen;
		level = J;
	}

	void getCoeff(vector<T> &appx, vector<T> &det)
	{
		appx = cA;
		det = cD;
	}

	void getDetailVec(vector<int> &detvec)
	{
		detvec = cD_length;
	}

	int getLevels()
	{
		return level;
	}

	virtual ~lwt()
	{
	}
};

template <typename T>
class ilwt
{
	vector<T> signal;

public:
	ilwt(vector<T> &sl, vector<T> &dl, liftscheme &lft)
	{
		vector<double> coeff;
		vector<int> lenv;
		string lat;
		double K;
		lft.getScheme(coeff, lenv, lat, K);
		// vector<T> sl,dl;
		// sl=cA;
		// dl=cD;
		double K1 = 1.0 / K;
		if (!IsInt<T>::value)
		{
			vecmult(sl, K);
			vecmult(dl, K1);
		}

		int N = lat.size();

		int cume_coeff = coeff.size();

		for (int i = N - 1; i >= 0; i--)
		{
			char lft_type = lat.at(i);
			vector<double> filt;
			int len_filt = lenv[2 * i];
			int max_pow = lenv[2 * i + 1];

			cume_coeff = cume_coeff - len_filt;

			for (int j = 0; j < len_filt; j++)
			{
				filt.push_back(coeff[cume_coeff + j]);
			}

			if (lft_type == 'd')
			{

				for (int len_dl = 0; len_dl < (int)dl.size(); len_dl++)
				{
					double temp = 0.0;
					for (int lf = 0; lf < len_filt; lf++)
					{
						if ((len_dl + max_pow - lf) >= 0 && (len_dl + max_pow - lf) < (int)sl.size())
						{
							temp = temp + filt[lf] * sl[len_dl + max_pow - lf];
						}
					}
					dl[len_dl] = dl[len_dl] + (T)temp;
				}
			}
			else if (lft_type == 'p')
			{

				for (int len_sl = 0; len_sl < (int)dl.size(); len_sl++)
				{
					double temp = 0.0;
					for (int lf = 0; lf < len_filt; lf++)
					{
						if ((len_sl + max_pow - lf) >= 0 && (len_sl + max_pow - lf) < (int)dl.size())
						{
							temp = temp + filt[lf] * dl[len_sl + max_pow - lf];
						}
					}
					sl[len_sl] = sl[len_sl] + (T)temp;
				}
			}
		}
		vector<T> idwt_oup;
		merge(idwt_oup, sl, dl);

		signal = idwt_oup;
	}

	ilwt(vector<T> &sl, vector<T> &dl, string &name)
	{
		liftscheme lft(name);
		ilwt<T> wavelift(sl, dl, lft);
		vector<T> sigx;
		wavelift.getSignal(sigx);
		signal = sigx;
	}

	ilwt(lwt<T> &wt, liftscheme &lft)
	{
		int J = wt.getLevels();
		vector<T> sl, dl;
		wt.getCoeff(sl, dl);
		vector<int> detv;

		wt.getDetailVec(detv);
		int total = 0;

		for (int i = 0; i < J; i++)
		{
			vector<T> temp, temp2;
			for (int j = 0; j < (int)detv[i]; j++)
			{
				temp.push_back(dl[total + j]);
			}
			total = total + (int)detv[i];
			ilwt<T> iwt(sl, temp, lft);
			iwt.getSignal(temp2);
			sl = temp2;
		}
		signal = sl;
	}

	ilwt(lwt<T> &wt, string &name)
	{
		liftscheme lft(name);
		ilwt<T> wavelift(wt, lft);
		vector<T> sigx;
		wavelift.getSignal(sigx);
		signal = sigx;
	}

	void getSignal(vector<T> &sig)
	{
		sig = signal;
	}

	virtual ~ilwt()
	{
	}
};

template <typename T>
class lwt2
{
	vector<T> cLL, cLH, cHL, cHH;
	int rowLL, colLL;
	int rowLH, colLH;
	int rowHL, colHL;
	int rowHH, colHH;
	int level;
	vector<int> coef_lengths;

public:
	lwt2(vector<T> &signal, int rows, int cols, liftscheme &lft)
	{
		vector<T> L, H;
		int rows_L, cols_L, rows_H, cols_H;
		// 行不变？ 列变 因为对一行进行处理
		rows_L = rows;
		rows_H = rows;
		for (int i = 0; i < rows; i++)
		{
			vector<T> temp;
			// 给temp赋值：signal 会替换掉条原有位置的值
			temp.assign(signal.begin() + i * cols, signal.begin() + (i + 1) * cols);
			lwt<T> lwt1(temp, lft); // 调用的一维的lwt的构造函数
			vector<T> a, d;
			lwt1.getCoeff(a, d); // 获得系数 拷贝到a d
			// 将向量 a 中的所有元素插入到向量 L 的末尾
			L.insert(L.end(), a.begin(), a.end());
			H.insert(H.end(), d.begin(), d.end());
			if (i == 0)
			{
				cols_L = a.size();
				cols_H = d.size();
			}
		}

		// 进行转置操作 得到LT
		vector<T> LT, HT;
		transpose(L, rows_L, cols_L, LT);
		transpose(H, rows_H, cols_H, HT);
		int rows_ll, cols_ll, rows_lh, cols_lh;

		// Remove cout
		vector<T> LL, LH;

		// Low Pass Stage
		cols_ll = cols_L;
		cols_lh = cols_L;

		for (int i = 0; i < cols_L; i++)
		{
			vector<T> temp;

			temp.assign(LT.begin() + i * rows_L, LT.begin() + (i + 1) * rows_L);
			// 对刚刚得到的数据再进行变换
			// 一个经过转置的二维数据集的每一列执行一维小波变换
			// 二维：先对行执行一维变换，再对结果的列执行一维变换。
			lwt<T> lwt1(temp, lft); // 调用的一维的lwt的构造函数
			vector<T> a, d;
			lwt1.getCoeff(a, d);
			LL.insert(LL.end(), a.begin(), a.end());
			LH.insert(LH.end(), d.begin(), d.end());
			if (i == 0)
			{
				rows_ll = a.size();
				rows_lh = d.size();
			}
		}

		int rows_hl, cols_hl, rows_hh, cols_hh;

		vector<T> HL, HH;

		// High Pass Stage
		cols_hl = cols_H;
		cols_hh = cols_H;
		for (int i = 0; i < cols_H; i++)
		{
			vector<T> temp;
			temp.assign(HT.begin() + i * rows_H, HT.begin() + (i + 1) * rows_H);
			lwt<T> lwt1(temp, lft);
			vector<T> a, d;
			lwt1.getCoeff(a, d);
			HL.insert(HL.end(), a.begin(), a.end());
			HH.insert(HH.end(), d.begin(), d.end());
			if (i == 0)
			{
				rows_hl = a.size();
				rows_hh = d.size();
			}
		}

		// cLL=LL;cLH=LH;cHL=HL;cHH=HH;
		transpose(LL, cols_ll, rows_ll, cLL);
		transpose(LH, cols_lh, rows_lh, cLH);
		transpose(HL, cols_hl, rows_hl, cHL);
		transpose(HH, cols_hh, rows_hh, cHH);
		rowLL = rows_ll;
		rowLH = rows_lh;
		rowHL = rows_hl;
		rowHH = rows_hh;
		colLL = cols_ll;
		colLH = cols_lh;
		colHL = cols_hl;
		colHH = cols_hh;
		level = 1;
		coef_lengths.push_back(rowLL);
		coef_lengths.push_back(colLL);
		coef_lengths.push_back(rowLH);
		coef_lengths.push_back(colLH);
		coef_lengths.push_back(rowHL);
		coef_lengths.push_back(colHL);
		coef_lengths.push_back(rowHH);
		coef_lengths.push_back(colHH);
	}

	lwt2(const vector<T> &signal, int rows, int cols, liftscheme &lft, int J)
	{
		vector<T> A1, B1, C1, D1;
		vector<int> siglen;
		vector<T> temp_signal = signal;
		for (int i = 0; i < J; i++)
		{
			lwt2<T> wt2(temp_signal, rows, cols, lft);
			vector<T> tempA, tempB, tempC, tempD;
			vector<int> temp_siglen;
			wt2.getCoef(tempA, tempB, tempC, tempD);
			temp_signal = tempA;
			A1 = tempA;
			B1.insert(B1.begin(), tempB.begin(), tempB.end());
			C1.insert(C1.begin(), tempC.begin(), tempC.end());
			D1.insert(D1.begin(), tempD.begin(), tempD.end());

			wt2.getDim(temp_siglen);
			rows = temp_siglen[0];
			cols = temp_siglen[1];
			if (i == J - 1)
			{
				siglen.insert(siglen.begin(), temp_siglen.begin(), temp_siglen.end());
			}
			else
			{
				siglen.insert(siglen.begin(), temp_siglen.begin() + 2, temp_siglen.end());
			}
		}
		cLL = A1;
		cLH = B1;
		cHL = C1;
		cHH = D1;
		level = J;
		coef_lengths = siglen;
	}

	void getCoef(vector<T> &aLL, vector<T> &aLH, vector<T> &aHL, vector<T> &aHH)
	{
		aLL = cLL;
		aLH = cLH;
		aHL = cHL;
		aHH = cHH;
	}

	void getDim(vector<int> &dimvec)
	{

		dimvec = coef_lengths;
	}

	int getLevels()
	{
		return level;
	}

	void getDetails(string align, int slevel, vector<T> &det_vec, vector<int> &det_len)
	{
		int J = level;
		int lev;
		if (slevel > J)
		{
			cout << " Decomposition has only " << J << " levels" << endl;
			exit(1);
		}
		else
		{
			lev = J - slevel;
		}

		vector<int> sig_vec = coef_lengths;
		vector<T> A1, B1, C1, D1;
		A1 = cLL;
		B1 = cLH;
		C1 = cHL;
		D1 = cHH;
		int total = 0;

		if (align == "LH" || align == "lh")
		{

			for (int i = 0; i < lev; i++)
			{
				total = total + sig_vec[2 + i * 6] * sig_vec[3 + i * 6];
			}
			det_vec.assign(B1.begin() + total, B1.begin() + total + sig_vec[2 + lev * 6] * sig_vec[3 + lev * 6]);
			det_len.push_back(sig_vec[2 + lev * 6]);
			det_len.push_back(sig_vec[3 + lev * 6]);
		}
		else if (align == "HL" || align == "hl")
		{

			for (int i = 0; i < lev; i++)
			{
				total = total + sig_vec[4 + i * 6] * sig_vec[5 + i * 6];
			}
			det_vec.assign(C1.begin() + total, C1.begin() + total + sig_vec[4 + lev * 6] * sig_vec[5 + lev * 6]);
			det_len.push_back(sig_vec[4 + lev * 6]);
			det_len.push_back(sig_vec[5 + lev * 6]);
		}
		else if (align == "HH" || align == "hh")
		{

			for (int i = 0; i < lev; i++)
			{
				total = total + sig_vec[6 + i * 6] * sig_vec[7 + i * 6];
			}
			det_vec.assign(D1.begin() + total, D1.begin() + total + sig_vec[6 + lev * 6] * sig_vec[7 + lev * 6]);
			det_len.push_back(sig_vec[6 + lev * 6]);
			det_len.push_back(sig_vec[7 + lev * 6]);
		}
		else
		{
			cout << "Accepted filter stages are LH or lh, HL or hl and HH or hh" << endl;
			exit(1);
		}
	}

	virtual ~lwt2()
	{
	}
};

template <typename T>
class ilwt2
{
	vector<T> signal;
	int oup_row, oup_col;

public:
	ilwt2(vector<T> &A, vector<T> &H, vector<T> &V, vector<T> &D, vector<int> &length, liftscheme &lft)
	{
		int rows_LL = length[0];
		int cols_L = length[1];
		int rows_LH = length[2];

		// 先对y
		vector<T> AT, HT, VT, DT;

		transpose(A, length[0], length[1], AT);
		transpose(H, length[2], length[3], HT);
		transpose(V, length[4], length[5], VT);
		transpose(D, length[6], length[7], DT);

		// 计算LL LH->L
		vector<T> L;
		int rows_L;

		for (int i = 0; i < cols_L; i++)
		{
			vector<T> temp1, temp2;
			temp1.assign(AT.begin() + i * rows_LL, AT.begin() + (i + 1) * rows_LL);
			temp2.assign(HT.begin() + i * rows_LH, HT.begin() + (i + 1) * rows_LH);
			ilwt<T> iwt(temp1, temp2, lft);
			vector<T> sig;
			iwt.getSignal(sig);
			L.insert(L.end(), sig.begin(), sig.end());
			if (i == 0)
			{
				rows_L = (int)sig.size();
			}
		}

		// 计算HL HH->H
		int cols_H = length[5];

		int rows_HL = length[4];
		int rows_HH = length[6];
		vector<T> H1;
		int rows_H;

		for (int i = 0; i < cols_H; i++)
		{
			vector<T> temp1, temp2;
			temp1.assign(VT.begin() + i * rows_HL, VT.begin() + (i + 1) * rows_HL);
			temp2.assign(DT.begin() + i * rows_HH, DT.begin() + (i + 1) * rows_HH);
			ilwt<T> iwt(temp1, temp2, lft);
			vector<T> sig;
			iwt.getSignal(sig);
			H1.insert(H1.end(), sig.begin(), sig.end());
			if (i == 0)
			{
				rows_H = (int)sig.size();
			}
		}

		// 转置回去
		vector<T> L2;
		transpose(L, cols_L, rows_L, L2);
		vector<T> H2;
		transpose(H1, cols_H, rows_H, H2);

		vector<T> oup;
		int cx;

		for (int i = 0; i < rows_L; i++)
		{
			vector<T> temp1, temp2;
			temp1.assign(L2.begin() + i * cols_L, L2.begin() + (i + 1) * cols_L);
			temp2.assign(H2.begin() + i * cols_H, H2.begin() + (i + 1) * cols_H);
			ilwt<T> iwt(temp1, temp2, lft);
			vector<T> sig;
			iwt.getSignal(sig);
			oup.insert(oup.end(), sig.begin(), sig.end());
			if (i == 0)
			{
				cx = (int)sig.size();
			}
		}

		signal = oup;
		oup_row = rows_L;
		oup_col = cx;
	}

	/* ilwt2(lwt2<T> &wt,liftscheme &lft) {
		 vector<T> A,H,V,D;
		 wt.getCoef(A,H,V,D);
		 vector<int> length;
		 wt.getDim(length);

		 int cols_L=length[1];

		 int rows_LL=length[0];
		 int rows_LH=length[2];

		 vector<T> AT,HT,VT,DT;

		 transpose(A,length[0],length[1],AT);
		 transpose(H,length[2],length[3],HT);
		 transpose(V,length[4],length[5],VT);
		 transpose(D,length[6],length[7],DT);
		 vector<T> L;
		 int rows_L;

		 for (int i=0; i < cols_L; i++) {
			 vector<T> temp1,temp2;
			 temp1.assign(AT.begin()+i*rows_LL,AT.begin()+(i+1)*rows_LL);
			 temp2.assign(HT.begin()+i*rows_LH,HT.begin()+(i+1)*rows_LH);
			 ilwt<T> iwt(temp1,temp2,lft);
			 vector<T> sig;
			 iwt.getSignal(sig);
			 L.insert(L.end(),sig.begin(),sig.end());
			 if (i==0) {
				 rows_L=(int) sig.size();
			 }


		 }

		 vector<T> L2;
		 transpose(L,cols_L,rows_L,L2);

		 int cols_H=length[5];

		 int rows_HL=length[4];
		 int rows_HH=length[6];
		 vector<T> H1;
		 int rows_H;

		 for (int i=0; i < cols_H; i++) {
			 vector<T> temp1,temp2;
			 temp1.assign(VT.begin()+i*rows_HL,VT.begin()+(i+1)*rows_HL);
			 temp2.assign(DT.begin()+i*rows_HH,DT.begin()+(i+1)*rows_HH);
			 ilwt<T> iwt(temp1,temp2,lft);
			 vector<T> sig;
			 iwt.getSignal(sig);
			 H1.insert(H1.end(),sig.begin(),sig.end());
			 if (i==0) {
				 rows_H=(int) sig.size();
			 }


		 }

		 vector<T> H2;
		 transpose(H1,cols_H,rows_H,H2);

		 vector<T> oup;
		 int cx;

		 for (int i=0; i < rows_L; i++) {
			 vector<T> temp1,temp2;
			 temp1.assign(L2.begin()+i*cols_L,L2.begin()+(i+1)*cols_L);
			 temp2.assign(H2.begin()+i*cols_H,H2.begin()+(i+1)*cols_H);
			 ilwt<T> iwt(temp1,temp2,lft);
			 vector<T> sig;
			 iwt.getSignal(sig);
			 oup.insert(oup.end(),sig.begin(),sig.end());
			 if (i==0) {
				 cx=(int) sig.size();
			 }


		 }

		 signal=oup;
		 oup_row=rows_L;
		 oup_col=cx;


	 }	*/

	ilwt2(lwt2<T> &wt, liftscheme &lft)
	{
		int J = wt.getLevels();
		vector<T> A1, B1, C1, D1;
		wt.getCoef(A1, B1, C1, D1);
		vector<int> len_coef;
		wt.getDim(len_coef);
		int slevel;

		int count = 0;
		vector<T> A, H, V, D;
		A = A1;
		vector<int> rxcx;

		for (int i = 0; i < J; i++)
		{
			vector<int> temp_coef, detlenH, detlenV, detlenD;
			slevel = J - i;
			if (i == 0)
			{
				temp_coef.assign(len_coef.begin(), len_coef.begin() + 8);
				count = count + 8;
			}
			else
			{
				temp_coef.assign(len_coef.begin() + count, len_coef.begin() + count + 6);
				temp_coef.insert(temp_coef.begin(), rxcx.begin(), rxcx.end());
				count = count + 6;
			}
			// Get LH,HL and HH coefficients
			getDetails(wt, "LH", slevel, H, detlenH);
			getDetails(wt, "HL", slevel, V, detlenV);
			getDetails(wt, "HH", slevel, D, detlenD);

			ilwt2<T> iwt2(A, H, V, D, temp_coef, lft);
			vector<T> temp;
			iwt2.getSignal(temp);
			rxcx.clear();
			iwt2.getDim(rxcx);
			oup_row = rxcx[0];
			oup_col = rxcx[1];
			A = temp;
		}
		signal = A;
	}

	void getDetails(lwt2<T> &wt, string align, int slevel, vector<T> &det_vec, vector<int> det_len)
	{
		// slevel 是现在读取的层数
		int J = wt.getLevels();
		int lev;
		if (slevel > J) // 如果大于总的层数
		{
			cout << " Decomposition has only " << J << " levels" << endl;
			exit(1);
		}
		else
		{
			lev = J - slevel; // 相当i?
		}

		vector<int> sig_vec;
		wt.getDim(sig_vec);
		vector<T> A1, B1, C1, D1;
		wt.getCoef(A1, B1, C1, D1);
		int total = 0;

		if (align == "LH" || align == "lh")
		{

			for (int i = 0; i < lev; i++)
			{
				total = total + sig_vec[2 + i * 6] * sig_vec[3 + i * 6];
			}
			det_vec.assign(B1.begin() + total, B1.begin() + total + sig_vec[2 + lev * 6] * sig_vec[3 + lev * 6]);
			det_len.push_back(sig_vec[2 + lev * 6]);
			det_len.push_back(sig_vec[3 + lev * 6]);
		}
		else if (align == "HL" || align == "hl")
		{

			for (int i = 0; i < lev; i++)
			{
				total = total + sig_vec[4 + i * 6] * sig_vec[5 + i * 6];
			}
			det_vec.assign(C1.begin() + total, C1.begin() + total + sig_vec[4 + lev * 6] * sig_vec[5 + lev * 6]);
			det_len.push_back(sig_vec[4 + lev * 6]);
			det_len.push_back(sig_vec[5 + lev * 6]);
		}
		else if (align == "HH" || align == "hh")
		{

			for (int i = 0; i < lev; i++)
			{
				total = total + sig_vec[6 + i * 6] * sig_vec[7 + i * 6];
			}
			det_vec.assign(D1.begin() + total, D1.begin() + total + sig_vec[6 + lev * 6] * sig_vec[7 + lev * 6]);
			det_len.push_back(sig_vec[6 + lev * 6]);
			det_len.push_back(sig_vec[7 + lev * 6]);
		}
		else
		{
			cout << "Accepted filter stages are LH or lh, HL or hl and HH or hh" << endl;
			exit(1);
		}
	}

	void getSignal(vector<T> &ilwt_oup)
	{
		ilwt_oup = signal;
	}

	void getDim(vector<int> &sigdim)
	{
		sigdim.push_back(oup_row);
		sigdim.push_back(oup_col);
	}

	virtual ~ilwt2()
	{
	}
};

template <typename T>
class lwt3
{
	vector<T> cLLL, cLLH, cLHL, cLHH, cHLL, cHLH, cHHL, cHHH;
	int widthLLL, lengthLLL, depthLLL;
	int widthLLH, lengthLLH, depthLLH;
	int widthLHL, lengthLHL, depthLHL;
	int widthLHH, lengthLHH, depthLHH;
	int widthHLL, lengthHLL, depthHLL;
	int widthHLH, lengthHLH, depthHLH;
	int widthHHL, lengthHHL, depthHHL;
	int widthHHH, lengthHHH, depthHHH;
	int level;
	vector<int> coef_lengths;

public:
	lwt3(vector<T> &signal, int width, int length, int depth, liftscheme &lft)
	{
		// 第一次 先x方向
		vector<T> L, H;
		int width_L, length_L, depth_L;
		int width_H, length_H, depth_H;
		// width减半 其他不变
		length_L = length;
		length_H = length;
		depth_L = depth;
		depth_H = depth;
		// rows_L = rows;
		// rows_H = rows;
		for (int i = 0; i < length * depth; i++)
		{
			// 本次处理的signal
			vector<T> temp; //{172, 47}
			temp.assign(signal.begin() + i * width, signal.begin() + (i + 1) * width);
			// 一维dwt
			lwt<T> lwt1(temp, lft);
			// 得到系数
			vector<T> a, d; // 154
			lwt1.getCoeff(a, d);
			// 将向量 ab 中的所有元素插入到向量 L H 的末尾
			L.insert(L.end(), a.begin(), a.end());
			H.insert(H.end(), d.begin(), d.end());
			if (i == 0)
			{
				width_L = a.size();
				width_H = d.size();
			}
		}

		// 第二次 对y进行
		//  转置操作
		vector<T> LT, HT;
		transpose3_xy(L, width_L, length_L, depth_L, LT);
		transpose3_xy(H, width_H, length_H, depth_H, HT);

		// 对原本的L 计算H L操作  Low Pass Stage
		int width_LL, length_LL, depth_LL;
		int width_LH, length_LH, depth_LH;
		// length改变 其余不变
		width_LL = width_L;
		depth_LL = depth_L;
		width_LH = width_L;
		depth_LH = depth_L;
		// Remove cout
		vector<T> LL, LH;
		for (int i = 0; i < width_L * depth_L; i++)
		{
			vector<T> temp;
			// LT:1*2*4
			/*{154.85638507985391, 218.49599538664319,
			224.85995641732214, 210.71782079359119,
			155.56349186104046, 185.96908345206202,
			 86.974134085945352,  202.23253941935261}
			*/
			temp.assign(LT.begin() + i * length_L, LT.begin() + (i + 1) * length_L);
			/*temp:{154.85638507985391,  218.49599538664319} 第一行*/

			// 进行dwt
			lwt<T> lwt1(temp, lft);
			// 获取参数
			vector<T> a, d; // 264 45
			lwt1.getCoeff(a, d);
			LL.insert(LL.end(), a.begin(), a.end());
			LH.insert(LH.end(), d.begin(), d.end());
			if (i == 0)
			{
				length_LL = a.size();
				length_LH = d.size();
			}
		}

		// 在原本的H 上进行H L
		int width_HL, length_HL, depth_HL;
		int width_HH, length_HH, depth_HH;
		// length改变 其余不变
		width_HH = width_H;
		depth_HH = depth_H;
		width_HL = width_H;
		depth_HL = depth_H;
		vector<T> HL, HH;
		for (int i = 0; i < width_H * depth_H; i++)
		{
			vector<T> temp;
			temp.assign(HT.begin() + i * length_H, HT.begin() + (i + 1) * length_H);
			lwt<T> lwt1(temp, lft);
			vector<T> a, d;
			lwt1.getCoeff(a, d);
			HL.insert(HL.end(), a.begin(), a.end());
			HH.insert(HH.end(), d.begin(), d.end());
			if (i == 0)
			{
				length_HL = a.size();
				length_HH = d.size();
			}
		}
		// 转置回去
		vector<T> cLL, cLH, cHL, cHH;
		transpose3_xy(LL, length_LL, width_LL, depth_LL, cLL);
		transpose3_xy(LH, length_LH, width_LH, depth_LH, cLH);
		transpose3_xy(HL, length_HL, width_HL, depth_HL, cHL);
		transpose3_xy(HH, length_HH, width_HH, depth_HH, cHH);

		// 第三次 对z进行变换
		// 转置处理
		vector<T> cLLT, cLHT, cHLT, cHHT;
		transpose3_xz(cLL, width_LL, length_LL, depth_LL, cLLT);
		transpose3_xz(cLH, width_LH, length_LH, depth_LH, cLHT);
		transpose3_xz(cHL, width_HL, length_HL, depth_HL, cHLT);
		transpose3_xz(cHH, width_HH, length_HH, depth_HH, cHHT);

		// 对LL进行处理
		int width_LLL, length_LLL, depth_LLL;
		int width_LLH, length_LLH, depth_LLH;
		// depth改变 其余不变
		width_LLL = width_LL;
		length_LLL = length_LL;
		width_LLH = width_LL;
		length_LLH = length_LL;
		vector<T> LLL, LLH;
		for (int i = 0; i < width_LL * length_LL; i++)
		{
			// cllt:{264, 308.00000000000006,241.50000000000003, 204.50000000000003}
			/*p a
			$8 = std::vector of length 2, capacity 2 = {404.46507883870521,
			  315.36962440920024}
			(gdb) p d
			$9 = std::vector of length 2, capacity 2 = {31.112698372208129,
			  -26.162950903902257}
			*/
			vector<T> temp;
			temp.assign(cLLT.begin() + i * depth_LL, cLLT.begin() + (i + 1) * depth_LL);
			lwt<T> lwt1(temp, lft);
			vector<T> a, d;
			lwt1.getCoeff(a, d);
			LLL.insert(LLL.end(), a.begin(), a.end());
			LLH.insert(LLH.end(), d.begin(), d.end());
			if (i == 0)
			{
				depth_LLL = a.size();
				depth_LLH = d.size();
			}
		}

		// 对LH进行处理
		int width_LHL, length_LHL, depth_LHL;
		int width_LHH, length_LHH, depth_LHH;
		// depth改变 其余不变
		width_LHL = width_LH;
		length_LHL = length_LH;
		width_LHH = width_LH;
		length_LHH = length_LH;
		vector<T> LHL, LHH;
		for (int i = 0; i < width_LH * length_LH; i++)
		{
			vector<T> temp;
			temp.assign(cLHT.begin() + i * depth_LH, cLHT.begin() + (i + 1) * depth_LH);
			lwt<T> lwt1(temp, lft);
			vector<T> a, d;
			lwt1.getCoeff(a, d);
			LHL.insert(LHL.end(), a.begin(), a.end());
			LHH.insert(LHH.end(), d.begin(), d.end());
			if (i == 0)
			{
				depth_LHL = a.size();
				depth_LHH = d.size();
			}
		}

		// 对HL进行处理
		int width_HLL, length_HLL, depth_HLL;
		int width_HLH, length_HLH, depth_HLH;
		// depth改变 其余不变
		width_HLL = width_HL;
		length_HLL = length_HL;
		width_HLH = width_HL;
		length_HLH = length_HL;
		vector<T> HLL, HLH;
		for (int i = 0; i < width_HL * length_HL; i++)
		{
			vector<T> temp;
			temp.assign(cHLT.begin() + i * depth_HL, cHLT.begin() + (i + 1) * depth_HL);
			lwt<T> lwt1(temp, lft);
			vector<T> a, d;
			lwt1.getCoeff(a, d);
			HLL.insert(HLL.end(), a.begin(), a.end());
			HLH.insert(HLH.end(), d.begin(), d.end());
			if (i == 0)
			{
				depth_HLL = a.size();
				depth_HLH = d.size();
			}
		}

		// 对HH进行处理
		int width_HHL, length_HHL, depth_HHL;
		int width_HHH, length_HHH, depth_HHH;
		// depth改变 其余不变
		width_HHL = width_HH;
		length_HHL = length_HH;
		width_HHH = width_HH;
		length_HHH = length_HH;
		vector<T> HHL, HHH;
		for (int i = 0; i < width_HH * length_HH; i++)
		{
			vector<T> temp;
			temp.assign(cHHT.begin() + i * depth_HH, cHHT.begin() + (i + 1) * depth_HH);
			lwt<T> lwt1(temp, lft);
			vector<T> a, d;
			lwt1.getCoeff(a, d);
			HHL.insert(HHL.end(), a.begin(), a.end());
			HHH.insert(HHH.end(), d.begin(), d.end());
			if (i == 0)
			{
				depth_HHL = a.size();
				depth_HHH = d.size();
			}
		}

		// 转置回去
		transpose3_xz(LLL, depth_LLL, length_LLL, width_LLL, cLLL);
		transpose3_xz(LLH, depth_LLH, length_LLH, width_LLH, cLLH);
		transpose3_xz(LHL, depth_LHL, length_LHL, width_LHL, cLHL);
		transpose3_xz(LHH, depth_LHH, length_LHH, width_LHH, cLHH);
		transpose3_xz(HLL, depth_HLL, length_HLL, width_HLL, cHLL);
		transpose3_xz(HLH, depth_HLH, length_HLH, width_HLH, cHLH);
		transpose3_xz(HHL, depth_HHL, length_HHL, width_HHL, cHHL);
		transpose3_xz(HHH, depth_HHH, length_HHH, width_HHH, cHHH);
		level = 1;
		widthLLL = width_LLL;
		lengthLLL = length_LLL;
		depthLLL = depth_LLL;

		widthLLH = width_LLH;
		lengthLLH = length_LLH;
		depthLLH = depth_LLH;

		widthLHL = width_LHL;
		lengthLHL = length_LHL;
		depthLHL = depth_LHL;

		widthLHH = width_LHH;
		lengthLHH = length_LHH;
		depthLHH = depth_LHH;

		widthHLL = width_HLL;
		lengthHLL = length_HLL;
		depthHLL = depth_HLL;

		widthHLH = width_HLH;
		lengthHLH = length_HLH;
		depthHLH = depth_HLH;

		widthHHL = width_HHL;
		lengthHHL = length_HHL;
		depthHHL = depth_HHL;

		widthHHH = width_HHH;
		lengthHHH = length_HHH;
		depthHHH = depth_HHH;
		// 依次将每个变量放入 coef_lengths
		coef_lengths.push_back(widthLLL);
		coef_lengths.push_back(lengthLLL);
		coef_lengths.push_back(depthLLL);

		coef_lengths.push_back(widthLLH);
		coef_lengths.push_back(lengthLLH);
		coef_lengths.push_back(depthLLH);

		coef_lengths.push_back(widthLHL);
		coef_lengths.push_back(lengthLHL);
		coef_lengths.push_back(depthLHL);

		coef_lengths.push_back(widthLHH);
		coef_lengths.push_back(lengthLHH);
		coef_lengths.push_back(depthLHH);

		coef_lengths.push_back(widthHLL);
		coef_lengths.push_back(lengthHLL);
		coef_lengths.push_back(depthHLL);

		coef_lengths.push_back(widthHLH);
		coef_lengths.push_back(lengthHLH);
		coef_lengths.push_back(depthHLH);

		coef_lengths.push_back(widthHHL);
		coef_lengths.push_back(lengthHHL);
		coef_lengths.push_back(depthHHL);

		coef_lengths.push_back(widthHHH);
		coef_lengths.push_back(lengthHHH);
		coef_lengths.push_back(depthHHH);
	}
	lwt3(const vector<T> &signal, int width, int length, int depth, liftscheme &lft, int J)
	{

		vector<T> LLL, LLH, LHL, LHH, HLL, HLH, HHL, HHH;
		vector<int> siglen;
		vector<T> temp_signal = signal; // 存储中的L结果
		for (int i = 0; i < J; i++)
		{
			lwt3<T> wt3(temp_signal, width, length, depth, lft);
			vector<T> tempA, tempB, tempC, tempD, tempE, tempF, tempG, tempH;
			vector<int> temp_siglen;
			wt3.getCoef(tempA, tempB, tempC, tempD, tempE, tempF, tempG, tempH);
			temp_signal = tempA;
			LLL = tempA;
			LLH.insert(LLH.begin(), tempB.begin(), tempB.end()); // 每一次LLH的结果依次存储
			LHL.insert(LHL.begin(), tempC.begin(), tempC.end());
			LHH.insert(LHH.begin(), tempD.begin(), tempD.end());
			HLL.insert(HLL.begin(), tempE.begin(), tempE.end());
			HLH.insert(HLH.begin(), tempF.begin(), tempF.end());
			HHL.insert(HHL.begin(), tempG.begin(), tempG.end());
			HHH.insert(HHH.begin(), tempH.begin(), tempH.end());
			wt3.getDim(temp_siglen);
			width = temp_siglen[0];
			length = temp_siglen[1];
			depth = temp_siglen[2];
			if (i == J - 1) // 如果是最后一次
			{
				siglen.insert(siglen.begin(), temp_siglen.begin(), temp_siglen.end());
			}
			else // 不要考LLL
			{
				siglen.insert(siglen.begin(), temp_siglen.begin() + 3, temp_siglen.end());
			}
		}
		cLLL = LLL;
		cLLH = LLH;
		cLHL = LHL;
		cLHH = LHH;
		cHLL = HLL;
		cHLH = HLH;
		cHHL = HHL;
		cHHH = HHH;
		level = J;
		coef_lengths = siglen;
	}
	void getDetails_LL1(vector<T> &det_vec)
	{
		vector<T> A1;
		A1 = cLLL;
		det_vec.assign(A1.begin(), A1.begin() + coef_lengths[0] * coef_lengths[1] * coef_lengths[2]);
	}

	void getCoef(vector<T> &aLLL, vector<T> &aLLH, vector<T> &aLHL, vector<T> &aLHH,
				 vector<T> &aHLL, vector<T> &aHLH, vector<T> &aHHL, vector<T> &aHHH)
	{
		aLLL = cLLL;
		aLLH = cLLH;
		aLHL = cLHL;
		aLHH = cLHH;
		aHLL = cHLL;
		aHLH = cHLH;
		aHHL = cHHL;
		aHHH = cHHH;
	}

	void getLLL(vector<T> &aLLL)
	{
		aLLL = cLLL;
	}
	void getHHH(vector<T> &aHHH)
	{
		aHHH = cHHH;
	}
	void getLLH(vector<T> &aLLH)
	{
		aLLH = cLLH;
	}
	void getDim(vector<int> &dimvec)
	{

		dimvec = coef_lengths;
	}
	int getLevels()
	{
		return level;
	}
	virtual ~lwt3()
	{
	}
};

template <typename T>
class ilwt3
{
	vector<T> signal;
	int oup_width, oup_length, oup_depth;

public:
	ilwt3(vector<T> &LLL, vector<T> &LLH, vector<T> &LHL, vector<T> &LHH, vector<T> &HLL, vector<T> &HLH, vector<T> &HHL, vector<T> &HHH, vector<int> &length, liftscheme &lft)
	{
		// 先对z
		vector<T> cLLL, cLLH, cLHL, cLHH, cHLL, cHLH, cHHL, cHHH;
		// 转置 把z换到最前面
		transpose3_xz(LLL, length[0], length[1], length[2], cLLL);
		transpose3_xz(LLH, length[3], length[4], length[5], cLLH);

		transpose3_xz(LHL, length[6], length[7], length[8], cLHL);
		transpose3_xz(LHH, length[9], length[10], length[11], cLHH);

		transpose3_xz(HLL, length[12], length[13], length[14], cHLL);
		transpose3_xz(HLH, length[15], length[16], length[17], cHLH);

		transpose3_xz(HHL, length[18], length[19], length[20], cHHL);
		transpose3_xz(HHH, length[21], length[22], length[23], cHHH);

		// 计算LL
		vector<T> LL;
		int depth_LL;
		for (int i = 0; i < length[0] * length[1]; i++)
		{
			vector<T> temp1, temp2;
			temp1.assign(cLLL.begin() + i * length[2], cLLL.begin() + (i + 1) * length[2]);
			temp2.assign(cLLH.begin() + i * length[5], cLLH.begin() + (i + 1) * length[5]);
			ilwt<T> iwt(temp1, temp2, lft);
			vector<T> sig;
			iwt.getSignal(sig);
			LL.insert(LL.end(), sig.begin(), sig.end());
			if (i == 0)
			{
				depth_LL = (int)sig.size();
			}
		}
		// 计算LH
		vector<T> LH;
		int depth_LH;
		for (int i = 0; i < length[6] * length[7]; i++)
		{
			vector<T> temp1, temp2;
			temp1.assign(cLHL.begin() + i * length[8], cLHL.begin() + (i + 1) * length[8]);
			temp2.assign(cLHH.begin() + i * length[11], cLHH.begin() + (i + 1) * length[11]);
			ilwt<T> iwt(temp1, temp2, lft);
			vector<T> sig;
			iwt.getSignal(sig);
			LH.insert(LH.end(), sig.begin(), sig.end());
			if (i == 0)
			{
				depth_LH = (int)sig.size();
			}
		}

		// HL
		vector<T> HL;
		int depth_HL;
		for (int i = 0; i < length[12] * length[13]; i++)
		{
			vector<T> temp1, temp2;
			temp1.assign(cHLL.begin() + i * length[14], cHLL.begin() + (i + 1) * length[14]);
			temp2.assign(cHLH.begin() + i * length[17], cHLH.begin() + (i + 1) * length[17]);
			ilwt<T> iwt(temp1, temp2, lft);
			vector<T> sig;
			iwt.getSignal(sig);
			HL.insert(HL.end(), sig.begin(), sig.end());
			if (i == 0)
			{
				depth_HL = (int)sig.size();
			}
		}

		// HH
		vector<T> HH;
		int depth_HH;
		for (int i = 0; i < length[18] * length[19]; i++)
		{
			vector<T> temp1, temp2;
			temp1.assign(cHHL.begin() + i * length[20], cHHL.begin() + (i + 1) * length[20]);
			temp2.assign(cHHH.begin() + i * length[23], cHHH.begin() + (i + 1) * length[23]);
			ilwt<T> iwt(temp1, temp2, lft);
			vector<T> sig;
			iwt.getSignal(sig);
			HH.insert(HH.end(), sig.begin(), sig.end());
			if (i == 0)
			{
				depth_HH = (int)sig.size();
			}
		}

		// 转置回去
		vector<T> cLL, cLH, cHL, cHH;
		transpose3_xz(LL, depth_LL, length[1], length[0], cLL);
		transpose3_xz(LH, depth_LH, length[7], length[6], cLH);

		transpose3_xz(HL, depth_HL, length[13], length[12], cHL);
		transpose3_xz(HH, depth_HH, length[19], length[18], cHH);

		// 在对y

		// 转置xy
		vector<T> cLLT, cLHT, cHLT, cHHT;
		transpose3_xy(cLL, length[0], length[1], depth_LL, cLLT);
		transpose3_xy(cLH, length[6], length[7], depth_LH, cLHT);
		transpose3_xy(cHL, length[12], length[13], depth_HL, cHLT);
		transpose3_xy(cHH, length[18], length[19], depth_HH, cHHT);

		// 计算LL LH ->L
		vector<T> L;
		int length_L;
		for (int i = 0; i < length[0] * depth_LL; i++)
		{
			vector<T> temp1, temp2;
			temp1.assign(cLLT.begin() + i * length[1], cLLT.begin() + (i + 1) * length[1]);
			temp2.assign(cLHT.begin() + i * length[7], cLHT.begin() + (i + 1) * length[7]);
			ilwt<T> iwt(temp1, temp2, lft);
			vector<T> sig;
			iwt.getSignal(sig);
			L.insert(L.end(), sig.begin(), sig.end());
			if (i == 0)
			{
				length_L = (int)sig.size();
			}
		}

		// 计算HL HH->H
		vector<T> H;
		int length_H;
		for (int i = 0; i < length[12] * depth_HL; i++)
		{
			vector<T> temp1, temp2;
			temp1.assign(cHLT.begin() + i * length[13], cHLT.begin() + (i + 1) * length[13]);
			temp2.assign(cHHT.begin() + i * length[19], cHHT.begin() + (i + 1) * length[19]);
			ilwt<T> iwt(temp1, temp2, lft);
			vector<T> sig;
			iwt.getSignal(sig);
			H.insert(H.end(), sig.begin(), sig.end());
			if (i == 0)
			{
				length_H = (int)sig.size();
			}
		}

		// 转置回去
		vector<T> LT, HT;
		transpose3_xy(L, length_L, length[0], depth_LL, LT);
		transpose3_xy(H, length_H, length[12], depth_LL, HT);

		// 对x L H -> x
		vector<T> oup;
		int width_oup;
		for (int i = 0; i < length_L * depth_LL; i++)
		{
			vector<T> temp1, temp2;
			temp1.assign(LT.begin() + i * length[0], LT.begin() + (i + 1) * length[0]);
			temp2.assign(HT.begin() + i * length[12], HT.begin() + (i + 1) * length[12]);
			ilwt<T> iwt(temp1, temp2, lft);
			vector<T> sig;
			iwt.getSignal(sig);
			oup.insert(oup.end(), sig.begin(), sig.end());
			if (i == 0)
			{
				width_oup = (int)sig.size();
			}
		}

		signal = oup;
		oup_width = width_oup;
		oup_length = length_L;
		oup_depth = depth_LL;
	}

	ilwt3(lwt3<T> &wt, liftscheme &lft)
	{
		int J = wt.getLevels();
		vector<T> LLL, LLH, LHL, LHH, HLL, HLH, HHL, HHH;
		wt.getCoef(LLL, LLH, LHL, LHH, HLL, HLH, HHL, HHH);
		vector<int> len_coef;
		wt.getDim(len_coef);

		vector<T> tempA, tempB, tempC, tempD, tempE, tempF, tempG, tempH;
		tempA = LLL;

		int slevel;
		int count = 0;

		vector<int> length_tempA, length_tempB, length_tempC, length_tempD, length_tempE, length_tempF, length_tempG, length_tempH;
		vector<int> rxcx; // 反变换之后的长度
		vector<T> temp;
		for (int i = 0; i < J; i++)
		{
			vector<int> temp_coef;
			slevel = J - i;
			if (i == 0) // 如果是第一次 复制全部的3*8个
			{
				temp_coef.assign(len_coef.begin(), len_coef.begin() + 24);
				count = count + 24;
			}
			else // 如果不是第一次 复制3*7个
			{
				temp_coef.assign(len_coef.begin() + count, len_coef.begin() + count + 21);
				temp_coef.insert(temp_coef.begin(), rxcx.begin(), rxcx.end());
				count = count + 21;
			}

			// getDetails(wt, "LLL", slevel, tempA, length_tempA);
			getDetails(wt, "LLH", slevel, tempB, length_tempB);
			getDetails(wt, "LHL", slevel, tempC, length_tempC);
			getDetails(wt, "LHH", slevel, tempD, length_tempD);
			getDetails(wt, "HLL", slevel, tempE, length_tempE);
			getDetails(wt, "HLH", slevel, tempF, length_tempF);
			getDetails(wt, "HHL", slevel, tempG, length_tempG);
			getDetails(wt, "HHH", slevel, tempH, length_tempH);
			ilwt3<T> iwt3(tempA, tempB, tempC, tempD, tempE, tempF, tempG, tempH, temp_coef, lft);

			iwt3.getSignal(temp);
			rxcx.clear();
			iwt3.getDim(rxcx);

			oup_width = rxcx[0];
			oup_length = rxcx[1];
			oup_depth = rxcx[2];
			tempA = temp;
		}
		signal = temp;
	}
	// 只通过LLL进行重建 应该是错误的
	//  ilwt3(vector<T> &LLL,  vector<int> &len_coef, liftscheme &lft, int J)
	//  {
	//  	vector<T> tempA, tempB, tempC, tempD, tempE, tempF, tempG, tempH;
	//  	tempA = LLL;
	//  	int slevel;
	//  	int count = 0;

	// 	vector<int> length_tempA, length_tempB, length_tempC, length_tempD, length_tempE, length_tempF, length_tempG, length_tempH;
	// 	vector<int> rxcx; // 反变换之后的长度
	// 	vector<T> temp;
	// 	for (int i = 0; i < J; i++)
	// 	{
	// 		vector<int> temp_coef;
	// 		slevel = J - i;
	// 		if (i == 0) // 如果是第一次 复制全部的3*8个
	// 		{
	// 			temp_coef.assign(len_coef.begin(), len_coef.begin() + 24);
	// 			count = count + 24;
	// 		}
	// 		else // 如果不是第一次 复制3*7个
	// 		{
	// 			temp_coef.assign(len_coef.begin() + count, len_coef.begin() + count + 21);
	// 			temp_coef.insert(temp_coef.begin(), rxcx.begin(), rxcx.end());
	// 			count = count + 21;
	// 		}
	// 		// getDetails(wt, "LLH", slevel, tempB, length_tempB);
	// 		// getDetails(wt, "LHL", slevel, tempC, length_tempC);
	// 		// getDetails(wt, "LHH", slevel, tempD, length_tempD);
	// 		// getDetails(wt, "HLL", slevel, tempE, length_tempE);
	// 		// getDetails(wt, "HLH", slevel, tempF, length_tempF);
	// 		// getDetails(wt, "HHL", slevel, tempG, length_tempG);
	// 		// getDetails(wt, "HHH", slevel, tempH, length_tempH);
	// 		// std::vector<float> tempB(tempA.size(), 0);
	// 		std::vector<float> tempB(tempA.size(), 0);
	// 		std::vector<float> tempC(tempA.size(), 0);
	// 		std::vector<float> tempD(tempA.size(), 0);
	// 		std::vector<float> tempE(tempA.size(), 0);
	// 		std::vector<float> tempF(tempA.size(), 0);
	// 		std::vector<float> tempG(tempA.size(), 0);
	// 		std::vector<float> tempH(tempA.size(), 0);

	// 		ilwt3<T> iwt3(tempA, tempB, tempC, tempD, tempE, tempF, tempG, tempH, temp_coef, lft);
	// 		iwt3.getSignal(temp);
	// 		rxcx.clear();
	// 		iwt3.getDim(rxcx);

	// 		oup_width = rxcx[0];
	// 		oup_length = rxcx[1];
	// 		oup_depth = rxcx[2];
	// 		tempA = temp;
	// 	}
	// 	signal = temp;
	// }
	ilwt3(vector<T> &LLL, vector<T> &LLH, vector<T> &LHL, vector<T> &LHH,
		  vector<T> &HLL, vector<T> &HLH, vector<T> &HHL, vector<T> &HHH, 
		  const vector<int> &len_coef, liftscheme &lft, int J)
	{
		vector<T> tempA, tempB, tempC, tempD, tempE, tempF, tempG, tempH;
		tempA = LLL;
		tempB = LLH;
		tempC = LHL;
		tempD = LHH;
		tempE = HLL;
		tempF = HLH;
		tempG = HHL;
		tempH = HHH;
		int slevel;
		int count = 0;

		vector<int> rxcx; // 反变换之后的长度
		vector<T> temp;
		for (int i = 0; i < J; i++)
		{
			vector<int> temp_coef;
			slevel = J - i;
			if (i == 0) // 如果是第一次 复制全部的3*8个
			{
				temp_coef.assign(len_coef.begin(), len_coef.begin() + 24);
				count = count + 24;
			}
			else // 如果不是第一次 复制3*7个
			{
				temp_coef.assign(len_coef.begin() + count, len_coef.begin() + count + 21);
				temp_coef.insert(temp_coef.begin(), rxcx.begin(), rxcx.end());
				count = count + 21;
			}

			ilwt3<T> iwt3(tempA, tempB, tempC, tempD, tempE, tempF, tempG, tempH, temp_coef, lft);
			iwt3.getSignal(temp);
			rxcx.clear();
			iwt3.getDim(rxcx);

			oup_width = rxcx[0];
			oup_length = rxcx[1];
			oup_depth = rxcx[2];
			tempA = temp;
		}
		signal = temp;
	}
	void getDetails(lwt3<T> &wt, string align, int slevel, vector<T> &det_vec, vector<int> &det_len)
	{
		// slevel 是现在读取的层数
		int J = wt.getLevels();
		int lev;
		if (slevel > J) // 如果大于总的层数
		{
			cout << " Decomposition has only " << J << " levels" << endl;
			exit(1);
		}
		else
		{
			lev = J - slevel; // 相当i?
		}
		vector<int> sig_vec;
		wt.getDim(sig_vec);

		vector<T> tempA, tempB, tempC, tempD, tempE, tempF, tempG, tempH;
		wt.getCoef(tempA, tempB, tempC, tempD, tempE, tempF, tempG, tempH);
		int total = 0;

		if (align == "LLH" || align == "llh")
		{
			for (int i = 0; i < lev; i++) // 差14
			{
				total = total + sig_vec[3 + i * 21] * sig_vec[4 + i * 21] * sig_vec[5 + i * 21];
			}
			det_vec.assign(tempB.begin() + total, tempB.begin() + total + sig_vec[3 + lev * 21] * sig_vec[4 + lev * 21] * sig_vec[5 + lev * 21]);
			det_len.push_back(sig_vec[3 + lev * 21]);
			det_len.push_back(sig_vec[4 + lev * 21]);
			det_len.push_back(sig_vec[5 + lev * 21]);
		}
		else if (align == "LHL" || align == "lhl")
		{

			for (int i = 0; i < lev; i++)
			{
				total = total + sig_vec[6 + i * 21] * sig_vec[7 + i * 21] * sig_vec[8 + i * 21];
			}
			det_vec.assign(tempC.begin() + total, tempC.begin() + total + sig_vec[6 + lev * 21] * sig_vec[7 + lev * 21] * sig_vec[8 + lev * 21]);
			det_len.push_back(sig_vec[6 + lev * 21]);
			det_len.push_back(sig_vec[7 + lev * 21]);
			det_len.push_back(sig_vec[8 + lev * 21]);
		}
		else if (align == "LHH" || align == "lhh")
		{

			for (int i = 0; i < lev; i++)
			{
				total = total + sig_vec[9 + i * 21] * sig_vec[10 + i * 21] * sig_vec[11 + i * 21];
			}
			det_vec.assign(tempD.begin() + total, tempD.begin() + total + sig_vec[9 + lev * 21] * sig_vec[10 + lev * 21] * sig_vec[11 + lev * 21]);
			det_len.push_back(sig_vec[9 + lev * 21]);
			det_len.push_back(sig_vec[10 + lev * 21]);
			det_len.push_back(sig_vec[11 + lev * 21]);
		}
		else if (align == "HLL" || align == "hll")
		{

			for (int i = 0; i < lev; i++)
			{
				total = total + sig_vec[12 + i * 21] * sig_vec[13 + i * 21] * sig_vec[14 + i * 21];
			}
			det_vec.assign(tempE.begin() + total, tempE.begin() + total + sig_vec[12 + lev * 21] * sig_vec[13 + lev * 21] * sig_vec[14 + lev * 21]);
			det_len.push_back(sig_vec[12 + lev * 21]);
			det_len.push_back(sig_vec[13 + lev * 21]);
			det_len.push_back(sig_vec[14 + lev * 21]);
		}
		else if (align == "HLH" || align == "hlh")
		{

			for (int i = 0; i < lev; i++)
			{
				total = total + sig_vec[15 + i * 21] * sig_vec[16 + i * 21] * sig_vec[17 + i * 21];
			}
			det_vec.assign(tempF.begin() + total, tempF.begin() + total + sig_vec[15 + lev * 21] * sig_vec[16 + lev * 21] * sig_vec[17 + lev * 21]);
			det_len.push_back(sig_vec[15 + lev * 21]);
			det_len.push_back(sig_vec[16 + lev * 21]);
			det_len.push_back(sig_vec[17 + lev * 21]);
		}
		else if (align == "HHL" || align == "hhl")
		{

			for (int i = 0; i < lev; i++)
			{
				total = total + sig_vec[18 + i * 21] * sig_vec[19 + i * 21] * sig_vec[20 + i * 21];
			}
			det_vec.assign(tempG.begin() + total, tempG.begin() + total + sig_vec[18 + lev * 21] * sig_vec[19 + lev * 21] * sig_vec[20 + lev * 21]);
			det_len.push_back(sig_vec[18 + lev * 21]);
			det_len.push_back(sig_vec[19 + lev * 21]);
			det_len.push_back(sig_vec[20 + lev * 21]);
		}
		else if (align == "HHH" || align == "hhh")
		{

			for (int i = 0; i < lev; i++)
			{
				total = total + sig_vec[21 + i * 21] * sig_vec[22 + i * 21] * sig_vec[23 + i * 21];
			}
			det_vec.assign(tempG.begin() + total, tempG.begin() + total + sig_vec[21 + lev * 21] * sig_vec[22 + lev * 21] * sig_vec[23 + lev * 21]);
			det_len.push_back(sig_vec[21 + lev * 21]);
			det_len.push_back(sig_vec[22 + lev * 21]);
			det_len.push_back(sig_vec[23 + lev * 21]);
		}
		else
		{
			cout << "Accepted filter stages are LH or lh, HL or hl and HH or hh" << endl;
			exit(1);
		}
	}
	void getSignal(vector<T> &ilwt_oup)
	{
		ilwt_oup = signal;
	}
	void getDim(vector<int> &sigdim)
	{
		sigdim.push_back(oup_width);
		sigdim.push_back(oup_length);
		sigdim.push_back(oup_depth);
	}
	virtual ~ilwt3()
	{
	}
};

#endif
