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

#ifndef ALG_H
#define ALG_H
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

// 用于将一个信号（sig）分割成偶数索引的部分（even）和奇数索引的部分（odd
template <typename T>
void inline split(vector<T> &sig, vector<T> &even, vector<T> &odd)
{

	for (int i = 0; i < (int)sig.size(); i++)
	{
		if (i % 2 == 0)
		{ // 偶数索引
			even.push_back(sig[i]);
		}
		else
		{ // 奇数索引
			odd.push_back(sig[i]);
		}
	}
}

// 第 i 个元素乘以 x，然后将结果存回相同的位置 归一化？
template <typename T>
void inline vecmult(vector<T> &sig, double x)
{

	for (int i = 0; i < (int)sig.size(); i++)
	{
		sig[i] = (T)(x * sig[i]);
	}
}

template <typename T>
void inline merge(vector<T> &sig, vector<T> &even, vector<T> &odd)
{

	int N = even.size() + odd.size();

	for (int i = 0; i < N; i++)
	{
		if (i % 2 == 0)
		{
			sig.push_back(even[i / 2]);
		}
		else
		{
			sig.push_back(odd[i / 2]);
		}
	}
}

// 将按行的二维数组变成按列的数组
template <typename T>
void inline transpose(vector<T> &sig, int rows, int cols, vector<T> &col)
{

	for (int i = 0; i < cols; i++)
	{ // 行
		for (int j = 0; j < rows; j++)
		{									  // 列
			col.push_back(sig[i + j * cols]); // 转置
		}
	}
}

template <typename T>
void inline transpose3_xy(vector<T> &sig, int width, int length, int depth, vector<T> &col)
{

	for (int z = 0; z < depth; z++)
	{
		for (int x = 0; x < width; x++)
		{
			for (int y = 0; y < length; y++)
			{
				// 有可能导致效率低
				col.push_back(sig[x + y * width + z * width * length]);
			}
		}
	}
}

template <typename T>
void inline transpose3_xz(vector<T> &sig, int width, int length, int depth, vector<T> &col)
{

	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < length; y++)
		{
			for (int z = 0; z < depth; z++)
			{
				col.push_back(sig[x + y * width + z * width * length]);
			}
		}
	}
}
#endif