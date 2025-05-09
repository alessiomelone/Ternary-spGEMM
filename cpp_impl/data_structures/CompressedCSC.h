#pragma once
#include <vector>
#include <iostream>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <random>
using namespace std;

class CompressedCSC
{
public:
	vector<int> col_start_pos;
	vector<int> col_start_neg;
	vector<int> row_index_pos;
	vector<int> row_index_neg;

	CompressedCSC(int *matrix, int K, int N)
	{
		int column_start_pos = 0;
		int column_start_neg = 0;
		for (int n = 0; n < N; n++)
		{
			this->col_start_pos.push_back(column_start_pos);
			this->col_start_neg.push_back(column_start_neg);
			for (int k = 0; k < K; k++)
			{
				if (matrix[k * N + n] >= 1)
				{
					column_start_pos++;
					this->row_index_pos.push_back(k);
				}
				else if (matrix[k * N + n] <= -1)
				{
					column_start_neg++;
					this->row_index_neg.push_back(k);
				}
			}
		}
		this->col_start_pos.push_back(column_start_pos);
		this->col_start_neg.push_back(column_start_neg);
	}
};

template <typename T>
vector<T> generateCompressedCSC(int H, int W, int nonZero, bool uniformDistribution)
{
	vector<T> y = vector<T>(H * W, 0);
	if (uniformDistribution)
	{
		for (int h = 0; h < H; h++)
		{
			for (int w = 0; w < W; w += nonZero * 2)
			{
				// Assign +1, -1 to each 2 x nonZero slots
				int randomA = rand() % nonZero * 2;
				int randomB = rand() % nonZero * 2;
				y[h * W + w + randomA] = 1;
				while (randomA == randomB)
				{
					randomB = rand() % nonZero * 2;
				}
				y[h * W + w + randomB] = -1;
			}
		}
	}
	else
	{
		mt19937 generator(static_cast<unsigned int>(time(0)));
		uniform_int_distribution<int> range(0, W - 1);
		uniform_int_distribution<int> variRange(0, int(W / nonZero / 20 + 1)); // The variation among different columns
		for (int h = 0; h < H; h++)
		{
			int posVari = variRange(generator);
			int limitPos = (W / nonZero) / 2 + posVari;
			int limitNeg = (W / nonZero) / 2 - posVari;

			// Assign +1 to W / nonZero / 2 places
			int count = 0;
			while (count < limitPos)
			{
				int randomA = range(generator);
				if (y[h * W + randomA] == 0)
				{
					y[h * W + randomA] = 1;
					count++;
				}
			}

			// Assign -1 to W / nonZero / 2 places
			count = 0;
			while (count < limitNeg)
			{
				int randomA = range(generator);
				if (y[h * W + randomA] == 0)
				{
					y[h * W + randomA] = -1;
					count++;
				}
			}
		}
	}

	return y;
}