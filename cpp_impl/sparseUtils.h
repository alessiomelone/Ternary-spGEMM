#pragma once
#include <iostream>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <random>

template <typename T>
std::vector<T> initX(int LEN, int Range)
{
	std::vector<T> X(LEN, 0);
	std::mt19937 generator(static_cast<unsigned int>(time(0)));
	std::uniform_int_distribution<int> range(-Range, Range);
	for (int i = 0; i < LEN; i++)
	{
		X[i] = range(generator);
	}
	return X;
};

template <typename T>
std::vector<T> generateSparseMatrix(int H, int W, int nonZero, bool uniformDistribution)
{
	std::vector<T> y = std::vector<T>(H * W, 0);
	if (uniformDistribution)
	{
		for (int h = 0; h < H; h++)
		{
			for (int w = 0; w < W; w += nonZero * 2)
			{
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
		std::mt19937 generator(static_cast<unsigned int>(time(0)));
		std::uniform_int_distribution<int> range(0, W - 1);
		std::uniform_int_distribution<int> variRange(0, int(W / nonZero / 20 + 1));
		for (int h = 0; h < H; h++)
		{
			int posVari = variRange(generator);
			int limitPos = (W / nonZero) / 2 + posVari;
			int limitNeg = (W / nonZero) / 2 - posVari;

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

template <typename T>
void GEMM(T *X, T *W, T *b, T *Y, int M, int N, int K)
{
#pragma omp parallel for
	for (int m = 0; m < M; m++)
	{
		for (int n = 0; n < N; n++)
		{
			T y = 0;
			for (int k = 0; k < K; k++)
			{
				y += X[m * K + k] * W[k * N + n];
			}
			Y[m * N + n] = y + b[n];
		}
	}
}

template <typename T>
bool compare_results(T *result, T *groundTruth, int H, int W)
{
	for (int h = 0; h < H; h++)
	{
		for (int w = 0; w < W; w++)
		{
			int i = h * W + w;
			if (std::abs(result[i] - groundTruth[i]) > 10e-6)
			{
				std::cout << "Error at: H=" << h << ", W=" << w << ", result=" << result[i] << ", groundTruth=" << groundTruth[i] << std::endl;
				return false;
			}
		}
	}

	return true;
}