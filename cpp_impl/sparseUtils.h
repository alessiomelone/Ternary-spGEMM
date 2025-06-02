#pragma once
#include <iostream>
#include <vector>
#include <random>
using namespace std;
template <typename T>
vector<T> initX(int LEN, int Range, bool uniformDistribution = false)
{
    vector<T> X(LEN, 0);
    mt19937 generator(static_cast<unsigned int>(time(0)));
    uniform_int_distribution<int> range(-Range, Range);
    if (uniformDistribution)
    {
        for (int i = 0; i < LEN; i++)
            X[i] = rand() % Range;
    }
    else
    {
        for (int i = 0; i < LEN; i++)
            X[i] = range(generator);
    }
    return X;
};

template <typename T>
vector<T> generateSparseMatrix(int H, int W, int nonZero, bool uniformDistribution, int seed = -1)
{
    if (seed != -1)
    {
        srand(seed);
    }

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
        mt19937 generator(static_cast<unsigned int>((seed == -1) ? time(0) : seed));
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
void GEMM_PreLU(T *X, T *W, T *b, T *alpha, T *Y, int M, int N, int K)
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
            T pre_activation = y + b[n];

            // Apply PReLU
            if (pre_activation >= 0)
            {
                Y[m * N + n] = pre_activation;
            }
            else
            {
                // Each output neuron 'n' has its own alpha value
                Y[m * N + n] = alpha[n] * pre_activation;
            }
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
            if (abs(result[i] - groundTruth[i]) > 10e-6)
            {
                cout << "Error at: H=" << h << ", W=" << w << ", result=" << result[i] << ", groundTruth=" << groundTruth[i] << endl;
                return false;
            }
        }
    }

    return true;
}