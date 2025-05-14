#pragma once
#include <functional>

using comp_func = std::function<void(float *X, float *B, float *Y, int M, int N, int K)>;

double perf_test(comp_func f, int M, int K, int N, int nonZero);