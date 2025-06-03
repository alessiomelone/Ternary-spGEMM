#include <functional>

using comp_func = std::function<void(float *X, float *B, float *Y, int M, int N, int K)>;
using comp_func_prelu = std::function<void(float *X, float *B, float *alpha, float *Y, int M, int N, int K)>;

float perf_test(comp_func f, int M, int K, int N, int nonZero);
float perf_test_prelu(comp_func_prelu f, int M, int K, int N, int nonZero);
// perf_test still generates data based on M, K, N, nonZero,
// but the comp_func 'f' it receives already has its specific sparse matrix data captured.