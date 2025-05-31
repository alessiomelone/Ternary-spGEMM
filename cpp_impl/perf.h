#include <functional> // <-- Add this

using comp_func = std::function<void(double *X, double *B, double *Y, int M, int N, int K)>;

// perf_test still generates data based on M, K, N, nonZero,
// but the comp_func 'f' it receives already has its specific sparse matrix data captured.
double perf_test(comp_func f, int M, int K, int N, int nonZero);