/*
 * debug.cpp
 *
 * Stand‑alone driver to debug a single sparse–dense matrix multiplication
 * implementation from comp.cpp.  Usage example:
 *
 *   ./debug -M 8 -K 8 -N 8 -s 20
 *
 * Flags
 *   -M      rows in X   (int > 0)
 *   -K      cols in X / rows in W   (int > 0)
 *   -N      cols in W   (int > 0)
 *   -s      number of non‑zero elements in W (int > 0)
 *
 * The program prints X, W, b, expected Y (FP32 GEMM) and the actual Y
 * returned by the chosen sparse kernel, then reports PASS/FAIL.
 */

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <functional>
#include <cstdlib>
#include <cstring>

#include "common.h"
#include "sparseUtils.h"

// --- Prototypes for implementations defined & explicitly instantiated in comp.cpp ---
template <typename T>
void BaseCSC(T *X, const BaseTCSC &W_csc, T *b, T *Y, int M, int N, int K);
template <typename T>
void CCSC_base(T *X, const CompressedCSC &W_csc, T *b, T *Y, int M, int N, int K);
template <typename T, int UNROLL_FACTOR>
void BaseCSC_unr(T *X, const BaseTCSC &W_csc, T *b, T *Y, int M, int N, int K);
// -------------------------------------------------------------------------------

using comp_func = std::function<void(float *, float *, float *, int, int, int)>;

/* -------------------------------------------------------------------------- */
static void print_usage(const char *prog)
{
    std::cerr << "Usage: " << prog
              << " -M <int> -K <int> -N <int> -s <int>\n";
}

/* -------------------------------------------------------------------------- */
int main(int argc, char **argv)
{
    if (argc < 9)
    { /* need all 5 flags and their arguments             */
        print_usage(argv[0]);
        return 1;
    }

    /* --- Parse command‑line ------------------------------------------------ */
    int M = 0, K = 0, N = 0, nonZero = 0;

    for (int i = 1; i < argc; ++i)
    {
        if (strcmp(argv[i], "-M") == 0 && i + 1 < argc)
            M = std::atoi(argv[++i]);
        else if (strcmp(argv[i], "-K") == 0 && i + 1 < argc)
            K = std::atoi(argv[++i]);
        else if (strcmp(argv[i], "-N") == 0 && i + 1 < argc)
            N = std::atoi(argv[++i]);
        else if (strcmp(argv[i], "-s") == 0 && i + 1 < argc)
            nonZero = std::atoi(argv[++i]);
    }

    if (M <= 0 || K <= 0 || N <= 0 || nonZero <= 0)
    {
        std::cerr << "ERROR: All dimensions must be positive.\n";
        print_usage(argv[0]);
        return 1;
    }

    /* --- Generate matrices ------------------------------------------------- */
    std::vector<int> W_raw = generateSparseMatrix<int>(K, N, nonZero, false, 0);
    std::vector<float> W_FP32(W_raw.begin(), W_raw.end());

    std::vector<float> X = initX<float>(M * K, 512);
    std::vector<float> B(N, 0.0f);
    std::vector<float> Y(M * N, 0.0f);
    std::vector<float> refY(M * N, 0.0f);

    GEMM(X.data(), W_FP32.data(), B.data(), refY.data(), M, N, K);

    /* --- Prepare sparse formats once -------------------------------------- */
    BaseTCSC sf(W_raw.data(), K, N);
    CompressedCSC ccsc(W_raw.data(), K, N);

    /* --- Dispatch to requested kernel ------------------------------------- */
    comp_func kernel;

    kernel = [&ccsc](float *Xv, float *Bv, float *Yv, int m, int n, int k)
    {
        CCSC_base<float>(Xv, ccsc, Bv, Yv, m, n, k);
    };

    /* --- Print matrices ---------------------------------------------------- */
    std::cout << "X matrix  |  W matrix\n";
    int maxRows = (M > K) ? M : K;
    for (int i = 0; i < maxRows; ++i)
    {
        if (i < M)
        {
            for (int j = 0; j < K; ++j)
                std::cout << std::setw(5) << X[i * K + j] << ' ';
        }
        else
        {
            for (int j = 0; j < K; ++j)
                std::cout << std::setw(5) << ' ' << ' ';
        }
        std::cout << "| ";
        if (i < K)
        {
            for (int j = 0; j < N; ++j)
                std::cout << std::setw(5) << W_FP32[i * N + j] << ' ';
        }
        std::cout << '\n';
    }

    std::cout << "B vector:\n";
    for (float b : B)
        std::cout << b << ' ';
    std::cout << "\nExpected Y matrix:\n";
    for (int i = 0; i < M; ++i)
    {
        for (int j = 0; j < N; ++j)
            std::cout << std::setw(5) << refY[i * N + j] << ' ';
        std::cout << '\n';
    }

    /* --- Run kernel -------------------------------------------------------- */
    kernel(X.data(), B.data(), Y.data(), M, N, K);

    std::cout << "Actual Y matrix:\n";
    for (int i = 0; i < M; ++i)
    {
        for (int j = 0; j < N; ++j)
            std::cout << std::setw(5) << Y[i * N + j] << ' ';
        std::cout << '\n';
    }

    /* --- Correctness check ------------------------------------------------- */
    bool ok = compare_results(Y.data(), refY.data(), M, N);
    std::cout << "Correctness: " << (ok ? "PASS" : "FAIL") << std::endl;

    return ok ? 0 : 1;
}