#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <functional>
#include <memory>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <random>  // For std::mt19937, std::random_device, std::uniform_int_distribution
#include <cstdint> // For uint8_t
#include <numeric> // For std::iota (not used for random, but good for sequential)
#include <chrono>  // For seeding random number generator

// #include "perf.h"
// #include "common.h"
#include "sparseUtils.h"

// --- Prototypes for implementations in comp.cpp ---
// These are now declarations of explicitly instantiated templates
template <typename T>
void CSR_base(T *X, const SparseFormatCSC &W_csc, T *b, T *Y, int M, int N, int K);
template <typename T, int UNROLL_FACTOR> // Provide default for UNROLL_FACTOR if used in declaration
void CSR_unrolled(T *X, const SparseFormatCSC &W_csc, T *b, T *Y, int M, int N, int K);
// --- End Prototypes ---

#ifndef M
#define M 32
#endif
#ifndef K
#define K 1024
#endif
#ifndef N
#define N 4096
#endif

#ifndef NONZERO
#define NONZERO 4
#endif

// --- Configuration (Ideally, get these programmatically or for your specific CPU) ---
// Example: 8MB LLC
const size_t TARGET_CACHE_SIZE_BYTES = 8 * 1024 * 1024;
// Common cache line size
const size_t CACHE_LINE_SIZE_BYTES = 64;

void fill_cache_region_with_random_data(std::vector<uint8_t> &buffer)
{
    if (buffer.empty() || CACHE_LINE_SIZE_BYTES == 0)
    {
        std::cerr << "Buffer is empty or cache line size is zero." << std::endl;
        return;
    }

    // Initialize a good quality random number generator
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed);
    std::uniform_int_distribution<uint8_t> distribution(0, 255); // For random bytes

    std::cout << "Attempting to fill a region of " << buffer.size()
              << " bytes with random data..." << std::endl;

    // Iterate through the buffer, writing to each conceptual cache line
    for (size_t i = 0; i < buffer.size(); i += CACHE_LINE_SIZE_BYTES)
    {
        // Write random data to fill one cache line starting at buffer[i]
        for (size_t j = 0; j < CACHE_LINE_SIZE_BYTES; ++j)
        {
            if (i + j < buffer.size())
            { // Ensure we don't write out of bounds
                buffer[i + j] = distribution(generator);
            }
            else
            {
                break; // Reached end of buffer
            }
        }
        // At this point, the cache line corresponding to memory at &buffer[i]
        // has likely been loaded into the cache hierarchy due to the writes.
    }

    // Optional: Perform a "read pass" to ensure data is in cache and to potentially
    // prevent the compiler from optimizing writes away if it thinks the data is unused.
    // A volatile read is stronger, but this often suffices.
    volatile uint8_t dummy_sum = 0; // volatile to hint to compiler not to optimize out
    for (size_t i = 0; i < buffer.size(); i += CACHE_LINE_SIZE_BYTES)
    {
        dummy_sum += buffer[i]; // Read the first byte of each cache line we wrote
    }

    std::cout << "Cache region fill attempt complete. Dummy sum (to ensure usage): "
              << static_cast<int>(dummy_sum) << std::endl;
}

void print_compile_flags()
{
#ifdef COMPULSORY
    printf("Profiling: COMPULSORY\n");
#endif
#ifdef BASE
    printf("Profiling: BASE\n");
#endif
#ifdef GEMM_BASE
    printf("Profiling: GEMM\n");
#endif
#ifdef CSR_BASE
    printf("Profiling: CSR_BASE\n");
#endif
#ifdef CSR_LU2
    printf("Profiling: CSR_LU2\n");
#endif
#ifdef CSR_LU12
    printf("Profiling: CSR_LU12\n");
#endif
}

void compulsory_only()
{
#ifndef COMPULSORY
    return
#else
    // Compute variables for MM computations
    vector<int> W_raw = generateSparseMatrix<int>(K, N, NONZERO, true); // For SparseFormatCSC
    vector<float> X_main = initX<float>(M * K, 512, true);              // Renamed X
    vector<float> B_main(N, 2);                                         // Renamed B
    vector<float> Y_main(M * N, 0);                                     // Renamed refY

#ifdef GEMM_BASE
    vector<float> W_FP32_main(W_raw.begin(), W_raw.end()); // Renamed W_FP32
    int y = W_FP32_main.data()[rand() % W_FP32_main.size()];
#else
    auto sf_csc_data = std::make_shared<SparseFormatCSC>(W_raw.data(), K, N);
    int y = sf_csc_data->col_start_pos[rand() % sf_csc_data->col_start_pos.size()];
#endif

    volatile int x = X_main.data()[rand() % X_main.size()] + B_main.data()[rand() % B_main.size()] + Y_main.data()[rand() % Y_main.size()] + y + W_raw.data()[rand() % W_raw.size()];
    exit(x % 255);
#endif
}

int main(int argc, char **argv)
{
    // Make our experiment deterministic by comparing different methods on the same matrices
    srand(10);

    print_compile_flags();
    compulsory_only();

    // Compute variables for MM computations
    vector<int> W_raw = generateSparseMatrix<int>(K, N, NONZERO, true); // For SparseFormatCSC
    auto sf_csc_data = std::make_shared<SparseFormatCSC>(W_raw.data(), K, N);

    vector<float> X_main = initX<float>(M * K, 512, true); // Renamed X
    vector<float> B_main(N, 2);                            // Renamed B
    vector<float> W_FP32_main(W_raw.begin(), W_raw.end()); // Renamed W_FP32
    vector<float> Y_main(M * N, 0);                        // Renamed refY

    // Calculcate sums to make sure that compiler doesn't optimize out
    // the float vectors on the base (no matrix mult) implementation.
    float sum = 0;
    for (auto i : W_FP32_main)
    {
        sum += i;
    }
    std::cout << "Dummy sum #1 is: " << sum << std::endl;
    CSR_unrolled<float, 2>(X_main.data(), *sf_csc_data, B_main.data(), Y_main.data(), M, N, K);
    sum = 0;
    for (auto i : Y_main)
    {
        sum += i;
    }
    std::cout << "Dummy sum #2 is: " << sum << std::endl;

    std::vector<uint8_t> memory_region(TARGET_CACHE_SIZE_BYTES);
    fill_cache_region_with_random_data(memory_region);
    // At this point, we did our best to ensure a cold cache

#ifndef BASE
#ifdef GEMM_BASE
    GEMM(X_main.data(), W_FP32_main.data(), B_main.data(), Y_main.data(), M, N, K);
#endif
#ifdef CSR_BASE
    CSR_base<float>(X_main.data(), *sf_csc_data, B_main.data(), Y_main.data(), M, N, K);
#endif
#ifdef CSR_LU2
    CSR_unrolled<float, 2>(X_main.data(), *sf_csc_data, B_main.data(), Y_main.data(), M, N, K);
#endif
#ifdef CSR_LU12
    CSR_unrolled<float, 12>(X_main.data(), *sf_csc_data, B_main.data(), Y_main.data(), M, N, K);
#endif
#endif

    // Access elements to again try to make the compiler not optimize out computations
    volatile int x = X_main.data()[rand() % X_main.size()] + B_main.data()[rand() % B_main.size()] + W_FP32_main.data()[rand() % W_FP32_main.size()] + Y_main.data()[rand() % Y_main.size()] + sf_csc_data->col_start_pos[rand() % sf_csc_data->col_start_pos.size()] + W_raw.data()[rand() % W_raw.size()];
    return (x % 255) > 128 ? 0 : 1;
}
