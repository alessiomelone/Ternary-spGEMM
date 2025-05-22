// #error Please comment out the next two lines under linux, then comment this error
// #include "stdafx.h"
// #include <windows.h>
#include "perf.h" // Should now include the std::function based comp_func

#ifndef WIN32
#include <sys/time.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdbool.h> // Can use bool from C++ directly
#include <string.h>
#include <vector>      // <-- Add this for creating X, B, Y in perf_test
// #include <functional> // Not strictly needed here if perf.h includes it for comp_func

#include "SparseGEMM.h" // For initX, generateSparseMatrix (though generateSparseMatrix not used directly in perf_test now)

#ifdef __x86_64__
#include "../include/tsc_x86.h" // Ensure path is correct
#endif

#ifdef __aarch64__
#include "../include/vct_arm.h" // Ensure path is correct
#ifdef PMU
#include "../include/kperf.h"   // Ensure path is correct
#endif
// #include <vector> // Already included above
#endif

#define NUM_RUNS 1
#define CYCLES_REQUIRED 1e8 // This is also defined in main.cpp, define once (e.g. in common.h or a config.h)
#define FREQUENCY 3.2e9
#define CALIBRATE

using namespace std; // Generally fine for .cpp files, avoid in headers

/*
 * Timing function based on the TimeStep Counter of the CPU.
 */
#ifdef __x86_64__
// Old: double rdtsc(comp_func func_to_test, float *X, SparseFormat *sparse_W, float *B, float *Y, int M, int N, int K)
// New: comp_func is std::function, doesn't need SparseFormat *sparse_W. X, B, Y are buffers for the test.
double rdtsc(comp_func func_to_test, float *X_buf, float *B_buf, float *Y_buf, int M_arg, int N_arg, int K_arg)
{
    int i, num_runs_actual; // Renamed num_runs
    myInt64 cycles_val;     // Renamed cycles
    myInt64 start_val;      // Renamed start
    num_runs_actual = NUM_RUNS;

#ifdef CALIBRATE
    while (num_runs_actual < (1 << 14))
    {
        start_val = start_tsc();
        for (i = 0; i < num_runs_actual; ++i)
        {
            // Y_buf is modified by func_to_test. For calibration, this is usually fine.
            // If strict identical calls are needed, Y_buf might need resetting or to be a scratchpad.
            func_to_test(X_buf, B_buf, Y_buf, M_arg, N_arg, K_arg);
        }
        cycles_val = stop_tsc(start_val);

        if (cycles_val >= CYCLES_REQUIRED)
            break;

        num_runs_actual *= 2;
    }
#endif

    start_val = start_tsc();
    for (i = 0; i < num_runs_actual; ++i)
    {
        func_to_test(X_buf, B_buf, Y_buf, M_arg, N_arg, K_arg);
    }

    cycles_val = stop_tsc(start_val) / num_runs_actual;
    return (double)cycles_val;
}
#endif

#ifdef __aarch64__
// Old: double rdvct(comp_func func_to_test, float *X, SparseFormat *sparse_W, float *B, float *Y, int M, int N, int K, int nonZero)
// New: nonZero might not be needed by rdvct itself if it was for sparse_W. Let's remove it from rdvct's signature for now.
// If any timer intrinsically needs nonZero (density), it can be added back.
double rdvct(comp_func func_to_test, float *X_buf, float *B_buf, float *Y_buf, int M_arg, int N_arg, int K_arg /*, int nonZero_arg (if needed by timer) */)
{
    int i, num_runs_actual;
    TIMESTAMP cycles_val;
    TIMESTAMP start_val;
    num_runs_actual = NUM_RUNS;

#ifdef CALIBRATE
    while (num_runs_actual < (1 << 14))
    {
        start_val = start_vct();
        for (i = 0; i < num_runs_actual; ++i)
        {
            func_to_test(X_buf, B_buf, Y_buf, M_arg, N_arg, K_arg);
        }
        cycles_val = stop_vct(start_val);

        if (cycles_val >= CYCLES_REQUIRED)
            break;

        num_runs_actual *= 2;
    }
#endif

    start_val = start_vct();
    for (i = 0; i < num_runs_actual; ++i)
    {
        func_to_test(X_buf, B_buf, Y_buf, M_arg, N_arg, K_arg);
    }

    cycles_val = stop_vct(start_val) / num_runs_actual;
    return (double)cycles_val;
}

#ifdef PMU
// Old: struct performance_counters rdpmu(comp_func func_to_test, float *X, SparseFormat *sparse_W, float *B, float *Y, int M, int N, int K)
// New:
struct performance_counters rdpmu(comp_func func_to_test, float *X_buf, float *B_buf, float *Y_buf, int M_arg, int N_arg, int K_arg)
{
    kperf_init();
    int i, num_runs_actual;
    struct performance_counters startperf, endperf, result;
    num_runs_actual = NUM_RUNS;

#ifdef CALIBRATE
    while (num_runs_actual < (1 << 14))
    {
        startperf = kperf_get_counters();
        for (i = 0; i < num_runs_actual; ++i)
        {
            func_to_test(X_buf, B_buf, Y_buf, M_arg, N_arg, K_arg);
        }
        endperf = kperf_get_counters();
        double cycles_pmu = endperf.cycles - startperf.cycles; // Renamed cycles
        if (cycles_pmu >= CYCLES_REQUIRED)
            break;

        num_runs_actual *= 2;
    }
#endif
    startperf = kperf_get_counters();
    for (i = 0; i < num_runs_actual; ++i)
    {
        func_to_test(X_buf, B_buf, Y_buf, M_arg, N_arg, K_arg);
    }

    endperf = kperf_get_counters();
    result.cycles = (endperf.cycles - startperf.cycles) / num_runs_actual;
    result.instructions = (endperf.instructions - startperf.instructions) / num_runs_actual;
    result.branches = (endperf.branches - startperf.branches) / num_runs_actual;
    result.branch_misses = (endperf.branch_misses - startperf.branch_misses) / num_runs_actual;
    result.retired_uops = (endperf.retired_uops - startperf.retired_uops) / num_runs_actual;
    result.int_uops = (endperf.int_uops - startperf.int_uops) / num_runs_actual;
    result.simdfp_uops = (endperf.simdfp_uops - startperf.simdfp_uops) / num_runs_actual;
    result.loadstore_uops = (endperf.loadstore_uops - startperf.loadstore_uops) / num_runs_actual;

    return result;
}
#endif // PMU
#endif // __aarch64__

// Old: double c_clock(comp_func func_to_test, float *X, SparseFormat *sparse_W, float *B, float *Y, int M, int N, int K, int nonZero)
// New: nonZero might be unneeded by c_clock itself.
double c_clock(comp_func func_to_test, float *X_buf, float *B_buf, float *Y_buf, int M_arg, int N_arg, int K_arg /*, int nonZero_arg if needed */)
{
    int i, num_runs_actual;
    double cycles_val;
    clock_t start_clk, end_clk; // Renamed start, end

    num_runs_actual = NUM_RUNS;
#ifdef CALIBRATE
    while (num_runs_actual < (1 << 14))
    {
        start_clk = clock();
        for (i = 0; i < num_runs_actual; ++i)
        {
            func_to_test(X_buf, B_buf, Y_buf, M_arg, N_arg, K_arg);
        }
        end_clk = clock();

        cycles_val = (double)(end_clk - start_clk);
        if (cycles_val >= CYCLES_REQUIRED / (FREQUENCY / CLOCKS_PER_SEC))
            break;

        num_runs_actual *= 2;
    }
#endif

    start_clk = clock();
    for (i = 0; i < num_runs_actual; ++i)
    {
        func_to_test(X_buf, B_buf, Y_buf, M_arg, N_arg, K_arg);
    }
    end_clk = clock();

    return (double)(end_clk - start_clk) / num_runs_actual;
}

#ifndef WIN32
// Old: double timeofday(comp_func func_to_test, float *X, SparseFormat *sparse_W, float *B, float *Y, int M, int N, int K, int nonZero)
// New: nonZero might be unneeded.
double timeofday(comp_func func_to_test, float *X_buf, float *B_buf, float *Y_buf, int M_arg, int N_arg, int K_arg /*, int nonZero_arg if needed */)
{
    int i, num_runs_actual;
    double cycles_val;
    struct timeval start_tv, end_tv; // Renamed start, end

    num_runs_actual = NUM_RUNS;
#ifdef CALIBRATE
    while (num_runs_actual < (1 << 14))
    {
        gettimeofday(&start_tv, NULL);
        for (i = 0; i < num_runs_actual; ++i)
        {
            func_to_test(X_buf, B_buf, Y_buf, M_arg, N_arg, K_arg);
        }
        gettimeofday(&end_tv, NULL);

        cycles_val = (double)((end_tv.tv_sec - start_tv.tv_sec) + (end_tv.tv_usec - start_tv.tv_usec) / 1e6) * FREQUENCY;
        if (cycles_val >= CYCLES_REQUIRED)
            break;

        num_runs_actual *= 2;
    }
#endif

    gettimeofday(&start_tv, NULL);
    for (i = 0; i < num_runs_actual; ++i)
    {
        func_to_test(X_buf, B_buf, Y_buf, M_arg, N_arg, K_arg);
    }
    gettimeofday(&end_tv, NULL);

    return (double)((end_tv.tv_sec - start_tv.tv_sec) + (end_tv.tv_usec - start_tv.tv_usec) / 1e6) / num_runs_actual;
}

#else // For WIN32

// Old: double gettickcount(comp_func func_to_test, float *X, SparseFormat *sparse_W, float *B, float *Y, int M, int N, int K, int nonZero)
// New: nonZero might be unneeded.
double gettickcount(comp_func func_to_test, float *X_buf, float *B_buf, float *Y_buf, int M_arg, int N_arg, int K_arg /*, int nonZero_arg if needed */)
{
    int i, num_runs_actual;
    double cycles_val, start_tc, end_tc; // Renamed start, end

    num_runs_actual = NUM_RUNS;
#ifdef CALIBRATE
    while (num_runs_actual < (1 << 14))
    {
        start_tc = (double)GetTickCount();
        for (i = 0; i < num_runs_actual; ++i)
        {
            func_to_test(X_buf, B_buf, Y_buf, M_arg, N_arg, K_arg);
        }
        end_tc = (double)GetTickCount();

        cycles_val = (end_tc - start_tc) * FREQUENCY / 1e3;
        if (cycles_val >= CYCLES_REQUIRED)
            break;

        num_runs_actual *= 2;
    }
#endif

    start_tc = (double)GetTickCount();
    for (i = 0; i < num_runs_actual; ++i)
    {
        func_to_test(X_buf, B_buf, Y_buf, M_arg, N_arg, K_arg);
    }
    end_tc = (double)GetTickCount();

    return (end_tc - start_tc) / num_runs_actual;
}

// Old: double queryperfcounter(comp_func func_to_test, float *X, SparseFormat *sparse_W, float *B, float *Y, int M, int N, int K, int nonZero, LARGE_INTEGER f)
// New: nonZero might be unneeded by the timer itself.
double queryperfcounter(comp_func func_to_test, float *X_buf, float *B_buf, float *Y_buf, int M_arg, int N_arg, int K_arg, /* int nonZero_arg,*/ LARGE_INTEGER f)
{
    int i, num_runs_actual;
    double cycles_val;
    LARGE_INTEGER start_pc, end_pc; // Renamed start, end

    num_runs_actual = NUM_RUNS;
#ifdef CALIBRATE
    while (num_runs_actual < (1 << 14))
    {
        QueryPerformanceCounter(&start_pc);
        for (i = 0; i < num_runs_actual; ++i)
        {
            func_to_test(X_buf, B_buf, Y_buf, M_arg, N_arg, K_arg);
        }
        QueryPerformanceCounter(&end_pc);

        cycles_val = (double)(end_pc.QuadPart - start_pc.QuadPart);
        if (cycles_val >= CYCLES_REQUIRED / (FREQUENCY / f.QuadPart))
            break;

        num_runs_actual *= 2;
    }
#endif

    QueryPerformanceCounter(&start_pc);
    for (i = 0; i < num_runs_actual; ++i)
    {
        func_to_test(X_buf, B_buf, Y_buf, M_arg, N_arg, K_arg);
    }
    QueryPerformanceCounter(&end_pc);

    return (double)(end_pc.QuadPart - start_pc.QuadPart) / num_runs_actual;
}

#endif // WIN32

// perf_test function signature remains the same, as 'f' (comp_func) now encapsulates the sparse data.
// 'nonZero' is used here if needed by initX or generateSparseMatrix.
// However, since generateSparseMatrix is called in main.cpp before lambdas are created,
// perf_test actually doesn't need to generate the W matrix itself.
// It only needs to generate X, B, Y for the benchmark run.
double perf_test(comp_func f, int M_param, int K_param, int N_param, int nonZero_param [[maybe_unused]]) // nonZero_param might be unused now by perf_test directly
{
    // Mark nonZero_param as potentially unused if SparseGEMM.h's initX or other local utilities don't need it.
    // It was primarily for generateSparseMatrix, which is now done in main.cpp before lambda creation.
    srand((unsigned)time(NULL));

    // Create X, B, Y buffers for the benchmark run.
    // The sparse matrix W is already captured within the std::function 'f'.
    vector<float> X_perf = initX<float>(M_param * K_param, 512);
    vector<float> Y_perf(M_param * N_param, 0); // Output buffer, will be written into by 'f'
    vector<float> B_perf(N_param, 2);

    Y_perf.insert(Y_perf.end(), 10, 0);
    X_perf.insert(X_perf.end(), 10, 0);

#ifdef __x86_64__
    return rdtsc(f, X_perf.data(), B_perf.data(), Y_perf.data(), M_param, N_param, K_param);
#elif defined(__aarch64__) && defined(PMU)
    struct performance_counters p = rdpmu(f, X_perf.data(), B_perf.data(), Y_perf.data(), M_param, N_param, K_param);
    // Original printf for PMU counters
    printf("\n"
           "Instructions     : %.3e\n"
           "Branches         : %.3e\n"
           "Branch Misses    : %.3e\n"
           "Retired Uops     : %.3e\n"
           "Int Uops         : %.3e\n"
           "SIMD FP Uops     : %.3e\n"
           "Load/Store Uops  : %.3e",
           p.instructions,
           p.branches,
           p.branch_misses,
           p.retired_uops,
           p.int_uops,
           p.simdfp_uops,
           p.loadstore_uops);
    return p.cycles;
#elif defined(__aarch64__)
    // Pass nonZero_param to rdvct only if rdvct intrinsically needs it, otherwise remove.
    // Assuming rdvct doesn't need nonZero if it was for sparse_W details.
    return rdvct(f, X_perf.data(), B_perf.data(), Y_perf.data(), M_param, N_param, K_param /*, nonZero_param */);
#elif defined(_WIN32) || defined(WIN32)
    LARGE_INTEGER freq;
    QueryPerformanceFrequency(&freq);
    // Pass nonZero_param to queryperfcounter only if it intrinsically needs it.
    return queryperfcounter(f, X_perf.data(), B_perf.data(), Y_perf.data(), M_param, N_param, K_param, /* nonZero_param,*/ freq) * (FREQUENCY / (double)freq.QuadPart);
#else // Generic Unix
    // Pass nonZero_param to timeofday only if it intrinsically needs it.
    return timeofday(f, X_perf.data(), B_perf.data(), Y_perf.data(), M_param, N_param, K_param /*, nonZero_param */) * FREQUENCY;
#endif
}