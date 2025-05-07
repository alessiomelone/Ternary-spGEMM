// #error Please comment out the next two lines under linux, then comment this error
// #include "stdafx.h"  //Visual studio expects this line to be the first one, comment out if different compiler
// #include <windows.h> // Include if under windows
#include "perf.h"

#ifndef WIN32
#include <sys/time.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include "SparseGEMM.h"

#ifdef __x86_64__
#include "../include/tsc_x86.h"
#endif

#ifdef __aarch64__
#include "../include/vct_arm.h"
#ifdef PMU
#include "../include/kperf.h"
#endif
#include <vector>
#endif

#define NUM_RUNS 1
#define CYCLES_REQUIRED 1e8
#define FREQUENCY 3.2e9
#define CALIBRATE

using namespace std;

/*
 * Timing function based on the TimeStep Counter of the CPU.
 */
#ifdef __x86_64__
double rdtsc(comp_func func_to_test, float *X, SparseFormat *sparse_W, float *B, float *Y, int M, int N, int K)
{
    int i, num_runs;
    myInt64 cycles;
    myInt64 start;
    num_runs = NUM_RUNS;

#ifdef CALIBRATE
    while (num_runs < (1 << 14))
    {
        start = start_tsc();
        for (i = 0; i < num_runs; ++i)
        {
            func_to_test(X,
                         sparse_W->col_start_pos.data(), sparse_W->col_start_neg.data(),
                         sparse_W->row_index_pos.data(), sparse_W->row_index_neg.data(),
                         B, Y, M, N, K);
        }
        cycles = stop_tsc(start);

        if (cycles >= CYCLES_REQUIRED)
            break;

        num_runs *= 2;
    }
#endif

    start = start_tsc();
    for (i = 0; i < num_runs; ++i)
    {
        func_to_test(X,
                     sparse_W->col_start_pos.data(), sparse_W->col_start_neg.data(),
                     sparse_W->row_index_neg.data(), sparse_W->row_index_pos.data(),
                     B, Y, M, N, K);
    }

    cycles = stop_tsc(start) / num_runs;
    return (double)cycles;
}
#endif

#ifdef __aarch64__
double rdvct(comp_func func_to_test, float *X, SparseFormat *sparse_W, float *B, float *Y, int M, int N, int K, int nonZero)
{
    int i, num_runs;
    TIMESTAMP cycles;
    TIMESTAMP start;
    num_runs = NUM_RUNS;

#ifdef CALIBRATE
    while (num_runs < (1 << 14))
    {
        start = start_vct();
        for (i = 0; i < num_runs; ++i)
        {
            func_to_test(X,
                         sparse_W->col_start_pos.data(), sparse_W->col_start_neg.data(),
                         sparse_W->row_index_pos.data(), sparse_W->row_index_neg.data(),
                         B, Y, M, N, K);
        }
        cycles = stop_vct(start);

        if (cycles >= CYCLES_REQUIRED)
            break;

        num_runs *= 2;
    }
#endif

    start = start_vct();
    for (i = 0; i < num_runs; ++i)
    {
        func_to_test(X,
                     sparse_W->col_start_pos.data(), sparse_W->col_start_neg.data(),
                     sparse_W->row_index_neg.data(), sparse_W->row_index_pos.data(),
                     B, Y, M, N, K);
    }

    cycles = stop_vct(start) / num_runs;
    return (double)cycles;
}

#ifdef PMU
struct performance_counters rdpmu(comp_func func_to_test, float *X, SparseFormat *sparse_W, float *B, float *Y, int M, int N, int K)
{
    kperf_init();
    int i, num_runs;
    struct performance_counters startperf, endperf, result;
    num_runs = NUM_RUNS;

#ifdef CALIBRATE
    while (num_runs < (1 << 14))
    {
        startperf = kperf_get_counters();
        for (i = 0; i < num_runs; ++i)
        {
            func_to_test(X,
                         sparse_W->col_start_pos.data(), sparse_W->col_start_neg.data(),
                         sparse_W->row_index_pos.data(), sparse_W->row_index_neg.data(),
                         B, Y, M, N, K);
        }
        endperf = kperf_get_counters();
        double cycles = endperf.cycles - startperf.cycles;
        if (cycles >= CYCLES_REQUIRED)
            break;

        num_runs *= 2;
    }
#endif
    startperf = kperf_get_counters();
    for (i = 0; i < num_runs; ++i)
    {
        func_to_test(X,
                     sparse_W->col_start_pos.data(), sparse_W->col_start_neg.data(),
                     sparse_W->row_index_neg.data(), sparse_W->row_index_pos.data(),
                     B, Y, M, N, K);
    }

    endperf = kperf_get_counters();
    result.cycles = (endperf.cycles - startperf.cycles) / num_runs;
    result.instructions = (endperf.instructions - startperf.instructions) / num_runs;
    result.branches = (endperf.branches - startperf.branches) / num_runs;
    result.branch_misses = (endperf.branch_misses - startperf.branch_misses) / num_runs;
    result.retired_uops = (endperf.retired_uops - startperf.retired_uops) / num_runs;
    result.int_uops = (endperf.int_uops - startperf.int_uops) / num_runs;
    result.simdfp_uops = (endperf.simdfp_uops - startperf.simdfp_uops) / num_runs;
    result.loadstore_uops = (endperf.loadstore_uops - startperf.loadstore_uops) / num_runs;


    return result;
}
#endif

#endif

double c_clock(comp_func func_to_test, float *X, SparseFormat *sparse_W, float *B, float *Y, int M, int N, int K, int nonZero)
{
    int i, num_runs;
    double cycles;
    clock_t start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while (num_runs < (1 << 14))
    {
        start = clock();
        for (i = 0; i < num_runs; ++i)
        {
            func_to_test(X,
                         sparse_W->col_start_pos.data(), sparse_W->col_start_neg.data(),
                         sparse_W->row_index_pos.data(), sparse_W->row_index_neg.data(),
                         B, Y, M, N, K);
        }
        end = clock();

        cycles = (double)(end - start);
        if (cycles >= CYCLES_REQUIRED / (FREQUENCY / CLOCKS_PER_SEC))
            break;

        num_runs *= 2;
    }
#endif

    start = clock();
    for (i = 0; i < num_runs; ++i)
    {
        func_to_test(X,
                     sparse_W->col_start_pos.data(), sparse_W->col_start_neg.data(),
                     sparse_W->row_index_neg.data(), sparse_W->row_index_pos.data(),
                     B, Y, M, N, K);
    }
    end = clock();

    return (double)(end - start) / num_runs;
}

#ifndef WIN32
double timeofday(comp_func func_to_test, float *X, SparseFormat *sparse_W, float *B, float *Y, int M, int N, int K, int nonZero)
{
    int i, num_runs;
    double cycles;
    struct timeval start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while (num_runs < (1 << 14))
    {
        gettimeofday(&start, NULL);
        for (i = 0; i < num_runs; ++i)
        {
            func_to_test(X,
                         sparse_W->col_start_pos.data(), sparse_W->col_start_neg.data(),
                         sparse_W->row_index_pos.data(), sparse_W->row_index_neg.data(),
                         B, Y, M, N, K);
        }
        gettimeofday(&end, NULL);

        cycles = (double)((end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1e6) * FREQUENCY;
        if (cycles >= CYCLES_REQUIRED)
            break;

        num_runs *= 2;
    }
#endif

    gettimeofday(&start, NULL);
    for (i = 0; i < num_runs; ++i)
    {
        func_to_test(X,
                     sparse_W->col_start_pos.data(), sparse_W->col_start_neg.data(),
                     sparse_W->row_index_neg.data(), sparse_W->row_index_pos.data(),
                     B, Y, M, N, K);
    }
    gettimeofday(&end, NULL);

    return (double)((end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1e6) / num_runs;
}

#else

double gettickcount(comp_func func_to_test, float *X, SparseFormat *sparse_W, float *B, float *Y, int M, int N, int K, int nonZero)
{
    int i, num_runs;
    double cycles, start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while (num_runs < (1 << 14))
    {
        start = (double)GetTickCount();
        for (i = 0; i < num_runs; ++i)
        {
            func_to_test(X,
                         sparse_W->col_start_pos.data(), sparse_W->col_start_neg.data(),
                         sparse_W->row_index_pos.data(), sparse_W->row_index_neg.data(),
                         B, Y, M, N, K);
        }
        end = (double)GetTickCount();

        cycles = (end - start) * FREQUENCY / 1e3; // end-start provides a measurement in the order of milliseconds

        if (cycles >= CYCLES_REQUIRED)
            break;

        num_runs *= 2;
    }
#endif

    start = (double)GetTickCount();
    for (i = 0; i < num_runs; ++i)
    {
        func_to_test(X,
                     sparse_W->col_start_pos.data(), sparse_W->col_start_neg.data(),
                     sparse_W->row_index_neg.data(), sparse_W->row_index_pos.data(),
                     B, Y, M, N, K);
    }
    end = (double)GetTickCount();

    return (end - start) / num_runs;
}

double queryperfcounter(comp_func func_to_test, float *X, SparseFormat *sparse_W, float *B, float *Y, int M, int N, int K, int nonZero, LARGE_INTEGER f)
{
    int i, num_runs;
    double cycles;
    LARGE_INTEGER start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while (num_runs < (1 << 14))
    {
        QueryPerformanceCounter(&start);
        for (i = 0; i < num_runs; ++i)
        {
            func_to_test(X,
                         sparse_W->col_start_pos.data(), sparse_W->col_start_neg.data(),
                         sparse_W->row_index_pos.data(), sparse_W->row_index_neg.data(),
                         B, Y, M, N, K);
        }
        QueryPerformanceCounter(&end);

        cycles = (double)(end.QuadPart - start.QuadPart);

        if (cycles >= CYCLES_REQUIRED / (FREQUENCY / f.QuadPart))
            break;

        num_runs *= 2;
    }
#endif

    QueryPerformanceCounter(&start);
    for (i = 0; i < num_runs; ++i)
    {
        func_to_test(X,
                     sparse_W->col_start_pos.data(), sparse_W->col_start_neg.data(),
                     sparse_W->row_index_neg.data(), sparse_W->row_index_pos.data(),
                     B, Y, M, N, K);
    }
    QueryPerformanceCounter(&end);

    return (double)(end.QuadPart - start.QuadPart) / num_runs;
}

#endif

double perf_test(comp_func f, int M, int K, int N, int nonZero)
{
    srand((unsigned)time(NULL));

    vector<float> X = initX<float>(M * K, 512);
    vector<int> W = generateSparseMatrix<int>(K, N, nonZero, false);
    vector<float> Y(M * N, 0);
    vector<float> B(N, 2);

    SparseFormat sf = SparseFormat(W.data(), K, N);

#ifdef __x86_64__
    // Su Intel/AMD x86_64: già restituisce i cicli
    return rdtsc(f, X.data(), &sf, B.data(), Y.data(), M, N, K);
#elif defined(__aarch64__) && defined(PMU)
    // Su ARM64 con PMU: usa i cicli hardware
    struct performance_counters p = rdpmu(f, X.data(), &sf, B.data(), Y.data(), M, N, K);
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
    // Su ARM64 senza PMU: assume che rdvct restituisca cicli (se fosse in Hz, va adattato)
    return rdvct(f, X.data(), &sf, B.data(), Y.data(), M, N, K, nonZero);
#elif defined(_WIN32) || defined(WIN32)
    // Su Windows: QueryPerformanceCounter -> cicli stimati tramite frequenza
    LARGE_INTEGER freq;
    QueryPerformanceFrequency(&freq);
    double ticks = queryperfcounter(f, X.data(), &sf, B.data(), Y.data(), M, N, K, nonZero, freq);
    // Conversione: ticks * (FREQUENCY / freq.QuadPart) = cicli
    return ticks * (FREQUENCY / (double)freq.QuadPart);
#else
    // Su Unix generico: timeofday -> secondi, converto in cicli
    double seconds = timeofday(f, X.data(), &sf, B.data(), Y.data(), M, N, K, nonZero);
    // Conversione: secondi * FREQUENCY = cicli
    return seconds * FREQUENCY;
#endif
}
