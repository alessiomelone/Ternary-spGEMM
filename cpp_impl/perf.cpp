#include "perf.h"
#include "common.h"
// #error Please comment out the next two lines under linux, then comment this error
// #include "stdafx.h"
// #include <windows.h>

#ifndef WIN32
#include <sys/time.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <vector>

#include "sparseUtils.h"

#ifdef __x86_64__
#include "../include/tsc_x86.h"
#endif

#ifdef __aarch64__
#include "../include/vct_arm.h"
#ifdef PMU
#include "../include/kperf.h"
#endif
#endif

#define NUM_RUNS 1
#define CYCLES_REQUIRED 1e8
#define FREQUENCY 3.2e9

using namespace std;

/*
 * Timing function based on the TimeStep Counter of the CPU.
 */
#ifdef __x86_64__
float rdtsc(comp_func func_to_test, float *X_buf, float *B_buf, float *Y_buf, int M_arg, int N_arg, int K_arg)
{
    int i, num_runs_actual;
    myInt64 cycles_val;
    myInt64 start_val;
    num_runs_actual = NUM_RUNS;

#ifdef CALIBRATE
    while (num_runs_actual < (1 << 14))
    {
        start_val = start_tsc();
        for (i = 0; i < num_runs_actual; ++i)
        {
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
    return (float)cycles_val;
}
#endif

#ifdef __aarch64__
float rdvct(comp_func func_to_test, float *X_buf, float *B_buf, float *Y_buf, int M_arg, int N_arg, int K_arg)
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
    return (float)cycles_val;
}

#ifdef PMU
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
        float cycles_pmu = endperf.cycles - startperf.cycles;
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

float c_clock(comp_func func_to_test, float *X_buf, float *B_buf, float *Y_buf, int M_arg, int N_arg, int K_arg)
{
    int i, num_runs_actual;
    float cycles_val;
    clock_t start_clk, end_clk;

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

        cycles_val = (float)(end_clk - start_clk);
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

    return (float)(end_clk - start_clk) / num_runs_actual;
}

#ifndef WIN32
float timeofday(comp_func func_to_test, float *X_buf, float *B_buf, float *Y_buf, int M_arg, int N_arg, int K_arg)
{
    int i, num_runs_actual;
    float cycles_val;
    struct timeval start_tv, end_tv;

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

        cycles_val = (float)((end_tv.tv_sec - start_tv.tv_sec) + (end_tv.tv_usec - start_tv.tv_usec) / 1e6) * FREQUENCY;
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

    return (float)((end_tv.tv_sec - start_tv.tv_sec) + (end_tv.tv_usec - start_tv.tv_usec) / 1e6) / num_runs_actual;
}

#else // For WIN32

float gettickcount(comp_func func_to_test, float *X_buf, float *B_buf, float *Y_buf, int M_arg, int N_arg, int K_arg)
{
    int i, num_runs_actual;
    float cycles_val, start_tc, end_tc;

    num_runs_actual = NUM_RUNS;
#ifdef CALIBRATE
    while (num_runs_actual < (1 << 14))
    {
        start_tc = (float)GetTickCount();
        for (i = 0; i < num_runs_actual; ++i)
        {
            func_to_test(X_buf, B_buf, Y_buf, M_arg, N_arg, K_arg);
        }
        end_tc = (float)GetTickCount();

        cycles_val = (end_tc - start_tc) * FREQUENCY / 1e3;
        if (cycles_val >= CYCLES_REQUIRED)
            break;

        num_runs_actual *= 2;
    }
#endif

    start_tc = (float)GetTickCount();
    for (i = 0; i < num_runs_actual; ++i)
    {
        func_to_test(X_buf, B_buf, Y_buf, M_arg, N_arg, K_arg);
    }
    end_tc = (float)GetTickCount();

    return (end_tc - start_tc) / num_runs_actual;
}

float queryperfcounter(comp_func func_to_test, float *X_buf, float *B_buf, float *Y_buf, int M_arg, int N_arg, int K_arg, LARGE_INTEGER f)
{
    int i, num_runs_actual;
    float cycles_val;
    LARGE_INTEGER start_pc, end_pc;

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

        cycles_val = (float)(end_pc.QuadPart - start_pc.QuadPart);
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

    return (float)(end_pc.QuadPart - start_pc.QuadPart) / num_runs_actual;
}

#endif // WIN32

float perf_test(comp_func f, int M_param, int K_param, int N_param, int nonZero)
{

    srand((unsigned)time(NULL));

    vector<float> X_perf = initX<float>(M_param * K_param, 512);
    vector<float> Y_perf(M_param * N_param, 0);
    vector<float> B_perf(N_param, 2);

    Y_perf.insert(Y_perf.end(), 10, 0);
    X_perf.insert(X_perf.end(), 10, 0);

#ifdef __x86_64__
    return rdtsc(f, X_perf.data(), B_perf.data(), Y_perf.data(), M_param, N_param, K_param);
#elif defined(__aarch64__) && defined(PMU)
    struct performance_counters p = rdpmu(f, X_perf.data(), B_perf.data(), Y_perf.data(), M_param, N_param, K_param);
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
    return rdvct(f, X_perf.data(), B_perf.data(), Y_perf.data(), M_param, N_param, K_param);
#elif defined(_WIN32) || defined(WIN32)
    LARGE_INTEGER freq;
    QueryPerformanceFrequency(&freq);
    return queryperfcounter(f, X_perf.data(), B_perf.data(), Y_perf.data(), M_param, N_param, K_param, freq) * (FREQUENCY / (float)freq.QuadPart);
#else // Generic Unix
    return timeofday(f, X_perf.data(), B_perf.data(), Y_perf.data(), M_param, N_param, K_param) * FREQUENCY;
#endif
}