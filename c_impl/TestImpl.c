// #error Please comment out the next two lines under linux, then comment this error
// #include "stdafx.h"  //Visual studio expects this line to be the first one, comment out if different compiler
// #include <windows.h> // Include if under windows

#ifndef WIN32
#include <sys/time.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include "sparse_format.h"

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
#define FREQUENCY 3.22e9
#define CALIBRATE

#define RED "\033[0;31m"
#define GREEN "\033[0;32m"
#define RESET "\033[0m"

/*
 * Timing function based on the TimeStep Counter of the CPU.
 */
#ifdef __x86_64__
float rdtsc(float *X, ternarySparseFormat *sparse_W, float *B, float *Y, int M, int N, int K)
{

    int i, num_runs;
    myInt64 cycles;
    myInt64 start;
    num_runs = NUM_RUNS;

    /*
     * The CPUID instruction serializes the pipeline.
     * Using it, we can create execution barriers around the code we want to time.
     * The calibrate section is used to make the computation large enough so as to
     * avoid measurements bias due to the timing overhead.
     */
#ifdef CALIBRATE
    while (num_runs < (1 << 14))
    {
        start = start_tsc();
        for (i = 0; i < num_runs; ++i)
        {
            sparseGEMM(X, sparse_W, B, Y, M, N, K);
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
        sparseGEMM(X, sparse_W, B, Y, M, N, K);
    }

    cycles = stop_tsc(start) / num_runs;
    return (float)cycles;
}
#endif

#ifdef __aarch64__
float rdvct(float *X, ternarySparseFormat *sparse_W, float *B, float *Y, int M, int N, int K, int nonZero)
{
    int i, num_runs;
    TIMESTAMP cycles;
    TIMESTAMP start;
    num_runs = NUM_RUNS;

    /*
     * The CPUID instruction serializes the pipeline.
     * Using it, we can create execution barriers around the code we want to time.
     * The calibrate section is used to make the computation large enough so as to
     * avoid measurements bias due to the timing overhead.
     */
#ifdef CALIBRATE
    while (num_runs < (1 << 14))
    {
        start = start_vct();
        for (i = 0; i < num_runs; ++i)
        {
            sparseGEMM(X, sparse_W, B, Y, M, N, K);
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
        sparseGEMM(X, sparse_W, B, Y, M, N, K);
    }

    cycles = stop_vct(start) / num_runs;
    return (float)cycles;
}

#ifdef PMU
struct performance_counters rdpmu(float *X, ternarySparseFormat *sparse_W, float *B, float *Y, int M, int N, int K)
{
    kperf_init();
    int i, num_runs;
    struct performance_counters startperf, endperf, result;
    num_runs = NUM_RUNS;

    /*
     * The CPUID instruction serializes the pipeline.
     * Using it, we can create execution barriers around the code we want to time.
     * The calibrate section is used to make the computation large enough so as to
     * avoid measurements bias due to the timing overhead.
     */
#ifdef CALIBRATE
    while (num_runs < (1 << 14))
    {
        startperf = kperf_get_counters();
        for (i = 0; i < num_runs; ++i)
        {
            sparseGEMM(X, sparse_W, B, Y, M, N, K);
        }
        endperf = kperf_get_counters();
        float cycles = endperf.cycles - startperf.cycles;
        if (cycles >= CYCLES_REQUIRED)
            break;

        num_runs *= 2;
    }
#endif

    startperf = kperf_get_counters();
    for (i = 0; i < num_runs; ++i)
    {
        sparseGEMM(X, sparse_W, B, Y, M, N, K);
    }

    endperf = kperf_get_counters();
    result.cycles = (endperf.cycles - startperf.cycles) / num_runs;
    result.instructions = (endperf.instructions - startperf.instructions) / num_runs;
    result.branches = (endperf.branches - startperf.branches) / num_runs;
    result.branch_misses = (endperf.branch_misses - startperf.branch_misses) / num_runs;

    return result;
}
#endif

#endif

float c_clock(float *X, ternarySparseFormat *sparse_W, float *B, float *Y, int M, int N, int K, int nonZero)
{
    int i, num_runs;
    float cycles;
    clock_t start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while (num_runs < (1 << 14))
    {
        start = clock();
        for (i = 0; i < num_runs; ++i)
        {
            sparseGEMM(X, sparse_W, B, Y, M, N, K);
        }
        end = clock();

        cycles = (float)(end - start);

        // Same as in c_clock: CYCLES_REQUIRED should be expressed accordingly to the order of magnitude of CLOCKS_PER_SEC
        if (cycles >= CYCLES_REQUIRED / (FREQUENCY / CLOCKS_PER_SEC))
            break;

        num_runs *= 2;
    }
#endif

    start = clock();
    for (i = 0; i < num_runs; ++i)
    {
        sparseGEMM(X, sparse_W, B, Y, M, N, K);
    }
    end = clock();

    return (float)(end - start) / num_runs;
}

#ifndef WIN32
float timeofday(float *X, ternarySparseFormat *sparse_W, float *B, float *Y, int M, int N, int K, int nonZero)
{
    int i, num_runs;
    float cycles;
    struct timeval start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while (num_runs < (1 << 14))
    {
        gettimeofday(&start, NULL);
        for (i = 0; i < num_runs; ++i)
        {
            sparseGEMM(X, sparse_W, B, Y, M, N, K);
        }
        gettimeofday(&end, NULL);

        cycles = (float)((end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1e6) * FREQUENCY;

        if (cycles >= CYCLES_REQUIRED)
            break;

        num_runs *= 2;
    }
#endif

    gettimeofday(&start, NULL);
    for (i = 0; i < num_runs; ++i)
    {
        sparseGEMM(X, sparse_W, B, Y, M, N, K);
    }
    gettimeofday(&end, NULL);

    return (float)((end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1e6) / num_runs;
}

#else

float gettickcount(float *X, ternarySparseFormat *sparse_W, float *B, float *Y, int M, int N, int K, int nonZero)
{
    int i, num_runs;
    float cycles, start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while (num_runs < (1 << 14))
    {
        start = (float)GetTickCount();
        for (i = 0; i < num_runs; ++i)
        {
            sparseGEMM(X, sparse_W, B, Y, M, N, K);
        }
        end = (float)GetTickCount();

        cycles = (end - start) * FREQUENCY / 1e3; // end-start provides a measurement in the order of milliseconds

        if (cycles >= CYCLES_REQUIRED)
            break;

        num_runs *= 2;
    }
#endif

    start = (float)GetTickCount();
    for (i = 0; i < num_runs; ++i)
    {
        sparseGEMM(X, sparse_W, B, Y, M, N, K);
    }
    end = (float)GetTickCount();

    return (end - start) / num_runs;
}

float queryperfcounter(float *X, ternarySparseFormat *sparse_W, float *B, float *Y, int M, int N, int K, int nonZero, LARGE_INTEGER f)
{
    int i, num_runs;
    float cycles;
    LARGE_INTEGER start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while (num_runs < (1 << 14))
    {
        QueryPerformanceCounter(&start);
        for (i = 0; i < num_runs; ++i)
        {
            sparseGEMM(X, sparse_W, B, Y, M, N, K);
        }
        QueryPerformanceCounter(&end);

        cycles = (float)(end.QuadPart - start.QuadPart);

        // Same as in c_clock: CYCLES_REQUIRED should be expressed accordingly to the order of magnitude of f
        if (cycles >= CYCLES_REQUIRED / (FREQUENCY / f.QuadPart))
            break;

        num_runs *= 2;
    }
#endif

    QueryPerformanceCounter(&start);
    for (i = 0; i < num_runs; ++i)
    {
        sparseGEMM(X, sparse_W, B, Y, M, N, K);
    }
    QueryPerformanceCounter(&end);

    return (float)(end.QuadPart - start.QuadPart) / num_runs;
}

#endif

int main(int argc, char **argv)
{
    srand((unsigned)time(NULL));

    int M = 0, K = 0, N = 0, nonZero = 0;
    
    // Check all arguments passed
    if (argc < 9) {
        fprintf(stderr, "Usage: %s -M <int> -K <int> -N <int> -s <int>\n", argv[0]);
        return 1;
    }
    
    M = atoi(argv[2]);
    K = atoi(argv[4]);
    N = atoi(argv[6]);
    nonZero = atoi(argv[8]);

    // Check for valid parameters
    if (M <= 0 || K <= 0 || N <= 0 || nonZero <= 0)
    {
        fprintf(stderr, "ERROR: All dimensions must be positive integers.\n");
        fprintf(stderr, "Usage: %s -M <int> -K <int> -N <int> -s <int>\n", argv[0]);
        return 1;
    }

    float *X = initX(M * K, 512);
    int *W = generateSparseMatrix(K, N, nonZero, true);

    // Initialize output and reference matrices
    float *Y = (float *)calloc(M * N, sizeof(float));
    float *B = (float *)calloc(N, sizeof(float)); // Initialize bias to zeros
    float *refY = (float *)calloc(M * N, sizeof(float));

    // Run sparse and dense implementations
    ternarySparseFormat *sparse_W = convertTernaryToSparseFormat(W, K, N, nonZero);
    sparseGEMM(X, sparse_W, B, Y, M, N, K);

#ifdef VALIDATE
    GEMM(X, W, B, refY, M, N, K);

    // Compare results
    if (compare_results(Y, refY, M, N))
    {
        printf("%sTest case passed!%s\n", GREEN, RESET);
    }
    else
    {
        printf("%sTest case failed!%s\n", RED, RESET);
        destroyTernarySparseFormat(sparse_W);
        free(refY);
        free(B);
        free(Y);
        free(X);
        free(W);
        return 0;
    }
#endif

#ifdef __x86_64__
    float r = rdtsc(X, sparse_W, B, Y, M, N, K);
    printf("# RDTSC instruction\n");
    printf("rdtsc_cycles=%.0lf\n", r);
    printf("rdtsc_seconds=%.8lf\n", r / (FREQUENCY));
    printf("rdtsc_freq_mhz=%.2lf\n", (FREQUENCY) / 1e6);
#endif

    float c = c_clock(X, sparse_W, B, Y, M, N, K, nonZero);
    printf("\n# C clock() function\n");
    printf("clock_cycles=%.0lf\n", c);
    printf("clock_freq_mhz=%.2lf\n", (float)CLOCKS_PER_SEC / 1e6);
    printf("clock_seconds=%.8lf\n", c / CLOCKS_PER_SEC);

#ifndef WIN32
    float t = timeofday(X, sparse_W, B, Y, M, N, K, nonZero);
    printf("\n# C gettimeofday() function\n");
    printf("timeofday_seconds=%.8lf\n", t);
#else
    LARGE_INTEGER f;
    float t = gettickcount(X, sparse_W, B, Y, M, N, K, nonZero);
    printf("# Windows getTickCount() function\n");
    printf("gettickcount_milliseconds=%.3lf\n", t);
    QueryPerformanceFrequency(&f);
    float p = queryperfcounter(X, sparse_W, B, Y, M, N, K, nonZero, f);
    printf("\n# Windows QueryPerformanceCounter() function\n");
    printf("queryperfcounter_cycles=%.0lf\n", p);
    printf("queryperfcounter_seconds=%.8lf\n", p / (float)f.QuadPart);
    printf("queryperfcounter_freq_mhz=%.2lf\n", (float)f.QuadPart / 1000);
#endif

#ifdef __aarch64__
    float v = rdvct(X, sparse_W, B, Y, M, N, K, nonZero);
    printf("\n# VCT instruction\n");
    printf("vct_cycles=%.0lf\n", v);
    printf("vct_seconds=%.8lf\n", v / (get_vct_freq()));
    printf("vct_freq_mhz=%.2lf\n", (get_vct_freq()) / 1e6);

#ifdef PMU
    struct performance_counters p = rdpmu(X, sparse_W, B, Y, M, N, K);
    printf("\n# PMU instruction\n");
    printf("pmu_cycles=%.0lf\n", p.cycles);
    printf("pmu_seconds=%.8lf\n", p.cycles / (FREQUENCY));
    printf("pmu_freq_mhz=%.2lf\n", (FREQUENCY) / 1e6);
    printf("pmu_instructions=%.0lf\n", p.instructions);
    printf("pmu_branches=%.0lf\n", p.branches);
    printf("pmu_branch_misses=%.0lf\n", p.branch_misses);
    printf("pmu_ipc=%.2lf\n", p.instructions / p.cycles);
#endif

#endif
    return 0;
}
