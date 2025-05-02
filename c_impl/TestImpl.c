//#error Please comment out the next two lines under linux, then comment this error
//#include "stdafx.h"  //Visual studio expects this line to be the first one, comment out if different compiler
//#include <windows.h> // Include if under windows

#ifndef WIN32
#include <sys/time.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
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

// TODO: Run valgrind on it to verify no unsafe mem accesses -- should be ok since tests pass
bool compare_results(int* result, int* groundTruth, int H, int W) {
	for (int h = 0; h < H; h++) {
		for (int w = 0; w < W; w++) {
			int i = h * W + w;
			if (abs(result[i] - groundTruth[i]) > 10e-6) {
				printf("Error at: H=%d, W=%d, result=%d, groundTruth=%d\n", h, w, result[i], groundTruth[i]);
				return false;
			}
		}
	}

	return true;
}

// Do basic GEMM (3-fold loop)
void GEMM(int* X, int* W, int* b, int* Y, int M, int N, int K) {
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < N; n++) {
			int y = 0;
			for (int k = 0; k < K; k++) {
				y += X[m * K + k] * W[k * N + n];
			}
			Y[m * N + n] = y + b[n];
		}
	}
}

int *generateSparseMatrix(int H, int W, int nonZero, bool uniformDistribution) {
    // TODO : Free y
    int *y = malloc(sizeof(int) * H * W);
    if (uniformDistribution) {
        for (int h = 0; h < H; h++) {
            for (int w = 0; w < W; w += nonZero * 2) {
                // Assign +1, -1 to each 2 x nonZero slots
                int randomA = rand() % nonZero * 2;
                int randomB = rand() % nonZero * 2;
                y[w + randomA] = 1;
                while (randomA==randomB) {
                    randomB = rand() % nonZero * 2;
                }
                y[w + randomB] = 1;
            }
        }
    }
    else {
        for (int h = 0; h < H; h++) {
            // Assign +1 to W / nonZero / 2 places
            int count = 0;
            int limit = (W / nonZero) / 2 - 1;
            while (count < limit) {
                int randomA = rand() % W;
                if (y[randomA] == 0) {
                    y[randomA] = 1;
                    count++;
                }
            }

            // Assign -1 to W / nonZero / 2 places
            count = 0;
            while (count < (W / nonZero / 2 - 1)) {
                int randomA = rand() % W;
                if (y[randomA] == 0) {
                    y[randomA] = -1;
                    count++;
                }
            }
        }
    }

    return y;
}

/*
 * Timing function based on the TimeStep Counter of the CPU.
 */
#ifdef __x86_64__
double rdtsc(int* X, ternarySparseFormat* sparse_W, int* B, int* Y, int M, int N, int K) {

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
    while(num_runs < (1 << 14)) {
        start = start_tsc();
        for (i = 0; i < num_runs; ++i) {
            sparseGEMM(X, sparse_W, B, Y, M, N, K);

        }
        cycles = stop_tsc(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    start = start_tsc();
    for (i = 0; i < num_runs; ++i) {
        sparseGEMM(X, sparse_W, B, Y, M, N, K);

    }

    cycles = stop_tsc(start)/num_runs;
    return (double) cycles;
}
#endif

#ifdef __aarch64__
double rdvct(int* X, ternarySparseFormat* sparse_W, int* B, int* Y, int M, int N, int K, int nonZero) {
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
    while(num_runs < (1 << 14)) {
        start = start_vct();
        for (i = 0; i < num_runs; ++i) {
            sparseGEMM(X, sparse_W, B, Y, M, N, K);

        }
        cycles = stop_vct(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    start = start_vct();
    for (i = 0; i < num_runs; ++i) {
        sparseGEMM(X, sparse_W, B, Y, M, N, K);

    }

    cycles = stop_vct(start)/num_runs;
    return (double) cycles;
}

#ifdef PMU
struct performance_counters rdpmu(int* X, ternarySparseFormat* sparse_W, int* B, int* Y, int M, int N, int K) {
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
    while(num_runs < (1 << 14)) {
        startperf = kperf_get_counters();
        for (i = 0; i < num_runs; ++i) {
            sparseGEMM(X, sparse_W, B, Y, M, N, K);

        }
        endperf = kperf_get_counters();
        double cycles = endperf.cycles - startperf.cycles;
        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    startperf = kperf_get_counters();
    for (i = 0; i < num_runs; ++i) {
        sparseGEMM(X, sparse_W, B, Y, M, N, K);

    }

    endperf = kperf_get_counters();
    result.cycles = (endperf.cycles - startperf.cycles)/num_runs;
    result.instructions = (endperf.instructions - startperf.instructions)/num_runs;
    result.branches = (endperf.branches - startperf.branches)/num_runs;
    result.branch_misses = (endperf.branch_misses - startperf.branch_misses)/num_runs;

    return result;
}
#endif

#endif

double c_clock(int* X, ternarySparseFormat* sparse_W, int* B, int* Y, int M, int N, int K, int nonZero) {
    int i, num_runs;
    double cycles;
    clock_t start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while(num_runs < (1 << 14)) {
        start = clock();
        for (i = 0; i < num_runs; ++i) {
            sparseGEMM(X, sparse_W, B, Y, M, N, K);

        }
        end = clock();

        cycles = (double)(end-start);

        // Same as in c_clock: CYCLES_REQUIRED should be expressed accordingly to the order of magnitude of CLOCKS_PER_SEC
        if(cycles >= CYCLES_REQUIRED/(FREQUENCY/CLOCKS_PER_SEC)) break;

        num_runs *= 2;
    }
#endif

    start = clock();
    for(i=0; i<num_runs; ++i) {
        sparseGEMM(X, sparse_W, B, Y, M, N, K);

    }
    end = clock();

    return (double)(end-start)/num_runs;
}

#ifndef WIN32
double timeofday(int* X, ternarySparseFormat* sparse_W, int* B, int* Y, int M, int N, int K, int nonZero) {
    int i, num_runs;
    double cycles;
    struct timeval start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while(num_runs < (1 << 14)) {
        gettimeofday(&start, NULL);
        for (i = 0; i < num_runs; ++i) {
            sparseGEMM(X, sparse_W, B, Y, M, N, K);

        }
        gettimeofday(&end, NULL);

        cycles = (double)((end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1e6)*FREQUENCY;

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    gettimeofday(&start, NULL);
    for(i=0; i < num_runs; ++i) {
        sparseGEMM(X, sparse_W, B, Y, M, N, K);

    }
    gettimeofday(&end, NULL);

    return (double)((end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1e6)/ num_runs;
}

#else

double gettickcount(int* X, ternarySparseFormat* sparse_W, int* B, int* Y, int M, int N, int K, int nonZero) {
    int i, num_runs;
    double cycles, start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while(num_runs < (1 << 14)) {
        start = (double)GetTickCount();
        for (i = 0; i < num_runs; ++i) {
            sparseGEMM(X, sparse_W, B, Y, M, N, K);

        }
        end = (double)GetTickCount();

        cycles = (end-start)*FREQUENCY/1e3; // end-start provides a measurement in the order of milliseconds

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    start = (double)GetTickCount();
    for(i=0; i < num_runs; ++i) {
        sparseGEMM(X, sparse_W, B, Y, M, N, K);

    }
    end = (double)GetTickCount();

    return (end-start)/num_runs;
}

double queryperfcounter(int* X, ternarySparseFormat* sparse_W, int* B, int* Y, int M, int N, int K, int nonZero, LARGE_INTEGER f) {
    int i, num_runs;
    double cycles;
    LARGE_INTEGER start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while(num_runs < (1 << 14)) {
        QueryPerformanceCounter(&start);
        for (i = 0; i < num_runs; ++i) {
            sparseGEMM(X, sparse_W, B, Y, M, N, K);

        }
        QueryPerformanceCounter(&end);

        cycles = (double)(end.QuadPart - start.QuadPart);

        // Same as in c_clock: CYCLES_REQUIRED should be expressed accordingly to the order of magnitude of f
        if(cycles >= CYCLES_REQUIRED/(FREQUENCY/f.QuadPart)) break;

        num_runs *= 2;
    }
#endif

    QueryPerformanceCounter(&start);
    for(i=0; i < num_runs; ++i) {
        sparseGEMM(X, sparse_W, B, Y, M, N, K);

    }
    QueryPerformanceCounter(&end);

    return (double)(end.QuadPart - start.QuadPart)/num_runs;
}

#endif



int main(int argc, char **argv) {
    int M = 0, K = 0, N = 0, nonZero = 0;
    int opt;

    // The colons after each letter mean the option requires an argument
    while ((opt = getopt(argc, argv, "M:K:N:s:")) != -1) {
        switch (opt) {
            case 'M':
                M = atoi(optarg);
                break;
            case 'K':
                K = atoi(optarg);
                break;
            case 'N':
                N = atoi(optarg);
                break;
            case 's':
                nonZero = atoi(optarg);
                break;
            default:
                fprintf(stderr, "Usage: %s -M <int> -K <int> -N <int> -z <int>\n", argv[0]);
                return 1;
        }
    }

    // Check for valid parameters
    if (M <= 0 || K <= 0 || N <= 0 || nonZero <= 0) { 
        fprintf(stderr, "ERROR: All dimensions must be positive integers.\n");
        fprintf(stderr, "Usage: %s -M <int> -K <int> -N <int> -s <int>\n", argv[0]);
        return 1;
    }
    
    int* X = generateSparseMatrix(M, K, 4, true);
    int* W = generateSparseMatrix(K, N, nonZero, true);
    
    // Initialize output and reference matrices
    int *Y = (int *)calloc(M * N, sizeof(int));
    int *B = (int *)calloc(N, sizeof(int));  // Initialize bias to zeros
    int *refY = (int *)calloc(M * N, sizeof(int));
    
    // Run sparse and dense implementations
    ternarySparseFormat *sparse_W = convertTernaryToSparseFormat(W, K, N, nonZero);
    sparseGEMM(X, sparse_W, B, Y, M, N, K);
    GEMM(X, W, B, refY, M, N, K);
    
    // Compare results
    if (compare_results(Y, refY, M, N)) {
        printf("%sTest case passed!%s\n", GREEN, RESET);
    } else {
        printf("%sTest case failed!%s\n", RED, RESET);
        destroyTernarySparceFormat(sparse_W);
        free(refY);
        free(B);
        free(Y);
        free(X);
        free(W);
        return 0;
    }
    
    #ifdef __x86_64__
    double r = rdtsc(X, sparse_W, B, Y, M, N, K);
    printf("# RDTSC instruction\n");
    printf("rdtsc_cycles=%.0lf\n", r);
    printf("rdtsc_seconds=%.8lf\n", r/(FREQUENCY));
    printf("rdtsc_freq_mhz=%.2lf\n", (FREQUENCY)/1e6);
#endif

    double c = c_clock(X, sparse_W, B, Y, M, N, K, nonZero);
    printf("# C clock() function\n");
    printf("clock_cycles=%.0lf\n", c);
    printf("clock_freq_mhz=%.2lf\n", (double) CLOCKS_PER_SEC/1e6);
    printf("clock_seconds=%.8lf\n", c/CLOCKS_PER_SEC);

#ifndef WIN32
    double t = timeofday(X, sparse_W, B, Y, M, N, K, nonZero);
    printf("# C gettimeofday() function\n");
    printf("timeofday_seconds=%.8lf\n", t);
#else
    LARGE_INTEGER f;
    double t = gettickcount(X, sparse_W, B, Y, M, N, K, nonZero);
    printf("# Windows getTickCount() function\n");
    printf("gettickcount_milliseconds=%.3lf\n", t);
    QueryPerformanceFrequency(&f);
    double p = queryperfcounter(X, sparse_W, B, Y, M, N, K, nonZero, f);
    printf("# Windows QueryPerformanceCounter() function\n");
    printf("queryperfcounter_cycles=%.0lf\n", p);
    printf("queryperfcounter_seconds=%.8lf\n", p/(double)f.QuadPart);
    printf("queryperfcounter_freq_mhz=%.2lf\n", (double)f.QuadPart/1000);
#endif

#ifdef __aarch64__
    double v = rdvct(X, sparse_W, B, Y, M, N, K, nonZero);
    printf("# VCT instruction\n");
    printf("vct_cycles=%.0lf\n", v);
    printf("vct_seconds=%.8lf\n", v/(get_vct_freq()));
    printf("vct_freq_mhz=%.2lf\n", (get_vct_freq())/1e6);

#ifdef PMU
    struct performance_counters p = rdpmu(X, sparse_W, B, Y, M, N, K);
    printf("# PMU instruction\n");
    printf("pmu_cycles=%.0lf\n", p.cycles);
    printf("pmu_seconds=%.8lf\n", p.cycles/(FREQUENCY));
    printf("pmu_freq_mhz=%.2lf\n", (FREQUENCY)/1e6);
    printf("pmu_instructions=%.0lf\n", p.instructions);
    printf("pmu_branches=%.0lf\n", p.branches);
    printf("pmu_branch_misses=%.0lf\n", p.branch_misses);
    printf("pmu_ipc=%.2lf\n", p.instructions / p.cycles);
#endif

#endif
    return 0;
}

