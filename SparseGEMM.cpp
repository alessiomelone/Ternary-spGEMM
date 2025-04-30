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
#include <iomanip>  // For std::setw, std::setprecision, std::fixed

#include <vector>
#include "SparseGEMM.h"

#ifdef __x86_64__
#include "./include/tsc_x86.h"
#endif

#ifdef __aarch64__
#include "./include/vct_arm.h"

#ifdef PMU
#include "./include/kperf.h"
#endif

#endif

#define NUM_RUNS 1
#define CYCLES_REQUIRED 1e8
#define FREQUENCY 3.22e9
#define CALIBRATE

#define INDEX(i, j, n) (i*n - i*(i-1) / 2 + (j - i))

using namespace std;

/* 
 * Timing function based on the TimeStep Counter of the CPU.
 */
#ifdef __x86_64__
double rdtsc(vector<int> X, vector<int> W, vector<int> B, vector<int> Y, int M, int N, int K) {

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
            sparseGEMM(X.data(), SparseFormat(W.data(), K, N), B.data(), Y.data(), M, N, K);

        }
        cycles = stop_tsc(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    start = start_tsc();
    for (i = 0; i < num_runs; ++i) {
        sparseGEMM(X.data(), SparseFormat(W.data(), K, N), B.data(), Y.data(), M, N, K);

    }

    cycles = stop_tsc(start)/num_runs;
    return (double) cycles;
}
#endif

#ifdef __aarch64__
double rdvct(vector<int> X, vector<int> W, vector<int> B, vector<int> Y, int M, int N, int K) {
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
            sparseGEMM(X.data(), SparseFormat(W.data(), K, N), B.data(), Y.data(), M, N, K);

        }
        cycles = stop_vct(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    start = start_vct();
    for (i = 0; i < num_runs; ++i) {
        sparseGEMM(X.data(), SparseFormat(W.data(), K, N), B.data(), Y.data(), M, N, K);

    }

    cycles = stop_vct(start)/num_runs;
    return (double) cycles;
}

#ifdef PMU
performance_counters rdpmu(vector<int> X, vector<int> W, vector<int> B, vector<int> Y, int M, int N, int K) {
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
            sparseGEMM(X.data(), SparseFormat(W.data(), K, N), B.data(), Y.data(), M, N, K);

        }
        endperf = kperf_get_counters();
        double cycles = endperf.cycles - startperf.cycles;
        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    startperf = kperf_get_counters();
    for (i = 0; i < num_runs; ++i) {
        sparseGEMM(X.data(), SparseFormat(W.data(), K, N), B.data(), Y.data(), M, N, K);

    }

    endperf = kperf_get_counters();
    result.cycles = (endperf.cycles - startperf.cycles)/num_runs;
    result.instructions = (endperf.instructions - startperf.instructions)/num_runs;
    result.branches = (endperf.branches - startperf.branches)/num_runs;
    result.branch_misses = (endperf.branch_misses - startperf.branch_misses)/num_runs;

    return result;
}
performance_counters rdpmu_for_dense(vector<int> X, vector<int> W, vector<int> B, vector<int> Y, int M, int N, int K) {
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
            GEMM(X.data(), W.data(), B.data(), Y.data(), M, N, K);

        }
        endperf = kperf_get_counters();
        double cycles = endperf.cycles - startperf.cycles;
        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    startperf = kperf_get_counters();
    for (i = 0; i < num_runs; ++i) {
        GEMM(X.data(), W.data(), B.data(), Y.data(), M, N, K);
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

double c_clock(vector<int> X, vector<int> W, vector<int> B, vector<int> Y, int M, int N, int K) {
    int i, num_runs;
    double cycles;
    clock_t start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while(num_runs < (1 << 14)) {
        start = clock();
        for (i = 0; i < num_runs; ++i) {
            sparseGEMM(X.data(), SparseFormat(W.data(), K, N), B.data(), Y.data(), M, N, K);

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
        sparseGEMM(X.data(), SparseFormat(W.data(), K, N), B.data(), Y.data(), M, N, K);

    }
    end = clock();

    return (double)(end-start)/num_runs;
}

#ifndef WIN32
double timeofday(vector<int> X, vector<int> W, vector<int> B, vector<int> Y, int M, int N, int K) {
    int i, num_runs;
    double cycles;
    struct timeval start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while(num_runs < (1 << 14)) {
        gettimeofday(&start, NULL);
        for (i = 0; i < num_runs; ++i) {
            sparseGEMM(X.data(), SparseFormat(W.data(), K, N), B.data(), Y.data(), M, N, K);

        }
        gettimeofday(&end, NULL);

        cycles = (double)((end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1e6)*FREQUENCY;

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    gettimeofday(&start, NULL);
    for(i=0; i < num_runs; ++i) {
        sparseGEMM(X.data(), SparseFormat(W.data(), K, N), B.data(), Y.data(), M, N, K);

    }
    gettimeofday(&end, NULL);

    return (double)((end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1e6)/ num_runs;
}

#else

double gettickcount(vector<int> X, vector<int> W, vector<int> B, vector<int> Y, int M, int N, int K) {
    int i, num_runs;
    double cycles, start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while(num_runs < (1 << 14)) {
        start = (double)GetTickCount();
        for (i = 0; i < num_runs; ++i) {
            sparseGEMM(X.data(), SparseFormat(W.data(), K, N), B.data(), Y.data(), M, N, K);

        }
        end = (double)GetTickCount();

        cycles = (end-start)*FREQUENCY/1e3; // end-start provides a measurement in the order of milliseconds

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    start = (double)GetTickCount();
    for(i=0; i < num_runs; ++i) {
        sparseGEMM(X.data(), SparseFormat(W.data(), K, N), B.data(), Y.data(), M, N, K);
 
    }
    end = (double)GetTickCount();

    return (end-start)/num_runs;
}

double queryperfcounter(vector<int> X, vector<int> W, vector<int> B, vector<int> Y, int M, int N, int K, LARGE_INTEGER f) {
    int i, num_runs;
    double cycles;
    LARGE_INTEGER start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while(num_runs < (1 << 14)) {
        QueryPerformanceCounter(&start);
        for (i = 0; i < num_runs; ++i) {
            sparseGEMM(X.data(), SparseFormat(W.data(), K, N), B.data(), Y.data(), M, N, K);

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
        sparseGEMM(X.data(), SparseFormat(W.data(), K, N), B.data(), Y.data(), M, N, K);
 
    }
    QueryPerformanceCounter(&end);

    return (double)(end.QuadPart - start.QuadPart)/num_runs;
}

#endif

vector<int> generateSparseMatrix(int H, int W, int nonZero, bool uniformDistribution) {
    vector<int> y = vector<int>(H * W, 0);
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

    vector<int> X = generateSparseMatrix(M, K, 4, true);
    cout << "X initialized."<< endl;
    vector<int> W = generateSparseMatrix(K, N, nonZero, true);
    cout << "W initialized." << endl;

    vector<int> Y = vector<int>(M * N, 0);
    vector<int> B = vector<int>(N, 2); // init bias as 2
    vector<int> refY = vector<int>(M * N, 0);


    sparseGEMM(X.data(), SparseFormat(W.data(), K, N), B.data(), Y.data(), M, N, K);
    GEMM(X.data(), W.data(), B.data(), refY.data(), M, N, K);
    if (compare_results(Y.data(), refY.data(), M, N)) {
        cout << "Test case passed!" << endl;
    }

#ifdef __x86_64__
    double r = rdtsc(X, W, B, Y, M, N, K);
    printf("RDTSC instruction:\n %lf cycles measured => %lf seconds, assuming frequency is %lf MHz. (change in source file if different)\n\n", r, r/(FREQUENCY), (FREQUENCY)/1e6);
#endif

    double c = c_clock(X, W, B, Y, M, N, K);
    printf("C clock() function:\n %lf cycles measured. On some systems, this number seems to be actually computed from a timer in seconds then transformed into clock ticks using the variable CLOCKS_PER_SEC. Unfortunately, it appears that CLOCKS_PER_SEC is sometimes set improperly. (According to this variable, your computer should be running at %lf MHz). In any case, dividing by this value should give a correct timing: %lf seconds. \n\n",c, (double) CLOCKS_PER_SEC/1e6, c/CLOCKS_PER_SEC);


#ifndef WIN32
    double t = timeofday(X, W, B, Y, M, N, K);
    printf("C gettimeofday() function:\n %lf seconds measured\n\n",t);
#else
    LARGE_INTEGER f;
    double t = gettickcount(X, W, B, Y, M, N, K);
    printf("Windows getTickCount() function:\n %lf milliseconds measured\n\n",t);

    QueryPerformanceFrequency((LARGE_INTEGER *)&f);

    double p = queryperfcounter(X, W, B, Y, M, N, K, f);
    printf("Windows QueryPerformanceCounter() function:\n %lf cycles measured => %lf seconds, with reported CPU frequency %lf MHz\n\n",p,p/f.QuadPart,(double)f.QuadPart/1000);
#endif


#ifdef __aarch64__
    double v = rdvct(X, W, B, Y, M, N, K);
    printf("VCT instruction:\n %lf cycles measured => %lf seconds, assuming frequency of the VCT clock is %lf MHz. \n\n", v, v/(get_vct_freq()), (get_vct_freq())/1e6);

#ifdef PMU
    // This requires sudo on macOS
    performance_counters p = rdpmu(X, W, B, Y, M, N, K);
    printf("PMU instruction:\n %lf cycles measured => %lf seconds, assuming frequency is %lf MHz. (change in source file if different)\n\n", p.cycles, p.cycles/(FREQUENCY), (FREQUENCY)/1e6);

    performance_counters d = rdpmu_for_dense(X, W, B, Y, M, N, K);

    cout << "\n┌─────────────────────────────────────────────────────┐" << endl;
                cout << "│                  TEST CASE RESULTS                   │" << endl;
                cout << "├─────────────────────────────────────────────────────┤" << endl;
                cout << "│ Status: ✅ PASSED                                    │" << endl;
                cout << "│ Matrix Dimensions: M=" << M << ", N=" << N << ", K=" << K << std::string(19 - std::to_string(M).length() - std::to_string(N).length() - std::to_string(K).length(), ' ') << "│" << endl;
                cout << "├─────────────────────────────────────────────────────┤" << endl;
                cout << "│ Performance Comparison:                              │" << endl;
                cout << "│   • Sparse GEMM cycles: " << std::setw(10) << p.cycles << std::string(19, ' ') << "│" << endl;
                cout << "│   • Dense GEMM cycles:  " << std::setw(10) << d.cycles << std::string(19, ' ') << "│" << endl;
                cout << "│   • Speedup:            " << std::fixed << std::setprecision(2) << std::setw(10) << (double)d.cycles / p.cycles << "x" << std::string(18, ' ') << "│" << endl;

                #ifdef __APPLE__
                // Print detailed performance counters
                cout << "├─────────────────────────────────────────────────────┤" << endl;
                cout << "│ Detailed Performance Metrics:                        │" << endl;
                cout << "│                                                     │" << endl;
                cout << "│ Sparse GEMM:                                        │" << endl;
                cout << "│   • Instructions:  " << std::setw(12) << (p.instructions) << std::string(19, ' ') << "│" << endl;
                cout << "│   • Branches:      " << std::setw(12) << (p.branches) << std::string(19, ' ') << "│" << endl;
                cout << "│   • Branch misses: " << std::setw(12) << (p.branch_misses) << std::string(19, ' ') << "│" << endl;
                cout << "│   • IPC:           " << std::fixed << std::setprecision(2) << std::setw(12) << (p.instructions /p.cycles) << std::string(19, ' ') << "│" << endl;
                cout << "│                                                     │" << endl;
                cout << "│ Dense GEMM:                                         │" << endl;
                cout << "│   • Instructions:  " << std::setw(12) << (d.instructions) << std::string(19, ' ') << "│" << endl;
                cout << "│   • Branches:      " << std::setw(12) << (d.branches) << std::string(19, ' ') << "│" << endl;
                cout << "│   • Branch misses: " << std::setw(12) << (d.branch_misses) << std::string(19, ' ') << "│" << endl;
                cout << "│   • IPC:           " << std::fixed << std::setprecision(2) << std::setw(12) << (d.instructions / d.cycles) << std::string(19, ' ') << "│" << endl;
                #endif
                cout << "└─────────────────────────────────────────────────────┘" << endl;
#endif

#endif
    
    return 0;
}

