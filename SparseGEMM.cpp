#include "SparseGEMM.h"
#include <iomanip>  // For std::setw, std::setprecision, std::fixed

// Include appropriate timing header based on architecture
#ifdef __aarch64__
#include "include/vct_arm.h"
#else
#include "include/tsc_x86.h"
#endif

// Include kperf for detailed performance counters (macOS only)
#ifdef __APPLE__
#include "include/kperf.h"
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

int main() {
    // Initialize timing infrastructure
    #ifdef __APPLE__
    // Initialize kperf for detailed performance counters
    if (!kperf_init()) {
        cout << "Failed to initialize kperf. Continuing without detailed performance counters." << endl;
    }
    #endif

    // Smaller test cases
    int TEST_CASES = 2;
    int M[] = { 1, 16};
    int K[] = { 512, 1024 };
    int N[] = { 2048, 4096 };
    int nonZero[] = { 2, 4, 8, 16 }; // 1/2, 1/4, 1/8, 1/16 non-zero values

    vector<int> X = generateSparseMatrix(M[1], K[1], 4, true);
    cout << "X initialized."<< endl;
    vector<int> W = generateSparseMatrix(K[1], N[1], nonZero[3], true);
    cout << "W initialized." << endl;

    for (int m = 0; m < TEST_CASES; m++) {
        for (int n = 0; n < TEST_CASES; n++) {
            vector<int> Y = vector<int>(M[m] * N[n], 0);
            vector<int> B = vector<int>(N[n], 2);
            vector<int> refY = vector<int>(M[m] * N[n], 0);
            // Time the sparse implementation
            #ifdef __APPLE__
            // Get detailed performance counters before sparse GEMM
            struct performance_counters before_sparse = kperf_get_counters();
            #endif

            #ifdef __aarch64__
            TIMESTAMP start_sparse = start_vct();
            #else
            myInt64 start_sparse = start_tsc();
            #endif

            sparseGEMM(X.data(), SparseFormat(W.data(), K[n], N[n]), B.data(), Y.data(), M[m], N[n], K[n]);

            #ifdef __aarch64__
            TIMESTAMP cycles_sparse = stop_vct(start_sparse);
            #else
            myInt64 cycles_sparse = stop_tsc(start_sparse);
            #endif

            #ifdef __APPLE__
            // Get detailed performance counters after sparse GEMM
            struct performance_counters after_sparse = kperf_get_counters();
            #endif

            // Time the dense implementation
            #ifdef __APPLE__
            // Get detailed performance counters before dense GEMM
            struct performance_counters before_dense = kperf_get_counters();
            #endif

            #ifdef __aarch64__
            TIMESTAMP start_dense = start_vct();
            #else
            myInt64 start_dense = start_tsc();
            #endif

            GEMM(X.data(), W.data(), B.data(), refY.data(), M[m], N[n], K[n]);

            #ifdef __aarch64__
            TIMESTAMP cycles_dense = stop_vct(start_dense);
            #else
            myInt64 cycles_dense = stop_tsc(start_dense);
            #endif

            #ifdef __APPLE__
            // Get detailed performance counters after dense GEMM
            struct performance_counters after_dense = kperf_get_counters();
            #endif

            if (compare_results(Y.data(), refY.data(), M[m], N[n])) {
                cout << "\n┌─────────────────────────────────────────────────────┐" << endl;
                cout << "│                  TEST CASE RESULTS                   │" << endl;
                cout << "├─────────────────────────────────────────────────────┤" << endl;
                cout << "│ Status: ✅ PASSED                                    │" << endl;
                cout << "│ Matrix Dimensions: M=" << M[m] << ", N=" << N[n] << ", K=" << K[n] << std::string(19 - std::to_string(M[m]).length() - std::to_string(N[n]).length() - std::to_string(K[n]).length(), ' ') << "│" << endl;
                cout << "├─────────────────────────────────────────────────────┤" << endl;
                cout << "│ Performance Comparison:                              │" << endl;
                cout << "│   • Sparse GEMM cycles: " << std::setw(10) << cycles_sparse << std::string(19, ' ') << "│" << endl;
                cout << "│   • Dense GEMM cycles:  " << std::setw(10) << cycles_dense << std::string(19, ' ') << "│" << endl;
                cout << "│   • Speedup:            " << std::fixed << std::setprecision(2) << std::setw(10) << (double)cycles_dense / cycles_sparse << "x" << std::string(18, ' ') << "│" << endl;

                #ifdef __APPLE__
                // Print detailed performance counters
                cout << "├─────────────────────────────────────────────────────┤" << endl;
                cout << "│ Detailed Performance Metrics:                        │" << endl;
                cout << "│                                                     │" << endl;
                cout << "│ Sparse GEMM:                                        │" << endl;
                cout << "│   • Instructions:  " << std::setw(12) << (after_sparse.instructions - before_sparse.instructions) << std::string(19, ' ') << "│" << endl;
                cout << "│   • Branches:      " << std::setw(12) << (after_sparse.branches - before_sparse.branches) << std::string(19, ' ') << "│" << endl;
                cout << "│   • Branch misses: " << std::setw(12) << (after_sparse.branch_misses - before_sparse.branch_misses) << std::string(19, ' ') << "│" << endl;
                cout << "│   • IPC:           " << std::fixed << std::setprecision(2) << std::setw(12) << (after_sparse.instructions - before_sparse.instructions) /
                                    (after_sparse.cycles - before_sparse.cycles) << std::string(19, ' ') << "│" << endl;
                cout << "│                                                     │" << endl;
                cout << "│ Dense GEMM:                                         │" << endl;
                cout << "│   • Instructions:  " << std::setw(12) << (after_dense.instructions - before_dense.instructions) << std::string(19, ' ') << "│" << endl;
                cout << "│   • Branches:      " << std::setw(12) << (after_dense.branches - before_dense.branches) << std::string(19, ' ') << "│" << endl;
                cout << "│   • Branch misses: " << std::setw(12) << (after_dense.branch_misses - before_dense.branch_misses) << std::string(19, ' ') << "│" << endl;
                cout << "│   • IPC:           " << std::fixed << std::setprecision(2) << std::setw(12) << (after_dense.instructions - before_dense.instructions) /
                                    (after_dense.cycles - before_dense.cycles) << std::string(19, ' ') << "│" << endl;
                #endif
                cout << "└─────────────────────────────────────────────────────┘" << endl;
            }
        }
    }
}