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
#include <iomanip>

#include "perf.h"
#include "common.h"
#include "SparseGEMM.h"
#include "data_structures/CompressedCSC.h"
#include "data_structures/TCSRMatrix.h"

using namespace std;

// --- Prototypes for implementations in comp.cpp ---
// These are now declarations of explicitly instantiated templates
template <typename T>
void CSC_base(T *X, const SparseFormat &W_csc, T *b, T *Y, int M, int N, int K);
template <typename T>
void CCSC_base(T *X, const CompressedCSC &W_csc, T *b, T *Y, int M, int N, int K);
template <typename T>
void TCSR_base(T *X, const TCSRMatrix &W_tcsr, T *b, T *Y, int M, int N, int K);
template <typename T, int UNROLL_FACTOR> // Provide default for UNROLL_FACTOR if used in declaration
void CSC_unrolled(T *X, const SparseFormat &W_csc, T *b, T *Y, int M, int N, int K);
// --- End Prototypes ---

vector<comp_func> userFuncs; // This is now vector<std::function<...>>
vector<string> funcNames;
int numFuncs = 0;

/*
 * Registers a user function to be tested by the driver program. Registers a
 * string description of the function as well
 */
void add_function(comp_func f, string name) // Signature of add_function itself doesn't change
{
    userFuncs.push_back(f);
    funcNames.emplace_back(name);
    numFuncs++;
}

// Example for a new format (ensure CustomMixedTypeFormat is defined and this prototype matches comp.cpp)
/*
template <typename T>
void sparseGEMM_custom_mixed_impl(T *X, const CustomMixedTypeFormat& W_custom, T *b, T *Y, int M, int N, int K);
*/
// --- End Prototypes ---

int main(int argc, char **argv)
{
    double perf_val; // Renamed perf
    int i_loop;      // Renamed i

    int M = 0, K = 0, N = 0, nonZero = 0;

    if (argc < 9)
    {
        fprintf(stderr, "Usage: %s -M <int> -K <int> -N <int> -s <int>\n", argv[0]);
        return 1;
    }

    M = atoi(argv[2]);
    K = atoi(argv[4]);
    N = atoi(argv[6]);
    nonZero = atoi(argv[8]);
    if (M <= 0 || K <= 0 || N <= 0 || nonZero <= 0)
    {
        fprintf(stderr, "ERROR: All dimensions must be positive integers.\n");
        fprintf(stderr, "Usage: %s -M <int> -K <int> -N <int> -s <int>\n", argv[0]);
        return 1;
    }

    // Create data structures that will be captured by lambdas.
    // Use std::shared_ptr to manage their lifetime.
    vector<int> W_raw = generateSparseMatrix<int>(K, N, nonZero, false, 0); // For SparseFormat
    auto sf_csc_data = std::make_shared<SparseFormat>(W_raw.data(), K, N);
    auto sf_ccsc_data = std::make_shared<CompressedCSC>(W_raw.data(), K, N);
    auto sf_tcsr_data = std::make_shared<TCSRMatrix>(W_raw.data(), K, N);

    // Example for a custom data structure:
    // You would need to generate or prepare data for CustomMixedTypeFormat here
    // auto custom_mixed_data = std::make_shared<CustomMixedTypeFormat>(/* constructor args for custom format */);

    // --- Register functions using lambdas ---
    add_function(
        [sf_csc_data](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            CSC_base<float>(X_arg, *sf_csc_data, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "CSC_base");

    // add_function(
    //     [sf_ccsc_data](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         CCSC_base<float>(X_arg, *sf_ccsc_data, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "CCSC_base");

    add_function(
        [sf_tcsr_data](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            TCSR_base<float>(X_arg, *sf_tcsr_data, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "TCSR_base");

    // add_function(
    //     [sf_csc_data](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg) {
    //         CSC_unrolled<float, 2>(X_arg, *sf_csc_data, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //         // Note: You can vary UNROLL_FACTOR here or make it part of the name if you test multiple unroll factors
    //     },
    //     "CSC_unrolled-16"
    // );

    // Example registration for a custom format:
    /*
    add_function(
        [custom_mixed_data](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg) {
            sparseGEMM_custom_mixed_impl<float>(X_arg, *custom_mixed_data, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "sparseGEMM_custom_mixed"
    );
    */
    // --- End Register functions ---

    if (numFuncs == 0)
    {
        cout << endl;
        cout << "No functions registered - nothing for driver to do" << endl;
        // The old message about register_funcs() in comp.cpp is no longer relevant here.
        cout << "Register functions in main.cpp by creating lambdas and using add_function." << endl;
        return 0;
    }
    cout << numFuncs << " functions registered." << endl;

    vector<float> X_main = initX<float>(M * K, 512);       // Renamed X
    vector<float> W_FP32_main(W_raw.begin(), W_raw.end()); // Renamed W_FP32
    vector<float> B_main(N, 0);                            // Renamed B
    vector<float> Y_main(M * N, 0);                        // Renamed Y
    vector<float> refY_main(M * N, 0);                     // Renamed refY

    GEMM(X_main.data(), W_FP32_main.data(), B_main.data(), refY_main.data(), M, N, K);

    for (i_loop = 0; i_loop < numFuncs; i_loop++)
    {
        fill(Y_main.begin(), Y_main.end(), 0);

        Y_main.insert(Y_main.end(), 10, 0); // extend Y so we can modify unused pad values without bounds checking
        comp_func func = userFuncs[i_loop]; // func is std::function

        Y_main.resize(Y_main.size() - 10);

        // Call the std::function directly. Sparse data is captured in the lambda.
        func(X_main.data(), B_main.data(), Y_main.data(), M, N, K);

        if (compare_results(Y_main.data(), refY_main.data(), M, N))
        {
            cout << "Test case " << funcNames[i_loop] << " passed!" << endl;
        }
        else
        {
            cout << "Test case " << funcNames[i_loop] << " failed!" << endl;
        }
    }

    // No need to explicitly clear X_main, W_raw, etc. if they are std::vectors,
    // they will clean up when they go out of scope.
    // The shared_ptrs (sf_csr_data, etc.) will also manage their memory.

    for (i_loop = 0; i_loop < numFuncs; i_loop++)
    {
        // perf_test takes the std::function.
        // Data for X, B, Y for the benchmark run is created inside perf_test.
        perf_val = perf_test(userFuncs[i_loop], M, K, N, nonZero);
        cout << endl
             << "Running: " << funcNames[i_loop] << endl;
        cout << perf_val << " cycles" << endl;
        // Ensure floating point division for performance metric
        cout << "Performance: " << static_cast<double>(M * N) * (1.0 + static_cast<double>(K) / nonZero) / perf_val << " flops/cycle" << endl;
    }

    return 0;
}
