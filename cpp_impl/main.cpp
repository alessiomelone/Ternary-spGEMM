#include "perf.h"
#include "common.h"
#include "sparseUtils.h"
#include "comp.h"

<<<<<<< HEAD
std::vector<comp_func> userFuncs;
std::vector<std::string> funcNames;
=======
using namespace std;

// --- Prototypes for implementations in comp.cpp ---
// These are now declarations of explicitly instantiated templates
template <typename T>
void CSC_base(T *X, const SparseFormat &W_csc, T *b, T *Y, int M, int N, int K);
template <typename T>
void CSC_base_testing(T *X, const SparseFormat &W_csc, T *b, T *Y, int M, int N, int K);
template <typename T>
void CCSC_base(T *X, const CompressedCSC &W_csc, T *b, T *Y, int M, int N, int K);
template <typename T>
void TCSR_base(T *X, const TCSRMatrix &W_tcsr, T *b, T *Y, int M, int N, int K);
template <typename T>
void TCSC_base(T *X, const TCSCMatrix &W_tcsc, T *b, T *Y, int M, int N, int K);
template <typename T, int UNROLL_FACTOR> // Provide default for UNROLL_FACTOR if used in declaration
void CSC_unrolled(T *X, const SparseFormat &W_csc, T *b, T *Y, int M, int N, int K);
template <typename T, int UNROLL_FACTOR>
void TCSR_unrolled(T *X, const TCSRMatrix &W_tcsr, T *b, T *Y, int M, int N, int K);
template <typename T, int UNROLL_FACTOR>
void TCSC_unrolled(T *X, const TCSCMatrix &W_tcsc, T *b, T *Y, int M, int N, int K);

template <typename T, int UNROLL_FACTOR, int TILE_M, int TILE_N>
void TCSC_unrolled_tiled(T *X_arg, const TCSCMatrix &W_tcsc, T *B_arg, T *Y_arg,
                         int M_dim, int N_dim, int K_dim);
// --- End Prototypes ---

vector<comp_func> userFuncs; // This is now vector<std::function<...>>
vector<string> funcNames;
>>>>>>> 00c090c (Add TCSR and TCSC base function prototypes; refactor CSC_base and CCSC_base implementations)
int numFuncs = 0;

void add_function(comp_func f, std::string name)
{
    userFuncs.push_back(f);
    funcNames.emplace_back(name);
    numFuncs++;
}

int main(int argc, char **argv)
{
    std::cout << "Starting program. ";
    double perf_val;
    int i_loop;

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

    // Generate sparse matrix to be converted
    std::vector<int> W_raw = generateSparseMatrix<int>(K, N, nonZero, false);

    // Initialize one instance per format
    auto sf_csc = std::make_shared<SparseFormatCSC>(W_raw.data(), K, N);
    auto sf_csr = std::make_shared<SparseFormatCSR>(W_raw.data(), K, N);
    auto sf_ccsc = std::make_shared<CompressedCSC>(W_raw.data(), K, N);
    auto sf_tcsr = std::make_shared<TCSRMatrix>(W_raw.data(), K, N);
    auto sf_tcsc = std::make_shared<TCSCMatrix>(W_raw.data(), K, N);

    // Register functions using the shared instances

    add_function(
        [sf_csc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            CSC_base<float>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "CSC_base");
    // add_function(
<<<<<<< HEAD
    //     [sf_ccsc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         CCSC_base<float>(X_arg, *sf_ccsc, B_arg, Y_arg, M_arg, N_arg, K_arg);
=======
    //     [sf_csc_data](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         CSC_base_testing<float>(X_arg, *sf_csc_data, B_arg, Y_arg, M_arg, N_arg, K_arg);
>>>>>>> 00c090c (Add TCSR and TCSC base function prototypes; refactor CSC_base and CCSC_base implementations)
    //     },
    //     "CSC_base_testing");

    add_function(
        [sf_ccsc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            CCSC_base<float>(X_arg, *sf_ccsc, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "CCSC_base");

    add_function(
        [sf_tcsr](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            TCSR_base<float>(X_arg, *sf_tcsr, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "TCSR_base");

<<<<<<< HEAD
    add_function(
        [sf_tcsr](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            TCSR_unrolled_tiled<float, 8, 8>(X_arg, *sf_tcsr, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
<<<<<<< HEAD
<<<<<<< HEAD
        "TCSR_unrolled-12");
=======
    // add_function(
    //     [sf_tcsr_unrolled_data](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         TCSR_unrolled<float, 12>(X_arg, *sf_tcsr_unrolled_data, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "TCSR_unrolled-12");
>>>>>>> 00c090c (Add TCSR and TCSC base function prototypes; refactor CSC_base and CCSC_base implementations)
=======
        "TCSR_unrolled_tiled-12-16-16");
>>>>>>> f2770dd (TCSC tiled is 1.5)
=======
        "TCSR_unrolled_tiled-12-4-4");
>>>>>>> 84e2726 (Tiling for TCSR improved but still needs works)

    add_function(
        [sf_tcsc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            TCSC_base<float>(X_arg, *sf_tcsc, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "TCSC_base");

<<<<<<< HEAD
<<<<<<< HEAD
    add_function(
        [sf_tcsc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            TCSC_unrolled<float, 12>(X_arg, *sf_tcsc, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "TCSC_unrolled-12");
=======
    // add_function(
    //     [sf_tcsc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         TCSC_unrolled<float, 12>(X_arg, *sf_tcsc, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "TCSC_unrolled-12");
>>>>>>> f2770dd (TCSC tiled is 1.5)

    add_function(
        [sf_tcsc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            TCSC_unrolled_tiled<float, 12, 8, 8>(X_arg, *sf_tcsc, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
<<<<<<< HEAD
        "TCSC_unrolled_tiled-12-32-32");
=======
    // add_function(
    //     [sf_tcsc_unrolled_data](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         TCSC_unrolled<float, 12>(X_arg, *sf_tcsc_unrolled_data, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "TCSC_unrolled-12");

    // add_function(
    //     [sf_tcsc_unrolled_tiled_data](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         TCSC_unrolled_tiled<float, 12, 32, 32>(X_arg, *sf_tcsc_unrolled_tiled_data, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "TCSC_unrolled_tiled-12-32-32");
>>>>>>> 00c090c (Add TCSR and TCSC base function prototypes; refactor CSC_base and CCSC_base implementations)
=======
        "TCSC_unrolled_tiled-12-16-16");
>>>>>>> f2770dd (TCSC tiled is 1.5)

    if (numFuncs == 0)
    {
        std::cout << std::endl;
        std::cout << "No functions registered - nothing for driver to do" << std::endl;
        std::cout << "Register functions in main.cpp by creating lambdas and using add_function." << std::endl;
        return 0;
    }
    std::cout << numFuncs << " functions registered." << std::endl;

    std::vector<float> X_main = initX<float>(M * K, 512);
    std::vector<float> W_FP32_main(W_raw.begin(), W_raw.end());
    std::vector<float> B_main(N, 2);
    std::vector<float> Y_main(M * N, 0);
    std::vector<float> refY_main(M * N, 0);

    GEMM(X_main.data(), W_FP32_main.data(), B_main.data(), refY_main.data(), M, N, K);

    for (i_loop = 0; i_loop < numFuncs; i_loop++)
    {
        fill(Y_main.begin(), Y_main.end(), 0);

<<<<<<< HEAD
<<<<<<< HEAD
=======
=======
        comp_func func = userFuncs[i_loop]; // func is std::function

>>>>>>> 6073638 (merge and fix bugs)
        Y_main.insert(Y_main.end(), 10, 0); // extend Y so we can modify unused pad values without bounds checking
        X_main.insert(X_main.end(), 10, 0); // extend Y so we can modify unused pad values without bounds checking

        // Call the std::function directly. Sparse data is captured in the lambda.
>>>>>>> 00c090c (Add TCSR and TCSC base function prototypes; refactor CSC_base and CCSC_base implementations)
        func(X_main.data(), B_main.data(), Y_main.data(), M, N, K);

        Y_main.resize(Y_main.size() - 10);

        if (compare_results(Y_main.data(), refY_main.data(), M, N))
        {
            std::cout << "Test case " << funcNames[i_loop] << " passed!" << std::endl;
        }
        else
        {
            std::cout << "Test case " << funcNames[i_loop] << " failed!" << std::endl;
        }
    }

    for (i_loop = 0; i_loop < numFuncs; i_loop++)
    {
        perf_val = perf_test(userFuncs[i_loop], M, K, N, nonZero);
        std::cout << std::endl
                  << "Running: " << funcNames[i_loop] << std::endl;
        std::cout << perf_val << " cycles" << std::endl;
        std::cout << "Performance: " << static_cast<double>(M * N) * (1.0 + static_cast<double>(K) / nonZero) / perf_val << " flops/cycle" << std::endl;
    }

    return 0;
}
