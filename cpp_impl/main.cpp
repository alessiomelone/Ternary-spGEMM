#include "perf.h"
#include "common.h"
#include "sparseUtils.h"
#include "comp.h"

std::vector<comp_func> userFuncs;
std::vector<std::string> funcNames;
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
            sparseGEMM_csc_base_impl<float>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "sparseGEMM_csc_base");

    add_function(
        [sf_csc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            sparseGEMM_csc_unrolled_impl<float, 16>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "sparseGEMM_csc_unrolled_16");

    add_function(
        [sf_csr](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            sparseGEMM_csr_base_impl<float>(X_arg, *sf_csr, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "sparseGEMM_csr_base");

    add_function(
        [sf_csr](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            sparseGEMM_csr_unrolled_impl<float, 16>(X_arg, *sf_csr, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "sparseGEMM_csr_unrolled_16");

    add_function(
        [sf_csc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            CSC_base<float>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "CSC_base");
    // add_function(
    //     [sf_ccsc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         CCSC_base<float>(X_arg, *sf_ccsc, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "CSC_base_testing");

    add_function(
        [sf_ccsc_data](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            CCSC_base<float>(X_arg, *sf_ccsc_data, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "CCSC_base");

    add_function(
        [sf_tcsr](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            TCSR_base<float>(X_arg, *sf_tcsr, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "TCSR_base");

    add_function(
        [sf_tcsr](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            TCSR_unrolled<float, 12>(X_arg, *sf_tcsr, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "TCSR_unrolled-12");

    add_function(
        [sf_tcsc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            TCSC_base<float>(X_arg, *sf_tcsc, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "TCSC_base");

    add_function(
        [sf_tcsc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            TCSC_unrolled<float, 12>(X_arg, *sf_tcsc, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "TCSC_unrolled-12");

    add_function(
        [sf_tcsc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            TCSC_unrolled_tiled<float, 12, 32, 32>(X_arg, *sf_tcsc, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "TCSC_unrolled_tiled-12-32-32");

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
        comp_func func = userFuncs[i_loop];


        Y_main.insert(Y_main.end(), 10, 0); // extend Y so we can modify unused pad values without bounds checking
        X_main.insert(X_main.end(), 10, 0); // extend Y so we can modify unused pad values without bounds checking
        comp_func func = userFuncs[i_loop]; // func is std::function

        Y_main.resize(Y_main.size() - 10);

        // Call the std::function directly. Sparse data is captured in the lambda.
        func(X_main.data(), B_main.data(), Y_main.data(), M, N, K);

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
