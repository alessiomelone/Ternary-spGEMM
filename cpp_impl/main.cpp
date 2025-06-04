#include "perf.h"
#include "common.h"
#include "sparseUtils.h"
#include "comp.h"
#include "comp_prelu.h"

#define BLOCK_SIZE 512
#define UNROLL_FACTOR 12

#define BENCHMARK_FUNCTION_NAME "BaseTCSC"

std::vector<comp_func> userFuncs;
std::vector<std::string> funcNames;
int numFuncs = 0;

// PrelU function infrastructure
std::vector<comp_func_prelu> userFuncs_prelu;
std::vector<std::string> funcNames_prelu;
int numFuncs_prelu = 0;

void add_function(comp_func f, std::string name)
{
    userFuncs.push_back(f);
    funcNames.emplace_back(name);
    numFuncs++;
}

void add_prelu_function(comp_func_prelu f, std::string name)
{
    userFuncs_prelu.push_back(f);
    funcNames_prelu.emplace_back(name);
    numFuncs_prelu++;
}

int main(int argc, char **argv)
{
    std::cout << "Starting program. ";
    float perf_val;
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
    bool check_correctness = false;
    if (argc > 9 && std::string(argv[9]) == "-correctness")
    {
        check_correctness = true;
    }

    // Generate sparse matrix to be converted
    std::vector<int> W_raw = generateSparseMatrix<int>(K, N, nonZero, false);

    // Initialize one instance per format
    auto sf_csc = std::make_shared<TCSC>(W_raw.data(), K, N);
    auto sf_csr = std::make_shared<TCSR>(W_raw.data(), K, N);
    // auto sf_blocked = std::make_shared<BlockedTCSC>(W_raw.data(), K, N);

    auto sf_vec_csc = std::make_shared<VectorTCSC>(W_raw.data(), K, N);

    auto sf_blocked = std::make_shared<BlockedTCSC<BLOCK_SIZE>>(W_raw.data(), K, N);
    // NOTE: BaseInterleavedBlockedTCSC and UnrolledInterleavedBlockedTCSC use different overloaded constructors
    auto sf_blocked_interleaved = std::make_shared<InterleavedBlockedTCSC<BLOCK_SIZE>>(W_raw.data(), K, N);
    auto sf_blocked_interleaved_unr = std::make_shared<InterleavedBlockedTCSC<BLOCK_SIZE>>(W_raw.data(), K, N, UNROLL_FACTOR);

    auto sf_interleaved = std::make_shared<InterleavedTCSC>(W_raw.data(), K, N);

    add_function(
        [sf_csc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            BaseTCSC<float>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "BaseTCSC");

    // add_function(
    //     [sf_blocked](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         BaseBlockedTCSC<float>(X_arg, *sf_blocked, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "BaseBlockedTCSC");

    // add_function(
    //     [sf_interleaved](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         BaseInterleavedTCSC<float>(X_arg, *sf_interleaved, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "BaseInterleavedTCSC");

    //     add_function(
    //     [sf_interleaved](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         UnrolledInterleavedTCSC<float, 12>(X_arg, *sf_interleaved, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "UnrolledInterleavedTCSC");

    // add_function(
    //     [sf_blocked_interleaved](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         BaseInterleavedBlockedTCSC<float>(X_arg, *sf_blocked_interleaved, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "BaseInterleavedBlockedTCSC");

    // add_function(
    //     [sf_csc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         UnrolledTCSC<float, 12>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "UnrolledTCSC_" + std::to_string(12));

    // add_function(
    //     [sf_csc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         DoubleUnrolledTCSC<float, 8, 4>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "DoubleUnrolledTCSC_K8_M4");

    add_function(
        [sf_csc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            DoubleUnrolledTCSC<float, 4, 4>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "DoubleUnrolledTCSC_K4_M4");

    // add_function(
    //     [sf_blocked_interleaved_unr](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         UnrolledInterleavedBlockedTCSC<float, BLOCK_SIZE, UNROLL_FACTOR>(X_arg, *sf_blocked_interleaved_unr, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "UnrolledInterleavedBlockedTCSC");

    // add_function(
    //     [sf_blocked](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         UnrolledBlockedTCSC<float, BLOCK_SIZE, 12>(X_arg, *sf_blocked, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "UnrolledBlockedTCSC_" + std::to_string(UNROLL_FACTOR));

    // add_function(
    //     [sf_csc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         UnrolledSimultaneousTCSC<float, 12>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "UnrolledSimultaneousTCSC_12");

    // add_function(
    //     [sf_csc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         UnrolledTCSC<float, UNROLL_FACTOR>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "UnrolledTCSC_" + std::to_string(UNROLL_FACTOR));

    // add_function(
    //     [sf_csc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         NeonTCSCHorizontalSimple<float>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "NeonTCSCHorizontalSimple");

    // add_function(
    //     [sf_vec_csc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         NeonTCSCVertical<float>(X_arg, *sf_vec_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "NeonTCSCVertical");

    // // Register PReLU function
    // add_prelu_function(
    //     [sf_csc](float *X_arg, float *B_arg, float *alpha_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         BaseTCSC_PreLU<float>(X_arg, *sf_csc, B_arg, alpha_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "BaseTCSC_PreLU");

    if (numFuncs == 0 && numFuncs_prelu == 0)
    {
        std::cout << std::endl;
        std::cout << "No functions registered - nothing for driver to do" << std::endl;
        std::cout << "Register functions in main.cpp by creating lambdas and using add_function or add_prelu_function." << std::endl;
        return 0;
    }

    std::cout << numFuncs << " regular functions and " << numFuncs_prelu << " PrelU functions registered." << std::endl;

    std::vector<float> X_main = initX<float>(M * K, 512);
    std::vector<float> W_FP32_main(W_raw.begin(), W_raw.end());
    std::vector<float> B_main(N, 2);
    std::vector<float> alpha_main(N, 0.1); // Alpha values for PrelU
    std::vector<float> Y_main(M * N, 0);
    std::vector<float> refY_main(M * N, 0);
    std::vector<float> refY_prelu_main(M * N, 0);

    // Compute reference for regular functions
    GEMM(X_main.data(), W_FP32_main.data(), B_main.data(), refY_main.data(), M, N, K);

    // Compute reference for PrelU functions
    GEMM_PreLU(X_main.data(), W_FP32_main.data(), B_main.data(), alpha_main.data(), refY_prelu_main.data(), M, N, K);

    if (check_correctness)
    {
        // Test regular functions
        for (i_loop = 0; i_loop < numFuncs; i_loop++)
        {
            fill(Y_main.begin(), Y_main.end(), 0);

            comp_func func = userFuncs[i_loop];
            func(X_main.data(), B_main.data(), Y_main.data(), M, N, K);

            if (compare_results(Y_main.data(), refY_main.data(), M, N))
            {
                std::cout << "Test case " << funcNames[i_loop] << " passed!" << std::endl;
            }
            else
            {
                std::cout << "Test case " << "\x1b[31m" << funcNames[i_loop] << " failed!" << "\x1b[0m" << std::endl;
                std::cout << "\n\n  Please fix the failing fn or comment out the invocaton from main.cpp.\n\nExiting...\n\n"
                          << std::endl;
                exit(1);
            }
        }

        // Test PrelU functions
        for (i_loop = 0; i_loop < numFuncs_prelu; i_loop++)
        {
            fill(Y_main.begin(), Y_main.end(), 0);

            comp_func_prelu func = userFuncs_prelu[i_loop];
            func(X_main.data(), B_main.data(), alpha_main.data(), Y_main.data(), M, N, K);

            if (compare_results(Y_main.data(), refY_prelu_main.data(), M, N))
            {
                std::cout << "Test case " << funcNames_prelu[i_loop] << " passed!" << std::endl;
            }
            else
            {
                std::cout << "Test case " << "\x1b[31m" << funcNames_prelu[i_loop] << " failed!" << "\x1b[0m" << std::endl;
                std::cout << "\n\n  Please fix the failing fn or comment out the invocaton from main.cpp.\n\nExiting...\n\n"
                          << std::endl;
                exit(1);
            }
        }
    }

    float base_cycles = 0;

    // Benchmark regular functions
    for (i_loop = 0; i_loop < numFuncs; i_loop++)
    {
        perf_val = perf_test(userFuncs[i_loop], M, K, N, nonZero);
        std::cout << "\nRunning: " << "\x1b[31m" << funcNames[i_loop] << "\x1b[0m" << std::endl;
        std::cout << perf_val << " cycles" << std::endl;
        if (funcNames[i_loop] == BENCHMARK_FUNCTION_NAME)
        {
            base_cycles = perf_val;
        }
        std::cout << "Speedup is: " << "\x1b[32m" << base_cycles / perf_val << "\x1b[0m" << std::endl;
#ifdef INSTRUMENTATION_RUN
        std::cout << "Flops: " << getTotalFlops() << std::endl;
        std::cout << "Performance: " << (float)getTotalFlops() / perf_val << " flops/cycle" << std::endl;
        float total_bytes = sizeof(float) * ((float)(M * K + M * N + N)) + (float)getDataStructureSizeInBytes();
        std::cout << "Total Input Size: " << (int)total_bytes << " Bytes" << std::endl;
        std::cout << "Operational Intensity: " << (float)getTotalFlops() / total_bytes << " Flops/Byte" << std::endl;
        std::cout << "Data Structure Size: " << getDataStructureSizeInBytes() << " Bytes" << std::endl;
#endif
    }

    // Benchmark PrelU functions
    float base_cycles_prelu = 0;
    for (i_loop = 0; i_loop < numFuncs_prelu; i_loop++)
    {
        perf_val = perf_test_prelu(userFuncs_prelu[i_loop], M, K, N, nonZero);
        std::cout << "\nRunning: " << "\x1b[31m" << funcNames_prelu[i_loop] << "\x1b[0m" << std::endl;
        std::cout << perf_val << " cycles" << std::endl;
        if (funcNames_prelu[i_loop] == "BaseTCSC_PreLU")
        {
            base_cycles_prelu = perf_val;
        }
        std::cout << "Speedup is: " << "\x1b[32m" << base_cycles_prelu / perf_val << "\x1b[0m" << std::endl;
#ifdef INSTRUMENTATION_RUN
        std::cout << "Flops: " << getTotalFlops() << std::endl;
        std::cout << "Performance: " << (float)getTotalFlops() / perf_val << " flops/cycle" << std::endl;
        float total_bytes = sizeof(float) * ((float)(M * K + M * N + N + N)) + (float)getDataStructureSizeInBytes(); // +N for alpha
        std::cout << "Total Input Size: " << (int)total_bytes << " Bytes" << std::endl;
        std::cout << "Operational Intensity: " << (float)getTotalFlops() / total_bytes << " Flops/Byte" << std::endl;
        std::cout << "Data Structure Size: " << getDataStructureSizeInBytes() << " Bytes" << std::endl;
#endif
    }

    return 0;
}
