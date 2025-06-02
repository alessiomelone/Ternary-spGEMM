#include "perf.h"
#include "common.h"
#include "sparseUtils.h"
#include "comp.h"

#define BLOCK_SIZE 512
#define UNROLL_FACTOR_IBTCSC 16

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

    auto sf_blocked = std::make_shared<BlockedTCSC<BLOCK_SIZE>>(W_raw.data(), K, N);
    auto sf_blocked_interleaved = std::make_shared<InterleavedBlockedTCSC<BLOCK_SIZE>>(W_raw.data(), K, N);

    auto sf_interleaved = std::make_shared<InterleavedTCSC>(W_raw.data(), K, N);

    add_function(
        [sf_csc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            BaseTCSC<float>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "BaseTCSC");

    add_function(
        [sf_csr](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            BaseTCSR<float>(X_arg, *sf_csr, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "BaseTCSR");

    add_function(
        [sf_blocked](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            BaseBlockedTCSC<float, BLOCK_SIZE>(X_arg, *sf_blocked, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "BaseBlockedTCSC");

    add_function(
        [sf_interleaved](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            BaseInterleavedTCSC<float>(X_arg, *sf_interleaved, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "BaseInterleavedTCSC");

    add_function(
        [sf_csc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            UnrolledModifiedTCSC<float, 12>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "UnrolledModifiedTCSC_12");

    add_function(
        [sf_csc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            UnrolledTCSC<float, 12>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "UnrolledTCSC_12");

    add_function(
        [sf_interleaved](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            UnrolledInterleavedTCSC<float, 8>(X_arg, *sf_interleaved, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "BaseInterleavedTCSC");

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

    if (check_correctness)
    {
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
    }

    float base_cycles = 0;
    for (i_loop = 0; i_loop < numFuncs; i_loop++)
    {
        perf_val = perf_test(userFuncs[i_loop], M, K, N, nonZero);
        std::cout << "\nRunning: " << "\x1b[31m" << funcNames[i_loop] << "\x1b[0m" << std::endl;
        std::cout << perf_val << " cycles" << std::endl;
        if (funcNames[i_loop] == "BaseTCSC")
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
        // std::cout << "Performance: " << static_cast<float>(M * N) * (1.0 + static_cast<float>(K) / nonZero) / perf_val << " flops/cycle" << std::endl;
    }

    return 0;
}
