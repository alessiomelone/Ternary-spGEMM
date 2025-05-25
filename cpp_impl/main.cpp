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
    bool check_correctness = false;
    if (argc > 9 && std::string(argv[9]) == "-correctness")
    {
        check_correctness = true;
    }

    // Generate sparse matrix to be converted
    std::vector<int> W_raw = generateSparseMatrix<int>(K, N, nonZero, false);

    // Initialize one instance per format
    auto sf_csc = std::make_shared<BaseTCSC>(W_raw.data(), K, N);
    auto sf_csr = std::make_shared<BaseTCSR>(W_raw.data(), K, N);
    auto sf_ccsc = std::make_shared<CompressedCSC>(W_raw.data(), K, N);
    auto sf_icsr = std::make_shared<ICSR>(W_raw.data(), K, N);
    auto sf_icsc = std::make_shared<ICSC>(W_raw.data(), K, N);
    auto sf_blocked = std::make_shared<BlockedTCSC<1024>>(W_raw.data(), K, N);
    auto sf_blocked_interleaved = std::make_shared<BlockedTCSC_interleaved<1024>>(W_raw.data(), K, N);

    // add_function(
    //     [sf_csc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         BaseCSC<float>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "BaseCSC_naive");

    add_function(
        [sf_icsc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            ICSC_base<float>(X_arg, *sf_icsc, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "TCSC_interleaf");

    add_function(
        [sf_csc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            BaseCSC_unr<float, 4>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "BaseCSC_unrolled_4");

    // add_function(
    //     [sf_blocked](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         BlockedCSC<float, 1024>(X_arg, *sf_blocked, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "BlockedCSC_1024");

    add_function(
        [sf_blocked](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            BlockedCSC_unr4<float, 1024>(X_arg, *sf_blocked, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "BlockedCSC_unr_1024");

    add_function(
        [sf_csr](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            BaseCSR<float>(X_arg, *sf_csr, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "BaseCSR_naive");

    add_function(
        [sf_icsr](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            ICSR_base<float>(X_arg, *sf_icsr, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "TCSR_interleaf");

    // add_function(
    //     [sf_csc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         BaseCSC_unr<float, 5>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "BaseCSC_unrolled_5");

    // add_function(
    //     [sf_csc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         BaseCSC_unr_tiled<float, 32, 32, 12>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "BaseCSC_unrolled_tiled_32x32x12");

    // add_function(
    //     [sf_csr](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         BaseCSR_unr<float, 8>(X_arg, *sf_csr, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "BaseCSR_unrolled_8");

    // add_function(
    //     [sf_csr](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         BaseCSR_unr<float, 16>(X_arg, *sf_csr, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "BaseCSR_unrolled_16");

    // add_function(
    //     [sf_ccsc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         CCSC_base<float>(X_arg, *sf_ccsc, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "CompressedCSC_naive");

    // add_function(
    //     [sf_ccsc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         CCSC_unr<float>(X_arg, *sf_ccsc, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "CompressedCSC_unrolled_5");

    // add_function(
    //     [sf_icsr](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         ICSR_base<float>(X_arg, *sf_icsr, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "ICSR_naive");

    // add_function(
    //     [sf_icsc](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         ICSC_base<float>(X_arg, *sf_icsc, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "ICSC_naive");

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
                std::cout << "Test case " << funcNames[i_loop] << " failed!" << std::endl;
            }
        }
    }

    for (i_loop = 0; i_loop < numFuncs; i_loop++)
    {
        perf_val = perf_test(userFuncs[i_loop], M, K, N, nonZero);
        std::cout << "\nRunning: " << "\x1b[31m" << funcNames[i_loop] << "\x1b[0m" << std::endl;
        std::cout << perf_val << " cycles" << std::endl;
#ifdef INSTRUMENTATION_RUN
        std::cout << "Flops: " << getTotalFlops() << std::endl;
        std::cout << "Performance: " << (double)getTotalFlops() / perf_val << " flops/cycle" << std::endl;
        double total_bytes = sizeof(float) * ((double)(M * K + M * N + N)) + (double)getDataStructureSizeInBytes();
        std::cout << "Total Input Size: " << (int)total_bytes << " Bytes" << std::endl;
        std::cout << "Operational Intensity: " << (double)getTotalFlops() / total_bytes << " Flops/Byte" << std::endl;
        std::cout << "Data Structure Size: " << getDataStructureSizeInBytes() << " Bytes" << std::endl;
#endif
        // std::cout << "Performance: " << static_cast<double>(M * N) * (1.0 + static_cast<double>(K) / nonZero) / perf_val << " flops/cycle" << std::endl;
    }

    return 0;
}
