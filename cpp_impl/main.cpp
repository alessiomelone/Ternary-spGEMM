#include "perf.h"
#include "common.h"
#include "sparseUtils.h"
#include "comp.h"

#define BLOCK_SIZE_IBTCSC 512
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
    auto sf_interleaved_baraq = std::make_shared<InterleavedTCSC_baraq>(W_raw.data(), K, N);

    auto sf_interleaved = std::make_shared<InterleavedTCSC>(W_raw.data(), K, N);
    auto sf_interleaved_padding = std::make_shared<InterleavedTCSCPadding>(W_raw.data(), K, N);
    auto sf_interleaved_blocked_harry = std::make_shared<BlockedTCSC_interleaved<BLOCK_SIZE_IBTCSC>>(W_raw.data(), K, N);
    auto sf_interleaved_blocked_unrolled_harry = std::make_shared<BlockedTCSC_interleaved<BLOCK_SIZE_IBTCSC>>(W_raw.data(), K, N, UNROLL_FACTOR_IBTCSC);    
    
    add_function(
        [sf_csc](double *X_arg, double *B_arg, double *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            BaseCSC<double>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "BaseCSC_naive");
    
        add_function(
        [sf_interleaved_baraq](double *X_arg, double *B_arg, double *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            InterleavedTCSC_baraq_comp<double>(X_arg, *sf_interleaved_baraq, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "InterleavedTCSC_baraq_comp");

        add_function(
        [sf_interleaved_baraq](double *X_arg, double *B_arg, double *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            InterleavedTCSC_baraq_comp_unr<double, 16>(X_arg, *sf_interleaved_baraq, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "InterleavedTCSC_baraq_comp_unr_16");

    // add_function(
    //     [sf_csc](double *X_arg, double *B_arg, double *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         TCSC_inter<double>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "TCSC_interleaf");

    // add_function(
    //     [sf_interleaved](double *X_arg, double *B_arg, double *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         TCSC_interleaved_ds<double>(X_arg, *sf_interleaved, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "TCSC_interleaf with DS 1/-1");

    // add_function(
    //     [sf_interleaved_padding](double *X_arg, double *B_arg, double *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         TCSC_interleaved_padding<double>(X_arg, *sf_interleaved_padding, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "TCSC_interleaf with padding");

    add_function(
        [sf_csc](double *X_arg, double *B_arg, double *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            BaseCSC_unr<double, 16>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "BaseCSC_unrolled_16");

    // add_function(
    //     [sf_csc](double *X_arg, double *B_arg, double *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         TCSC_inter_unr<double, 16>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "TCSC_interleaf_unrolled_16");

        std::string s = "BlockedTCSC_interleaved_base_Block_Size:" +  std::to_string(BLOCK_SIZE_IBTCSC) ;
            add_function(
        [sf_interleaved_blocked_harry](double *X_arg, double *B_arg, double *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            BlockedTCSC_interleaved_base<double, BLOCK_SIZE_IBTCSC>(X_arg, *sf_interleaved_blocked_harry, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
         s.data());
    
        s = "BlockedTCSC_interleaved_unrolled_Block_Size:" +  std::to_string(BLOCK_SIZE_IBTCSC)  + "_Unroll_Factor:" + std::to_string(UNROLL_FACTOR_IBTCSC);
            add_function(
        [sf_interleaved_blocked_unrolled_harry](double *X_arg, double *B_arg, double *Y_arg, int M_arg, int N_arg, int K_arg)
        {
            BlockedTCSC_interleaved_unr<double, BLOCK_SIZE_IBTCSC, UNROLL_FACTOR_IBTCSC>(X_arg, *sf_interleaved_blocked_unrolled_harry, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
         s.data());

    // add_function(
    //     [sf_csc](double *X_arg, double *B_arg, double *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         BaseCSC_unr_tiled<double, 12, 12, 12>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "BaseCSC_unrolled_tiled_12x12x12");

    // add_function(
    //     [sf_csc](double *X_arg, double *B_arg, double *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         TCSC_inter_unr_tiled<double, 12, 12, 12>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "TCSC_inter_unr_tiled_12x12x12");

    // add_function(
    //     [sf_blocked](double *X_arg, double *B_arg, double *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         BlockedCSC<double, 1024>(X_arg, *sf_blocked, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "BlockedCSC_1024");

    // add_function(
    //     [sf_blocked](double *X_arg, double *B_arg, double *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         BlockedCSC_unr4<double, 1024>(X_arg, *sf_blocked, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "BlockedCSC_unr_1024");

    // add_function(
    //     [sf_csr](double *X_arg, double *B_arg, double *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         BaseCSR<double>(X_arg, *sf_csr, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "BaseCSR_naive");

    // add_function(
    //     [sf_csr](double *X_arg, double *B_arg, double *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         TCSR_inter<double>(X_arg, *sf_csr, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "TCSR_interleaf");

    // add_function(
    //     [sf_csc](double *X_arg, double *B_arg, double *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         BaseCSC_unr<double, 5>(X_arg, *sf_csc, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "BaseCSC_unrolled_5");

    // add_function(
    //     [sf_csr](double *X_arg, double *B_arg, double *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         BaseCSR_unr<double, 8>(X_arg, *sf_csr, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "BaseCSR_unrolled_8");

    // add_function(
    //     [sf_csr](double *X_arg, double *B_arg, double *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         BaseCSR_unr<double, 16>(X_arg, *sf_csr, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "BaseCSR_unrolled_16");

    // add_function(
    //     [sf_ccsc](double *X_arg, double *B_arg, double *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         CCSC_base<double>(X_arg, *sf_ccsc, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "CompressedCSC_naive");

    // add_function(
    //     [sf_ccsc](double *X_arg, double *B_arg, double *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         CCSC_unr<double>(X_arg, *sf_ccsc, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "CompressedCSC_unrolled_5");

    // add_function(
    //     [sf_icsr](double *X_arg, double *B_arg, double *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         ICSR_base<double>(X_arg, *sf_icsr, B_arg, Y_arg, M_arg, N_arg, K_arg);
    //     },
    //     "ICSR_naive");

    // add_function(
    //     [sf_icsc](double *X_arg, double *B_arg, double *Y_arg, int M_arg, int N_arg, int K_arg)
    //     {
    //         ICSC_base<double>(X_arg, *sf_icsc, B_arg, Y_arg, M_arg, N_arg, K_arg);
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

    std::vector<double> X_main = initX<double>(M * K, 512);
    std::vector<double> W_FP32_main(W_raw.begin(), W_raw.end());
    std::vector<double> B_main(N, 2);
    std::vector<double> Y_main(M * N, 0);
    std::vector<double> refY_main(M * N, 0);

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
                std::cout << "Test case "  << "\x1b[31m" << funcNames[i_loop] << " failed!" << "\x1b[0m" <<  std::endl;
                std::cout << "\n\n  Please fix the failing fn or comment out the invocaton from main.cpp.\n\nExiting...\n\n" << std::endl;
                exit(1);
            }
        }
    }

    double base_cycles = 0;
    for (i_loop = 0; i_loop < numFuncs; i_loop++)
    {
        perf_val = perf_test(userFuncs[i_loop], M, K, N, nonZero);
        std::cout << "\nRunning: " << "\x1b[31m" << funcNames[i_loop] << "\x1b[0m" << std::endl;
        std::cout << perf_val << " cycles" << std::endl;
        if (funcNames[i_loop] == "BaseCSC_naive") {
            base_cycles = perf_val;
        }
        std::cout << "Speedup is: " << "\x1b[32m" << base_cycles / perf_val << "\x1b[0m" << std::endl;
#ifdef INSTRUMENTATION_RUN
        std::cout << "Flops: " << getTotalFlops() << std::endl;
        std::cout << "Performance: " << (double)getTotalFlops() / perf_val << " flops/cycle" << std::endl;
        double total_bytes = sizeof(double) * ((double)(M * K + M * N + N)) + (double)getDataStructureSizeInBytes();
        std::cout << "Total Input Size: " << (int)total_bytes << " Bytes" << std::endl;
        std::cout << "Operational Intensity: " << (double)getTotalFlops() / total_bytes << " Flops/Byte" << std::endl;
        std::cout << "Data Structure Size: " << getDataStructureSizeInBytes() << " Bytes" << std::endl;
#endif
        // std::cout << "Performance: " << static_cast<double>(M * N) * (1.0 + static_cast<double>(K) / nonZero) / perf_val << " flops/cycle" << std::endl;
    }

    return 0;
}
