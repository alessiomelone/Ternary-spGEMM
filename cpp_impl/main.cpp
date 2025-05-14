#include "perf.h"
#include "common.h"
#include "sparseUtils.h"
#include "comp.h"


// --- End Prototypes ---

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


    //Generate sparse matrix to be converted
    std::vector<int> W_raw = generateSparseMatrix<int>(K, N, nonZero, false);

    //SpraseFormatCSC
    auto sf_csc_data = std::make_shared<SparseFormatCSC>(W_raw.data(), K, N);

    add_function(
        [sf_csc_data](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg) {
            sparseGEMM_csr_base_impl<float>(X_arg, *sf_csc_data, B_arg, Y_arg, M_arg, N_arg, K_arg);
        },
        "sparseGEMM_csr_base"
    );

    add_function(
        [sf_csc_data](float *X_arg, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg) {
            sparseGEMM_csr_unrolled_impl<float, 16>(X_arg, *sf_csc_data, B_arg, Y_arg, M_arg, N_arg, K_arg);
            // Note: You can vary UNROLL_FACTOR here or make it part of the name if you test multiple unroll factors
        },
        "sparseGEMM_csr_unrolled_16"
    );



    if (numFuncs == 0)
    {
        std::cout << std::endl;
        std::cout << "No functions registered - nothing for driver to do" << std::endl;
        // The old message about register_funcs() in comp.cpp is no longer relevant here.
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

