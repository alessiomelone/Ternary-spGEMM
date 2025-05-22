/* 
 *  Compile: 
 *  $ g++ -O3 -DPMU -DNDEBUG -march=native -mtune=native -fstrict-aliasing  cpp_impl/main.cpp cpp_impl/comp.cpp cpp_impl/perf.cpp cpp_impl/data_structures/RSR/naive.cpp cpp_impl/data_structures/RSR/rsr*.cpp cpp_impl/data_structures/RSR/utils.cpp -o cpp_impl/SparseGEMM.out
 *  Run (K must equal N!!):
 *  $ sudo ./cpp_impl/SparseGEMM.out -M 128 -K 2048 -N 2048 -s 2
 */
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <assert.h>
#include "naive.h"
#include "rsr.h"
#include "rsrpp.h"
#include "utils.h"
#include "rsr_driver.h"


using namespace std;

static pair<vector<vector<int>>,vector<vector<int>>> convertSparseToTwoBinaryMatrices(vector<int> W_raw, int K, int N);
static vector<int> calc_indices(const vector<vector<int>>& permutations, const vector<vector<int>>& segments, int k);

// X: MxK
// W: KxN
RSR::RSR(vector<int> W_raw, int K, int N) {
    int val_k = (int) (ceil(log2(K) - log2(log2(K))));  // TODO: Experiment with value of k
    this->k = val_k;

    auto pr = convertSparseToTwoBinaryMatrices(W_raw, K, N);
    auto Bminus = preprocess(pr.first, k);
    auto Bplus = preprocess(pr.second, k);

    auto Bmin_ind2_vec = calc_indices(Bminus.first, Bminus.second, k);
    auto Bplu_ind2_vec = calc_indices(Bplus.first, Bplus.second, k);
    int elements = Bmin_ind2_vec.size();
    this->Bmin_ind2 = (int *) malloc(elements * sizeof(int));
    for (int i = 0; i < elements; i++) {
        this->Bmin_ind2[i] = Bmin_ind2_vec[i];
    }
    
    elements = Bplu_ind2_vec.size();
    this->Bplu_ind2 = (int *) malloc(elements * sizeof(int));
    for (int i = 0; i < elements; i++) {
        this->Bplu_ind2[i] = Bplu_ind2_vec[i];
    }

    this->perm0_size = Bminus.first[0].size();
    this->ssize = Bminus.second[0].size();
    this->perm_size = Bminus.first.size();
    
    int us_buf_size = (int) pow(2, k);
    this->us_buf_size_f = us_buf_size * sizeof(float);
    this->us_buf = (float *) calloc(us_buf_size, sizeof(float));
    this->result = (float *) malloc(perm0_size * sizeof(float));
}

RSR::~RSR() {
    free(us_buf);
    free(result);
    free(Bmin_ind2);
    free(Bplu_ind2);
}

void MMPlusB(float *X_arg, RSR& rsr, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg) {
    int k = rsr.k;
    auto Bmin_ind = rsr.Bmin_ind2;
    auto Bplu_ind = rsr.Bplu_ind2;
    int perm0_size = rsr.perm0_size;
    int ssize = rsr.ssize;
    int perm_size = rsr.perm_size;
    float *us_buf = rsr.us_buf;
    float *result = rsr.result;
    int us_buf_size_f = rsr.us_buf_size_f;
    int init_j0 = pow(2, k);
    int init_j2 = init_j0 >> 1;

    for (int row = 0; row < M_arg; row++) {
        float *v = X_arg + row * K_arg;
        float *res = Y_arg + row * N_arg;

        int indices_plu_idx = 0, indices_min_idx = 0;
        for (int i = 0; i < perm_size; i++) {
            for (int j = 0; j < ssize; j++) {
                float acc = 0;
                int len = Bplu_ind[indices_plu_idx++];
                for (int l = 0; l < len; l++) {
                    acc += v[Bplu_ind[indices_plu_idx++]];
                }
                len = Bmin_ind[indices_min_idx++];
                for (int l = 0; l < len; l++) {
                    acc -= v[Bmin_ind[indices_min_idx++]];
                }
                us_buf[j] = acc;
            }

            int sum;
            int init_max_j = init_j0;
            int init_max_j_2 = init_j2;
            for (int r = k; r > 0; r--) {
                sum = 0;
                for (int j = 1; j < init_max_j; j += 2) {
                    sum += us_buf[j];
                }
                
                res[i * k + r - 1] += sum;
                for (int j = 0; j < init_max_j_2; j++) {
                    us_buf[j] = us_buf[j * 2] + us_buf[j * 2 + 1];
                }
                init_max_j >>= 1;
                init_max_j_2 >>= 1;
            }
        }

        for (int i = 0; i < N_arg; i+=8) {
            res[i] += B_arg[i];
            res[i+1] += B_arg[i+1];
            res[i+2] += B_arg[i+2];
            res[i+3] += B_arg[i+3];
            res[i+4] += B_arg[i+4];
            res[i+5] += B_arg[i+5];
            res[i+6] += B_arg[i+6];
            res[i+7] += B_arg[i+7];
        }
    }
}


static pair<vector<vector<int>>,vector<vector<int>>> convertSparseToTwoBinaryMatrices(vector<int> W_raw, int K, int N) {
    vector<vector<int>> Bminus(K, vector<int>(N));
    vector<vector<int>> Bplus(K, vector<int>(N));
    for (size_t i = 0; i < K; i++) {
        for (size_t j = 0; j < N; j++) {
            int index = i * N + j;
            if (W_raw[index] == 1) {
                Bplus[i][j] = 1;
                Bminus[i][j] = 0;
            }
            else if (W_raw[index] == -1) {
                Bplus[i][j] = 0;
                Bminus[i][j] = 1;
            }
            else {
                Bplus[i][j] = 0;
                Bminus[i][j] = 0;
            }
        }
    }
    return make_pair(Bminus, Bplus);
}

static vector<int> calc_indices(const vector<vector<int>>& permutations, const vector<vector<int>>& segments, int k) {
    int n = permutations[0].size();
    vector<int> res;

    int start;
    int end;
    vector<int> segment;
    vector<int> permutation;
    for (int i = 0; i < permutations.size(); i++) {
        segment = segments[i];
        permutation = permutations[i];

        // Each block
        for (int j = 0; j < segment.size(); j++) {
            start = segment[j];
            if (j < segment.size() - 1) {
                end = segment[j + 1];
            } else {
                end = n;
            }
            // Segmented sum
            int length = (end - start) >= 0 ? (end - start) : 0;
            vector<int> sequence(length);
            res.push_back( length );
            for (int index = start; index < end; index++) {
                res.push_back( permutation[index] );
            }
        }
    }
    return res;
}