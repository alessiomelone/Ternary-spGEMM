<<<<<<< HEAD
/* 
 *  Compile: 
 *  $ g++ -O3 -DPMU -DNDEBUG -march=native -mtune=native -fstrict-aliasing  cpp_impl/main.cpp cpp_impl/comp.cpp cpp_impl/perf.cpp cpp_impl/data_structures/RSR/naive.cpp cpp_impl/data_structures/RSR/rsr*.cpp cpp_impl/data_structures/RSR/utils.cpp -o cpp_impl/SparseGEMM.out
 *  Run (K must equal N!!):
 *  $ sudo ./cpp_impl/SparseGEMM.out -M 128 -K 2048 -N 2048 -s 2
 */
=======

>>>>>>> e9a5eca (Add RSR data structure)
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
<<<<<<< HEAD
#include <assert.h>
=======
>>>>>>> e9a5eca (Add RSR data structure)
#include "naive.h"
#include "rsr.h"
#include "rsrpp.h"
#include "utils.h"
#include "rsr_driver.h"

<<<<<<< HEAD

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

=======
using namespace std;

static void rsr_inference_matrix(float *matrix, float *Y, float *B, int M, int K, int N, const RSR& rsr);
>>>>>>> e9a5eca (Add RSR data structure)

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

<<<<<<< HEAD
static vector<int> calc_indices(const vector<vector<int>>& permutations, const vector<vector<int>>& segments, int k) {
    int n = permutations[0].size();
    vector<int> res;
=======
static vector<vector<int>> generateBinaryMatrix2(int k) {
    int rows = pow(2, k);  // 2^k rows
    vector<vector<int>> matrix(rows, vector<int>(k, 0));  // Initialize matrix with 0s

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < k; ++j) {
            // Generate the binary value for each position
            matrix[i][k - j - 1] = (i >> j) & 1;  // Extract the j-th bit from i
        }
    }

    return matrix;
}

RSR::RSR(vector<int> W_raw, int K, int N) {
    int n = K; // > N ? K : N; // TODO : Experiment with that
    int k = static_cast<int>(ceil(log2(n) - log2(log2(n))));
    auto pr = convertSparseToTwoBinaryMatrices(W_raw, K, N);
    this->pre_Bminus = preprocess(pr.first, k);
    this->pre_Bplus = preprocess(pr.second, k);
    this->bin_k = generateBinaryMatrix(k);
    this->k = k;
}
// RSR *ConvertToRSR(vector<int> W_raw, int K, int N) {
//     RSR *rsr = new RSR();
//     int n = K; // > N ? K : N; // TODO : Experiment with that
//     int k = static_cast<int>(ceil(log2(n) - log2(log2(n))));
//     auto pr = convertSparseToTwoBinaryMatrices(W_raw, K, N);
//     rsr->pre_Bminus = preprocess(pr.first, k);
//     rsr->pre_Bplus = preprocess(pr.second, k);
//     rsr->bin_k = generateBinaryMatrix(k);
//     rsr->k = k;
//     return rsr;
// }



static void rsr_inference2(float *v, float *res, bool isMinus, const vector<vector<int>>& permutations, const vector<vector<int>>& segments, const vector<vector<int>>& bin_k, int k) {
    int n = permutations[0].size();

    // segmented sums
    vector<vector<float>> us(permutations.size(), vector<float>(pow(2, k)));
>>>>>>> e9a5eca (Add RSR data structure)

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
<<<<<<< HEAD
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
=======
            for (int index = start; index < end; index++) {
                us[i][j] += v[permutation[index]];
            }           
        }
    }

    // vector<float> result(n);


    // Block product to Bin_k
    // TODO: change from here for RSR++
    if (! isMinus) {
        vector<float> partial_result;
        for (int i = 0; i < us.size(); i++) {
            partial_result = vectorMatrixMultiply2(us[i], bin_k);
            for (int j = 0; j < k; j++) {
                res[i * k + j] += partial_result[j];
            }
        }
    } else {
        vector<float> partial_result;
        for (int i = 0; i < us.size(); i++) {
            partial_result = vectorMatrixMultiply2(us[i], bin_k);
            for (int j = 0; j < k; j++) {
                res[i * k + j] -= partial_result[j];
            }
        }
    }

    // return result;
}

static void rsr_inference_matrix2(float *matrix, float *Y, float *B, int M, int K, int N, const RSR& rsr) {
    for (int row = 0; row < M; row++) {
        float *v = matrix + row * K;
        int n = rsr.pre_Bminus.first[0].size();
        float *res = Y + row * N;
        rsr_inference2(v, res, 1, rsr.pre_Bminus.first, rsr.pre_Bminus.second, rsr.bin_k, rsr.k);
        rsr_inference2(v, res, 0, rsr.pre_Bplus.first, rsr.pre_Bplus.second, rsr.bin_k, rsr.k);
        for (int i = 0; i < n; i++) {
            res[i] += B[i];
        }
    }
}

void MMPlusB(float *X_arg, const RSR& rsr, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg) {
    rsr_inference_matrix2(X_arg, Y_arg, B_arg, M_arg, K_arg, N_arg, rsr);
    cout << "exititng..." << endl;
}


#if 0
int main() {
    for (int i = 2; i <= 12; i++) {
        int n = pow(2, i);
        cout << "n = " << n << endl;
        int k = static_cast<int>(ceil(log2(n) - log2(log2(n))));

        vector<vector<int>> mat = generateBinaryRandomMatrix(n);
        vector<int> v = generateRandomVector(n);
        vector<vector<int>> bin_k = generateBinaryMatrix(k);
        
        cout << "preprocessing..." << endl;
        auto per_segs = preprocess(mat, k);

        cout << "inference..." << endl;

        vector<int> actual_1 = rsr_pp_inference(v, per_segs.first, per_segs.second, k);
        vector<int> actual_2 = rsr_inference(v, per_segs.first, per_segs.second, bin_k, k);
        vector<int> expected = vectorMatrixMultiply(v, mat);

        for(int i = 0; i < n; i++) {
            if (actual_1[i] != expected[i]) {
                cout << "ERROR" << endl;
                return 1;
            }
            if (actual_2[i] != expected[i]) {
                cout << "ERROR" << endl;
                return 1;
            }
        }
    }
    cout << "PASSED" << endl;

    return 0;
}


// X is MxK
// W is KxN
static void rsr_inference_matrix(float *matrix, float *Y, float *B, int M, int K, int N, const RSR& rsr) {
    //const vector<vector<int>>& permutations, const vector<vector<int>>& segments, vector<vector<int>> bin_k, int k) {
    auto& permutations_Bminus = rsr.pre_Bminus.first;
    auto& permutations_Bplus = rsr.pre_Bplus.first;
    auto& segmenents_Bminus = rsr.pre_Bminus.second;
    auto& segmenents_Bplus = rsr.pre_Bplus.second;
    auto& bin_k = rsr.bin_k;
    int k = rsr.k;

    int n = permutations_Bminus[0].size();
    // int n_Bplus  = permutations_Bplus[0].size();

    // segmented sums
    int powk = pow(2, k);
    vector<vector<int>> us(permutations_Bminus.size(), vector<int>(powk));

    int start_Bminus, start_Bplus, end_Bminus, end_Bplus;
    vector<int> segmenent_Bminus, permutation_Bminus, segmenent_Bplus, permutation_Bplus;

    for (int row = 0; row < M; row++) {
        float *v = matrix + row * K;
    // TODO : Verify that both outer loops have the same size (hmm)
        for (int i = 0; i < permutations_Bminus.size(); i++) {
            segmenent_Bminus = segmenents_Bminus[i];
            permutation_Bminus = permutations_Bminus[i];
            segmenent_Bplus = segmenents_Bplus[i];
            permutation_Bplus = permutations_Bplus[i];

            // Each block
            for (int j = 0; j < segmenent_Bminus.size(); j++) {
                start_Bminus = segmenent_Bminus[j];
                start_Bplus = segmenent_Bplus[j];
                if (j < segmenent_Bminus.size() - 1) {
                    end_Bminus = segmenent_Bminus[j + 1];
                    end_Bplus = segmenent_Bplus[j+1];
                } else {
                    end_Bminus = n; // segmenent_Bminus[j + 1];
                    end_Bplus = n; // segmenent_Bplus[j+1];
                }
                // Segmented sum
                for (int index = start_Bminus; index < end_Bminus; index++) {
                    us[i][j] -= v[permutation_Bminus[index]];
                }
                for (int index = start_Bplus; index < end_Bplus; index++) {
                    us[i][j] += v[permutation_Bplus[index]];
                }           
            }
        }

        // vector<int> result(n);
        float *result = Y + row * N;

        // Block product to Bin_k
        // TODO: change from here for RSR++
        vector<int> partial_result;
        for (int i = 0; i < us.size(); i++) {
            partial_result = vectorMatrixMultiply(us[i], bin_k);
            for (int j = 0; j < k; j++) {
                result[i * k + j] = partial_result[j] + B[j];  // TODO: Verify B here
            }
        }
    }
}
#endif
>>>>>>> e9a5eca (Add RSR data structure)
