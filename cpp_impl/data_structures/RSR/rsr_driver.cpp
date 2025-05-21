
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "naive.h"
#include "rsr.h"
#include "rsrpp.h"
#include "utils.h"
#include "rsr_driver.h"

using namespace std;


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
    auto pre_Bminus = preprocess(pr.first, k);
    auto pre_Bplus = preprocess(pr.second, k);
    this->bin_k = generateBinaryMatrix2(k);
    this->k = k;
    this->n = pre_Bminus.first[0].size();
    for (int i = 0; i < pre_Bminus.second.size(); i++) {
        this->segment_sizes_Bminus.push_back(pre_Bminus.second[i].size());
    }
    for (int i = 0; i < pre_Bplus.second.size(); i++) {
        this->segment_sizes_Bplus.push_back(pre_Bplus.second[i].size());
    }
    this->permutations_size = pre_Bminus.first.size();
    
    auto& permutations = pre_Bminus.first;
    auto& segments = pre_Bminus.second;
    // vector<vector<int>> start_end_idx(permutations.size(), vector<int>(2 * segments.size()));
    for (int i = 0; i < permutations.size(); i++) {
        auto& segment = segments[i];
        auto& permutation = permutations[i];
        int start, end;
        // Each block
        for (int j = 0; j < segment.size(); j++) {
            start = segment[j];
            if (j < segment.size() - 1) {
                end = segment[j + 1];
            } else {
                end = n;
            }
            // Segmented sum
            int id = end - start;
            vector<int> local_indices(id > 0 ? id : 0);
            for (int index = start, k = 0; index < end; index++, k++) {
                local_indices[k] = permutation[index];
            }
            this->indices_Bminus.push_back(local_indices);
        }
    }
    permutations = pre_Bplus.first;
    segments = pre_Bplus.second;
    // vector<vector<int>> start_end_idx(permutations.size(), vector<int>(2 * segments.size()));
    for (int i = 0; i < permutations.size(); i++) {
        auto& segment = segments[i];
        auto& permutation = permutations[i];
        int start, end;
        // Each block
        for (int j = 0; j < segment.size(); j++) {
            start = segment[j];
            if (j < segment.size() - 1) {
                end = segment[j + 1];
            } else {
                end = n;
            }
            // Segmented sum
            int id = end - start;
            vector<int> local_indices(id > 0 ? id : 0);
            for (int index = start, k = 0; index < end; index++, k++) {
                local_indices[k] = permutation[index];
            }
            this->indices_Bplus.push_back(local_indices);
        }
    }

    this->powK = pow(2, k);
}

// static void rsr_inference2(float *v, float *res, bool isMinus, int permutations_size, 
//     const vector<int>& segment_sizes, const vector<vector<int>>& bin_k, int k, const vector<vector<int>>& idx, vector<vector<float>>& us) {
//    // segmented sums
//     int index = 0;
//     for (int i = 0; i < permutations_size; i++) {
//         int seg_size = segment_sizes[i];
//         cout << seg_size << " , k=" << k << endl;
//         for (int j = 0; j < seg_size; j++) {
//             auto& local_indices = idx[index++];
//             for (int l = 0; l < local_indices.size(); l++) {
//                 us[i][j] += v[local_indices[l]];
//             }
//             if (! isMinus) {
//                 for (int l=0; l < k; l++) {
//                     res[i*k + l] += us[i][j] * bin_k[j][l];
//                 }
//             } else {
//                 for (int l=0; l < k; l++) {
//                     res[i*k + l] -= us[i][j] * bin_k[j][l];
//                 }
//             }
//         }
//     }
// }


void MMPlusB(float *X_arg, RSR& rsr, float *B_arg, float *Y_arg, int M_arg, int N_arg, int K_arg) {
    int n = rsr.n;
    int perm_size = rsr.permutations_size;
    auto& bin_k = rsr.bin_k;
    auto& seg_size = rsr.segment_sizes_Bplus[0];
    int k = rsr.k;
    for (int row = 0; row < M_arg; row++) {
        float *v = X_arg + row * K_arg;
        float *res = Y_arg + row * N_arg;
        auto us = vector<vector<float>>(perm_size, vector<float>(rsr.powK));
        auto uss = vector<vector<float>>(perm_size, vector<float>(rsr.powK));
        int index = 0, index1 = 0;
        for (int i = 0; i < perm_size; i++) {
            for (int j = 0; j < seg_size; j++) {
                auto& local_indices = rsr.indices_Bplus[index++];
                auto& local_indices_min = rsr.indices_Bminus[index1++];

                for (int l = 0; l < local_indices.size(); l++) {
                    us[i][j] += v[local_indices[l]];
                }
                for (int l = 0; l < local_indices_min.size(); l++) {
                    uss[i][j] += v[local_indices_min[l]];
                }
                for (int l=0; l < k; l++) {
                    res[i*k + l] += (us[i][j] - uss[i][j]) * bin_k[j][l];
                }
            }
        }

        for (int i = 0; i < n; i++) {
            res[i] += B_arg[i];
        }
    }
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
